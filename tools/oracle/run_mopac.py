# SPDX-License-Identifier: GPL-3.0-or-later
"""Run MOPAC v23.2.5 on a geometry and parse the .aux/.out results.

Set ``MOPAC_EXE`` to an official OpenMOPAC executable. The binary is not
redistributed with xndo-rs.

Usage:
    python run_mopac.py <xyz-file> [--method MNDO]
        [--charge Q] [--mult M] [--mode 1scf|gradient|force] [--keep]

Prints a JSON dict with heat_of_formation_kcal, energies, charges, dipole,
gradients (kcal/mol/A) and frequencies (cm-1) when available.

"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
MOPAC_EXE = os.environ.get(
    "MOPAC_EXE",
    os.path.join(HERE, "mopac", "mopac-23.2.5-win", "bin", "mopac.exe"),
)


def read_xyz(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    nat = int(lines[0].split()[0])
    atoms = []
    for line in lines[2 : 2 + nat]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms


# MOPAC's SCF convergence, tightened past what `PRECISE` gives.
#
# `PRECISE` alone leaves MOPAC's own orbital energies uncertain at about 4e-4 eV,
# which was large enough to be mistaken for a model difference: CS2 under INDO
# disagreed with xndo-rs by 3.4e-4 eV per orbital, and MOPAC's *own* eigenvalues
# move by 3.8e-4 eV when the SCF is tightened. At `RELSCF=0.0001` the same
# comparison closes to 1.05e-5 eV, which is two units in the last printed digit.
#
# Smaller values gain nothing -- 1e-6 gives the identical answer -- so this is
# the point where MOPAC's SCF stops being the limit and its print format starts.
SCF_TIGHTENING = "RELSCF=0.0001"


def write_mop(path, atoms, method="MNDO", charge=0, mult=1, mode="1scf", extra_keywords=()):
    keywords = [method, "PRECISE", SCF_TIGHTENING, f"CHARGE={charge}", "AUX(PRECISION=9)"]
    if mode == "1scf":
        keywords.append("1SCF")
    elif mode == "gradient":
        keywords += ["1SCF", "GRADIENTS"]
    elif mode == "force":
        keywords.append("FORCE")
    elif mode != "optimize":
        raise ValueError("mode must be 1scf, gradient, force, or optimize")
    if mult > 1:
        keywords.append({2: "DOUBLET", 3: "TRIPLET", 4: "QUARTET", 5: "QUINTET"}.get(mult, f"MS={(mult-1)/2}"))
        keywords.append("UHF")
    keywords.extend(extra_keywords)
    with open(path, "w") as fh:
        fh.write(" ".join(keywords) + "\n")
        fh.write("xndo-rs independent oracle reference\n\n")
        for sym, x, y, z in atoms:
            fh.write(f"{sym:3s} {x:15.8f} 1 {y:15.8f} 1 {z:15.8f} 1\n")


def parse_aux(path):
    """Parse MOPAC .aux into {key: float|list-of-floats}."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    errors = re.findall(r"^ ERROR.*$", text, re.M)
    if errors:
        raise RuntimeError("MOPAC AUX reported: " + " | ".join(errors))
    out = {}
    # Scalar entries: KEY:UNITS=value  or KEY=value
    for m in re.finditer(r"^ ([A-Z_.0-9]+)(?::([A-Z0-9/().^-]*))?=(.+)$", text, re.M):
        key, _units, val = m.group(1), m.group(2), m.group(3).strip()
        try:
            out[key] = float(val.replace("D", "E"))
        except ValueError:
            out[key] = val
    # Vector entries: KEY:UNITS[n]= v1 v2 ... (block until next KEY)
    for m in re.finditer(
        r"^ ([A-Z_.0-9]+)(?::([A-Z0-9/().^-]*))?\[(\d+)\]=\s*\n?((?:[^\n=]*\n)*?)(?=^ [A-Z#])",
        text,
        re.M,
    ):
        key, _units, _n, block = m.groups()
        vals = []
        ok = True
        for tok in block.split():
            try:
                vals.append(float(tok.replace("D", "E")))
            except ValueError:
                ok = False
                break
        if ok and vals:
            out[key] = vals
    return out


def parse_out_eigenvalues(path):
    """Read the EIGENVALUES block from the printed output.

    The AUX block is not reliable for the INDO path: on BF3 it declares
    ``EIGENVALUES[014]`` and lists 14 values where the printed output carries the
    full 16, dropping the lowest orbital and one of a degenerate pair. Benzene
    and CCl4 lose 10 and 6. The printed block is complete, so it is what the
    eigenvalue comparison uses. See ORACLE_NOTES item 17.

    Returns None when the output has no EIGENVALUES block, so the caller can
    fall back to AUX for the methods where AUX is complete.
    """
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            lines = handle.read().splitlines()
    except OSError:
        return None
    values = []
    seen = False
    for index, line in enumerate(lines):
        if "EIGENVALUES" in line and "=" not in line:
            seen = True
            for follow in lines[index + 1:]:
                if not follow.strip():
                    if values:
                        break
                    continue
                fields = follow.split()
                try:
                    values.extend(float(f) for f in fields)
                except ValueError:
                    break
            break
    if not seen:
        return None
    return values or None


def parse_out_ci_count(path):
    """How many spin-adapted configurations MOPAC actually put in the CI.

    This is not `occupied * virtual + 1`. MOPAC's INDO CI filters the
    configuration list by spatial symmetry, and in a high-symmetry point group
    it can throw nearly all of them away: CCl4 in Td, asked for a 5x4 window,
    reports `CI excitations= 1` and diagonalises nothing at all.

    A CIS in this crate has no symmetry, so a symmetry-filtered CI is a
    different calculation wearing the same name. Comparing against one would not
    be a loose comparison, it would be a meaningless one -- hence the check.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "CI excitations=" in line:
                tail = line.split("CI excitations=", 1)[1]
                digits = ""
                for ch in tail.strip():
                    if ch.isdigit():
                        digits += ch
                    else:
                        break
                if digits:
                    return int(digits)
    return None


def parse_out_ci_states(path):
    """Read the INDO CI transition table from the printed output.

    MOPAC's INDO CI writes **two** tables, and they are easy to confuse. The
    first is headed

        sym   eV   cm**-1 -dets- dipole oscilator X FRAG ....Excitations named

    and lists the *spin-adapted configurations* -- the diagonal of the CI matrix
    before it is diagonalised, each labelled by the single excitation it came
    from. The second is headed

        CI trans.  energy frequency wavelength oscillator--- polarization---

    and lists the actual CI eigenstates.

    Only the second is comparable with a CIS calculation. ORACLE_NOTES item 12
    recorded a "5->7 root 1.0 eV low" and "oscillator strengths 20-35% low"
    which were entirely this confusion: for formaldehyde the 5->7 *configuration*
    sits at 9.887 eV with f = 0.547, while the CI root it dominates is at
    8.883 eV with f = 0.018.

    Returns a list of dicts, one per excited root (the ground state is not a
    root and is not in this table), or None if the job had no CI.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    start = text.find("CI trans.")
    if start < 0:
        return None
    # The table ends at the polarizability summary, or at the CI-contribution
    # listing if polarizabilities were not requested.
    end = len(text)
    for marker in ("Polarizability", "Major CI contributions"):
        at = text.find(marker, start)
        if at > 0:
            end = min(end, at)
    states = []
    for line in text[start:end].splitlines():
        fields = line.split()
        # Two shapes: with a polarization vector (bright) and without (dark).
        #   idx  eV  cm-1  nm  f  px py pz  |mu|  mux muy muz     -> 12 fields
        #   idx  eV  cm-1  nm  f              |mu|  mux muy muz   ->  9 fields
        if len(fields) not in (9, 12):
            continue
        try:
            values = [float(v.rstrip(".")) for v in fields]
        except ValueError:
            continue
        if values[0] != int(values[0]):
            continue
        polarization = values[5:8] if len(fields) == 12 else None
        states.append({
            "index": int(values[0]),
            "energy_ev": values[1],
            "energy_cm1": values[2],
            "wavelength_nm": values[3],
            "oscillator_strength": values[4],
            "polarization": polarization,
            "dipole_magnitude_debye": values[-4],
            "dipole_debye": values[-3:],
        })
    return states or None


def parse_out_ci_window(path):
    """The CI window MOPAC actually used, as (first_occ, last_occ, first_vir,
    last_vir) in its own 1-based MO numbering.

    MOPAC echoes the window it chose rather than the one that was asked for, so
    reading it back is what makes a window comparison honest.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "SINGLE excitations FROM orbs" in line:
                fields = line.split()
                # SINGLE excitations FROM orbs   1 to   6 INTO orbs   7 to  10
                return (int(fields[4]), int(fields[6]), int(fields[9]), int(fields[11]))
    return None


def parse_cartesian_hessian(path, aux):
    """Return the MOPAC FORCE Cartesian Hessian in eV/Bohr^2.

    MOPAC writes the lower triangle after a comment line, so the generic AUX
    vector parser intentionally does not consume it. The stored matrix is
    mass-weighted and expressed in millidynes/Angstrom.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    match = re.search(r"^ HESSIAN_MATRIX:[^\n=]+\[(\d+)\]=\s*$", text, re.M)
    if match is None:
        return None
    expected = int(match.group(1))
    values = []
    for line in text[match.end() :].splitlines():
        if line.startswith(" #") or not line.strip():
            continue
        if re.match(r"^ [A-Z][A-Z_.0-9]*(?::[^=]*)?(?:\[\d+\])?=", line):
            break
        for token in line.split():
            values.append(float(token.replace("D", "E")))
    if len(values) != expected:
        raise ValueError(f"expected {expected} packed Hessian values, found {len(values)}")

    masses = aux.get("ISOTOPIC_MASSES")
    if not isinstance(masses, list) or not masses:
        raise ValueError("MOPAC FORCE output did not contain ISOTOPIC_MASSES")
    ndof = 3 * len(masses)
    if expected != ndof * (ndof + 1) // 2:
        raise ValueError("packed Hessian size does not match the atom count")

    # 1 millidyne/Angstrom = 100 N/m. Convert to eV/Angstrom^2, then to
    # eV/Bohr^2. Undo MOPAC's sqrt(m_i*m_j) mass weighting at the same time.
    bohr_to_angstrom = 0.529177210903
    mdyn_per_angstrom_to_ev_per_bohr2 = (100.0 / 16.02176634) * bohr_to_angstrom**2
    dof_masses = [mass for mass in masses for _ in range(3)]
    matrix = [[0.0] * ndof for _ in range(ndof)]
    k = 0
    for i in range(ndof):
        for j in range(i + 1):
            value = (
                values[k]
                * (dof_masses[i] * dof_masses[j]) ** 0.5
                * mdyn_per_angstrom_to_ev_per_bohr2
            )
            matrix[i][j] = value
            matrix[j][i] = value
            k += 1
    return matrix


def run(
    xyz,
    method="MNDO",
    charge=0,
    mult=1,
    mode="1scf",
    keep=False,
    workdir=None,
    extra_keywords=(),
    displacements=(),
):
    atoms = read_xyz(xyz)
    for atom_index, axis, delta in displacements:
        if not 0 <= atom_index < len(atoms):
            raise ValueError(f"atom index {atom_index + 1} is outside the XYZ geometry")
        if axis not in (0, 1, 2):
            raise ValueError(f"axis {axis} is invalid")
        sym, x, y, z = atoms[atom_index]
        coord = [x, y, z]
        coord[axis] += delta
        atoms[atom_index] = (sym, *coord)
    tmp = workdir or tempfile.mkdtemp(prefix="xndors_oracle_")
    os.makedirs(tmp, exist_ok=True)
    name = os.path.splitext(os.path.basename(xyz))[0]
    mop = os.path.join(tmp, name + ".mop")
    write_mop(mop, atoms, method, charge, mult, mode, extra_keywords)
    if not os.path.isfile(MOPAC_EXE):
        # Say which binary is missing. Left to subprocess this surfaces as a bare
        # localised OS error ("file not found"), which names neither MOPAC nor the
        # variable that would fix it, and which the caller then prints once per
        # molecule.
        raise FileNotFoundError(
            "MOPAC_EXE does not point at an OpenMOPAC executable; set it to one "
            "(see tests/data/ORACLE_NOTES.md section (a)). OpenMOPAC is not bundled."
        )
    completed = subprocess.run(
        [MOPAC_EXE, mop], check=True, capture_output=True, text=True, cwd=tmp
    )
    aux_path = os.path.join(tmp, name + ".aux")
    out_path = os.path.join(tmp, name + ".out")
    aux = parse_aux(aux_path)
    result = {
        "method": method,
        "mode": mode,
        "mopac_version": aux.get("MOPAC_VERSION"),
        "heat_of_formation_kcal": aux.get("HEAT_OF_FORMATION"),
        "energy_electronic_ev": aux.get("ENERGY_ELECTRONIC"),
        "energy_nuclear_ev": aux.get("ENERGY_NUCLEAR"),
        "total_energy_ev": aux.get("TOTAL_ENERGY"),
        "dipole_debye": aux.get("DIP_VEC"),
        "charges": aux.get("ATOM_CHARGES"),
        "gradients_kcal_mol_ang": aux.get("GRADIENTS"),
        "frequencies_cm": aux.get("VIB._FREQ"),
        "hessian_ev_per_bohr2": parse_cartesian_hessian(aux_path, aux),
        "spin_sz": aux.get("SPIN_COMPONENT"),
        # AUX(PRECISION=9) carries far more digits than the printed output's
        # 3-4 decimals, so frontier orbital energies are read from here.
        "ionization_potential_ev": aux.get("IONIZATION_POTENTIAL"),
        "eigenvalues_ev": parse_out_eigenvalues(out_path) or aux.get("EIGENVALUES"),
        "mo_occupancies": aux.get("MOLECULAR_ORBITAL_OCCUPANCIES"),
        "n_electrons": aux.get("NUM_ELECTRONS"),
        "aux_keys": sorted(aux.keys()),
        # Present only for INDO CI jobs. Read from the CI transition table, not
        # from the spin-adapted configuration table above it -- see
        # `parse_out_ci_states`.
        "ci_states": parse_out_ci_states(out_path),
        "ci_window": parse_out_ci_window(out_path),
        "ci_count": parse_out_ci_count(out_path),
    }
    # The occupied count is what separates HOMO from LUMO in the eigenvalue list,
    # so it has to be counted on the *same* list the eigenvalues came from.
    #
    # The AUX occupancy vector is truncated exactly where the AUX eigenvalue
    # vector is (ORACLE_NOTES item 17): on BF3 it lists ten 2.0s where twelve
    # orbitals are occupied, because the two it dropped were occupied ones. Using
    # it against the complete printed eigenvalue list puts HOMO and LUMO two
    # orbitals too low. So the occupancies are used only when their length
    # matches the eigenvalue list; otherwise the electron count decides, which is
    # exact for a closed shell and gives the alpha count for an open one.
    occupancies = result["mo_occupancies"]
    eigenvalues = result["eigenvalues_ev"]
    usable = (
        isinstance(occupancies, list)
        and occupancies
        and isinstance(eigenvalues, list)
        and len(occupancies) == len(eigenvalues)
    )
    if usable:
        result["n_occupied"] = sum(1 for o in occupancies if o and float(o) > 0.5)
    elif isinstance(result["n_electrons"], (int, float)):
        electrons = int(result["n_electrons"])
        result["n_occupied"] = (electrons + (mult - 1)) // 2
    else:
        result["n_occupied"] = None
    if result["heat_of_formation_kcal"] is None:
        stderr = completed.stderr.strip()
        raise RuntimeError(f"MOPAC produced no heat of formation{': ' + stderr if stderr else ''}")
    if not keep and workdir is None:
        shutil.rmtree(tmp, ignore_errors=True)
    else:
        result["workdir"] = tmp
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xyz")
    ap.add_argument("--method", default="MNDO")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    ap.add_argument("--mode", default="1scf", choices=["1scf", "gradient", "force", "optimize"])
    ap.add_argument("--keyword", action="append", default=[], help="additional MOPAC keyword (repeatable)")
    ap.add_argument(
        "--displace",
        action="append",
        nargs=3,
        metavar=("ATOM", "AXIS", "DELTA_ANG"),
        default=[],
        help="displace one XYZ atom before the run (atom is 1-based; axis is x/y/z)",
    )
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()
    axis = {"x": 0, "y": 1, "z": 2}
    displacements = []
    for atom, direction, delta in args.displace:
        direction = direction.lower()
        if direction not in axis:
            ap.error(f"invalid displacement axis {direction!r}; use x, y, or z")
        if int(atom) < 1:
            ap.error("displacement atom indices are 1-based")
        displacements.append((int(atom) - 1, axis[direction], float(delta)))
    result = run(
        args.xyz,
        args.method,
        args.charge,
        args.mult,
        args.mode,
        args.keep,
        extra_keywords=args.keyword,
        displacements=displacements,
    )
    json.dump(result, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
