#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Drive MOPAC7 1.15 for MINDO/3 reference values.

MOPAC 23.2.5 does not implement MINDO/3 -- `molkst_C.F90:331-367` dispatches
twenty methods and MINDO/3 is not among them -- so the modern binary cannot be
the oracle here. MOPAC7 1.15 can: it is public domain, it is the program
`src/data/mindo3_parameters.csv` names as its upstream, and its `fortran/block.f`
is where those numbers come from. It is the source rather than a proxy.

The alternative was transcribing published MINDO/3 tables, which is strictly
weaker: the papers give values at MINDO/3-optimised geometries in internal
coordinates, so SCF and optimiser error are mixed together and nothing better
than about 0.1 kcal/mol is recoverable.

Two things about driving it:

* it takes no command-line arguments. File assignment is through the `FOR005`
  ... `FOR012` environment variables, the Fortran unit convention of its era.
* this build prints `== symtrz.f ... ==` debug lines into the output, left in by
  whoever last touched the symmetry code. They are filtered before parsing.

MOPAC7 prints more than MOPAC 23 does by default: the total energy, the
electronic energy and the core-core repulsion separately, which is what makes it
possible to say whether a disagreement is in the density or in the geometry
term.
"""

import math
import os
import re
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
# A MOPAC7 binary. On Windows this is a WSL path: MOPAC7 is Fortran 77 and the
# GNU toolchain is what Debian's packaging demonstrates works for it.
MOPAC7_EXE = os.environ.get("MOPAC7_EXE", "")
MOPAC7_WSL = os.environ.get("MOPAC7_WSL", "1" if os.name == "nt" else "0") != "0"

_DEBUG_LINE = re.compile(r"^\s*==\s*\w+\.f\s")


def _distance(a, b):
    return math.sqrt(sum((a[i + 1] - b[i + 1]) ** 2 for i in range(3)))


def _angle_degrees(a, b, c):
    """The a-b-c angle in degrees, b at the vertex."""
    u = [a[i + 1] - b[i + 1] for i in range(3)]
    v = [c[i + 1] - b[i + 1] for i in range(3)]
    nu = math.sqrt(sum(t * t for t in u))
    nv = math.sqrt(sum(t * t for t in v))
    dot = sum(u[i] * v[i] for i in range(3)) / (nu * nv)
    return math.degrees(math.acos(max(-1.0, min(1.0, dot))))


def first_three_are_collinear(atoms, tolerance_degrees=1.0e-4):
    """Whether MOPAC7 will reject these atoms as Cartesian input.

    `getgeo.f:373` converts Cartesian input to internal coordinates and then
    stops outright -- "DUE TO PROGRAM BUG, THE FIRST THREE ATOMS MUST NOT LIE IN
    A STRAIGHT LINE" -- if the 1-2-3 angle comes back as 0 or 180 degrees, to
    within 1e-4. The same tolerance is used here so this predicts the rejection
    rather than approximating it.
    """
    if len(atoms) < 3:
        return True
    theta = _angle_degrees(atoms[0], atoms[1], atoms[2])
    return abs(theta) < tolerance_degrees or abs(theta - 180.0) < tolerance_degrees


def write_input(path, atoms, keywords="MINDO3 1SCF PRECISE", comment="xndo-rs oracle"):
    """Write a MOPAC7 deck for these atoms.

    Which of MOPAC7's two input forms is used is not a preference; each is the
    only one that works in its own case, and `getgeo.f` says where the line is:

    * **Three atoms or fewer: internal coordinates, always.** `getgeo.f:256`
      reads `IF(NATOMS.GT.3) THEN INT=(NA(4).NE.0) ELSE INT=.TRUE.` -- with
      three atoms or fewer MOPAC7 assumes a Z-matrix unconditionally, and no
      keyword overrides it (`XYZ` does not; it is consulted later, at `:345`,
      and only affects which flags get cleared). Handing it Cartesians there
      does not fail, which is the dangerous part: it reads the three numbers as
      a bond length, an angle and a dihedral and computes a *different
      molecule*. Water came back as a straight O-H-H chain with a 0.7575 A bond
      that way, converged, and reported an energy.
    * **Four atoms or more: Cartesian**, which `NA(4) = 0` selects. A Z-matrix
      would do here too, but only if every atom has three non-collinear
      predecessors, and building one for an arbitrary molecule means inserting
      dummy atoms. Cartesian input needs none of that.

    A Z-matrix for three atoms or fewer never needs a dihedral and never has an
    undefined angle, so the linear-molecule problem that makes Z-matrices awkward
    in general does not arise at this size. `getgeo.f:373`'s collinearity check
    lives inside the Cartesian branch, so a linear triatomic goes through the
    internal-coordinate path untouched.

    MOPAC7 rebuilds its own Cartesian frame from a Z-matrix, so for these
    molecules the orientation is MOPAC7's and not the caller's. Every quantity
    the MINDO/3 suite compares is invariant under rigid motion -- energies,
    charges, orbital energies, and the dipole *magnitude* -- which is why that
    suite compares the magnitude rather than the components.
    """
    with open(path, "w", newline="\n") as fh:
        fh.write(keywords + "\n" + comment + "\n\n")
        if len(atoms) <= 3:
            for i, atom in enumerate(atoms):
                fields = [0.0, 0.0, 0.0]
                na = nb = nc = 0
                if i >= 1:
                    fields[0] = _distance(atom, atoms[0] if i == 1 else atoms[i - 1])
                    na = 1 if i == 1 else i
                if i == 2:
                    fields[1] = _angle_degrees(atom, atoms[1], atoms[0])
                    na, nb = 2, 1
                fh.write(
                    "  {:<3s} {:>18.12f} 0 {:>18.12f} 0 {:>18.12f} 0 {:>4d}{:>4d}{:>4d}\n".format(
                        atom[0], fields[0], fields[1], fields[2], na, nb, nc
                    )
                )
        else:
            for sym, x, y, z in atoms:
                fh.write(
                    "  {:<3s} {:>18.12f} 1 {:>18.12f} 1 {:>18.12f} 1    0    0    0\n".format(
                        sym, x, y, z
                    )
                )
        fh.write("\n")


def parse_output(path):
    """Parse a MOPAC7 `.OUT` into the quantities the suite compares."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [l for l in fh if not _DEBUG_LINE.match(l)]
    text = "".join(lines)
    if "JOB ENDED NORMALLY" not in text and "DONE" not in text:
        raise RuntimeError("MOPAC7 did not finish: " + text.strip()[-200:])
    for bad in ("UNABLE TO ACHIEVE SELF-CONSISTENCE", "NOT ALLOWED", "IS NOT AVAILABLE"):
        if bad in text:
            raise RuntimeError("MOPAC7 reported: " + bad)

    def scalar(label):
        m = re.search(re.escape(label) + r"\s*=\s*(-?[\d.]+)", text)
        return float(m.group(1)) if m else None

    eigenvalues = []
    for i, line in enumerate(lines):
        if "EIGENVALUES" in line:
            for follow in lines[i + 1:]:
                toks = follow.split()
                if not toks:
                    if eigenvalues:
                        break
                    continue
                try:
                    eigenvalues.extend(float(t) for t in toks)
                except ValueError:
                    break
            break

    charges = []
    m = re.search(r"NET ATOMIC CHARGES.*?\n(.*?)(?:\n\s*\n)", text, re.S)
    if m:
        for line in m.group(1).splitlines():
            fields = line.split()
            if len(fields) == 4 and fields[0].isdigit():
                charges.append(float(fields[2]))

    dipole = None
    dipole_total = None
    m = re.search(r"^ SUM\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", text, re.M)
    if m:
        dipole = [float(m.group(i)) for i in (1, 2, 3)]
        # The fourth column. Frame-independent, unlike the three before it, and
        # so the one a Z-matrix deck can still be compared on -- see
        # `write_input` on whose frame those components are in.
        dipole_total = float(m.group(4))

    return {
        "heat_of_formation_kcal": scalar("FINAL HEAT OF FORMATION"),
        "total_energy_ev": scalar("TOTAL ENERGY"),
        "electronic_energy_ev": scalar("ELECTRONIC ENERGY"),
        "core_repulsion_ev": scalar("CORE-CORE REPULSION"),
        "ionization_potential_ev": scalar("IONIZATION POTENTIAL"),
        "n_occupied": int(scalar("NO. OF FILLED LEVELS") or 0) or None,
        "eigenvalues_ev": eigenvalues or None,
        "charges": charges or None,
        "dipole_debye": dipole,
        "dipole_magnitude_debye": dipole_total,
    }


def _wsl_path(path):
    return subprocess.run(
        ["wsl", "-e", "wslpath", "-a", path.replace("\\", "/")],
        capture_output=True, text=True, check=True,
    ).stdout.strip()


#: Multiplicity keyword for MOPAC7's `UHF` open-shell path. MOPAC7 names the
#: multiplicity in words rather than as a number, and gives no keyword for a
#: closed-shell singlet -- that is the default.
MULTIPLICITY_KEYWORD = {
    1: "",
    2: "UHF DOUBLET",
    3: "UHF TRIPLET",
    4: "UHF QUARTET",
    5: "UHF QUINTET",
}


#: MOPAC7 stops on any interatomic distance below this, unless `GEO-OK` is
#: given (`moldat.f:638`). It is a guard against a mistyped geometry, not a
#: model limit -- H2's bond is 0.741 A and trips it.
MIN_DISTANCE_ANGSTROM = 0.8


def keywords_for(charge=0, mult=1, base="MINDO3 1SCF PRECISE", atoms=None):
    """The MOPAC7 keyword line for one calculation.

    Deliberately **without** `VECTORS`. It reads like the way to get more, but it
    replaces the plain `EIGENVALUES` block with `ALPHA`/`BETA EIGENVECTORS`
    tables -- so asking for the eigenvectors is how you stop being given the
    eigenvalues in the form anything parses. Without it MOPAC7 prints
    `EIGENVALUES` for closed and open shells alike, the open-shell block being
    the alpha channel, which is the one `mo_energies_ev` holds on this side too.
    """
    parts = [base]
    if mult not in MULTIPLICITY_KEYWORD:
        raise ValueError("MOPAC7 has no keyword for multiplicity {}".format(mult))
    word = MULTIPLICITY_KEYWORD[mult]
    if word:
        parts.append(word)
    if charge:
        parts.append("CHARGE={}".format(int(charge)))
    if atoms and _closest_pair(atoms) < MIN_DISTANCE_ANGSTROM:
        # Added only where MOPAC7's proximity guard would actually fire, so it
        # stays armed for every other molecule. Switching it on unconditionally
        # would be the easy thing and would also switch off the one check that
        # catches a mistyped coordinate.
        parts.append("GEO-OK")
    return " ".join(parts)


def _closest_pair(atoms):
    """The shortest interatomic distance, or infinity for a single atom."""
    return min(
        (_distance(atoms[i], atoms[j])
         for i in range(len(atoms)) for j in range(i + 1, len(atoms))),
        default=float("inf"),
    )


def run(atoms, keywords="MINDO3 1SCF PRECISE", workdir=None):
    if not MOPAC7_EXE:
        raise RuntimeError(
            "MOPAC7_EXE is not set to a MOPAC7 binary; see tests/data/ORACLE_NOTES.md item 25"
        )
    tmp = workdir or tempfile.mkdtemp(prefix="xndors_mopac7_")
    os.makedirs(tmp, exist_ok=True)
    dat = os.path.join(tmp, "job.dat")
    out = os.path.join(tmp, "job.OUT")
    write_input(dat, atoms, keywords)
    if MOPAC7_WSL:
        wsl_dir = _wsl_path(tmp)
        script = (
            "cd '{d}' && FOR005=job.dat FOR006=job.out FOR009=job.res "
            "FOR010=job.den FOR011=job.log FOR012=job.arc '{exe}' > job.OUT 2>&1"
        ).format(d=wsl_dir, exe=MOPAC7_EXE)
        subprocess.run(["wsl", "-e", "bash", "-c", script], check=True, capture_output=True)
    else:
        env = dict(os.environ)
        env.update({
            "FOR005": "job.dat", "FOR006": "job.out", "FOR009": "job.res",
            "FOR010": "job.den", "FOR011": "job.log", "FOR012": "job.arc",
        })
        with open(out, "w") as fh:
            subprocess.run([MOPAC7_EXE], stdout=fh, stderr=fh, cwd=tmp, env=env, check=True)
    return parse_output(out)


def main():
    import argparse
    import json

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("xyz")
    ap.add_argument("--keywords", default="MINDO3 1SCF PRECISE")
    args = ap.parse_args()
    with open(args.xyz) as fh:
        lines = fh.read().splitlines()
    atoms = []
    for line in lines[2:]:
        f = line.split()
        if len(f) >= 4:
            atoms.append((f[0], float(f[1]), float(f[2]), float(f[3])))
    print(json.dumps(run(atoms, args.keywords), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
