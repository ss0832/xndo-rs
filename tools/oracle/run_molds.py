#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Read MolDS CNDO/2 and INDO output, either bundled or freshly run.

xndo-rs's CNDO/2 and INDO engines are ports of MolDS 0.3.1, so MolDS is the
program those parameter tables came from and is the right oracle for them.

There are two ways to get its numbers, and this module handles both.

**Tier 1 -- the regression outputs MolDS ships.** `test/` in the 0.3.1 source
carries matched `.in`/`.dat` pairs written by the upstream author. They are the
implementer's own values, they carry no build risk at all, and they are licensed
GPL-3.0-or-later, the same as this crate. Five of them use CNDO/2 or INDO and
are vendored under `third_party/molds/tests/`.

**Tier 2 -- running a locally built MolDS.** Needed to choose the geometries and
so guarantee element coverage; the bundled set covers only what the author
happened to test. Set `MOLDS_EXE` and pass `run=True`.

Tier 1 comes first deliberately: MolDS is 2012 C++ against Boost and MPI, and
building it is the riskiest step in the whole oracle plan. Mining the bundled
outputs means a failed build costs coverage, not the whole method.

What a MolDS output carries, per calculation:

* the **whole** MO spectrum, occupied and virtual, in a.u. and eV;
* the electronic energy (which *includes* core repulsions -- MolDS says so on
  the line below it) and the core repulsion separately;
* the total, electronic and core dipole moments, in a.u. and Debye, as x/y/z
  and magnitude;
* Mulliken charges per atom.

That is far more per molecule than the MOPAC sets give, because MOPAC prints
only the frontier region by default.
"""

import os
import re
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
BUNDLED = os.path.join(ROOT, "third_party", "molds", "tests")
# Path to a MolDS binary. On Windows this is a WSL path and the run goes through
# `wsl`, because MolDS is a GNU/Linux build: Debian's packaging is the only thing
# that demonstrates a working toolchain for it, and trying to reproduce that on
# Windows buys nothing.
MOLDS_EXE = os.environ.get("MOLDS_EXE", "")
MOLDS_WSL = os.environ.get("MOLDS_WSL", "1" if os.name == "nt" else "0") != "0"

# MolDS writes atomic units and eV/Debye side by side. The eV and Debye columns
# are what the suite compares, because that is what xndo-rs reports; the a.u.
# columns are read as well so a unit mix-up shows up as a factor, not a fudge.
_MO = re.compile(
    r"Energy of MO:\s+(\d+)\s+(occ|unocc)\s+([-\d.e+]+)\s+([-\d.e+]+)"
)
_CHARGE = re.compile(
    r"Mulliken charge\(SCF\):\s+\d+\s+(\d+)\s+(\w+)\s+([-\d.e+]+)\s+([-\d.e+]+)"
)
_ATOM = re.compile(
    r"Atom coordinates:\s+(\d+)\s+(\w+)\s+"
    r"[-\d.e+]+\s+[-\d.e+]+\s+[-\d.e+]+\s+"
    r"([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)"
)


def _scalar(text, label):
    """Read one labelled `a.u. eV` pair, returning the eV value."""
    m = re.search(re.escape(label) + r":\s*([-\d.e+]+)\s+([-\d.e+]+)", text)
    return float(m.group(2)) if m else None


def parse_output(path):
    """Parse one MolDS output file into the quantities the suite compares."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    if "did not met convergence" in text or "Error" in text:
        raise RuntimeError("{}: MolDS reported an error".format(os.path.basename(path)))

    theory = None
    m = re.search(r"\*{10}\s+START:\s+([A-Za-z0-9/]+)-SCF", text)
    if m:
        theory = m.group(1)

    atoms = [
        (mm.group(2).capitalize(), float(mm.group(3)), float(mm.group(4)), float(mm.group(5)))
        for mm in _ATOM.finditer(text)
    ]
    orbitals = [
        (int(mm.group(1)), mm.group(2), float(mm.group(3)), float(mm.group(4)))
        for mm in _MO.finditer(text)
    ]
    charges = [
        (int(mm.group(1)), mm.group(2).capitalize(), float(mm.group(4)))
        for mm in _CHARGE.finditer(text)
    ]

    # The dipole block prints x, y, z, magnitude in a.u. then the same in Debye.
    dipole = None
    m = re.search(
        r"Total Dipole moment\(SCF\):\s+"
        r"([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+"
        r"([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)",
        text,
    )
    if m:
        dipole = [float(m.group(i)) for i in (5, 6, 7)]

    return {
        "theory": theory,
        "atoms": atoms,
        "n_valence_ao": _int_after(text, "Total number of valence AOs"),
        "n_valence_electrons": _int_after(text, "Total number of valence electrons"),
        # MolDS's "Electronic energy(SCF)" already contains the core repulsions;
        # the file says so on the next line. Subtracting the separately printed
        # core repulsion gives the purely electronic part, which is what
        # xndo-rs calls `electronic_ev`.
        "total_energy_ev": _scalar(text, "Electronic energy(SCF)"),
        "core_repulsion_ev": _scalar(text, "Core repulsion energy"),
        "eigenvalues_ev": [ev for _, _, _, ev in orbitals],
        "occupations": [occ for _, occ, _, _ in orbitals],
        "n_occupied": sum(1 for _, occ, _, _ in orbitals if occ == "occ"),
        "charges": [q for _, _, q in charges],
        "dipole_debye": dipole,
    }


def _int_after(text, label):
    m = re.search(re.escape(label) + r":\s*(\d+)", text)
    return int(m.group(1)) if m else None


def bundled_cases():
    """The vendored upstream regression cases, as {name: (theory, path)}."""
    out = {}
    if not os.path.isdir(BUNDLED):
        return out
    for entry in sorted(os.listdir(BUNDLED)):
        if not entry.endswith(".dat"):
            continue
        stem = entry[:-4]
        theory = "INDO" if stem.endswith("_indo") else "CNDO/2"
        out[stem] = (theory, os.path.join(BUNDLED, entry))
    return out


def _wsl_path(path):
    return subprocess.run(
        ["wsl", "-e", "wslpath", "-a", path.replace("\\", "/")],
        capture_output=True, text=True, check=True,
    ).stdout.strip()


def write_input(path, atoms, theory="cndo/2", rms_density="0.000000001"):
    """Write a MolDS input deck.

    `rms_density` is tightened well past the 1e-6 the bundled regression cases
    use. The same lesson as ORACLE_NOTES item 22: an oracle converged only as far
    as its default leaves a residual that looks like a model difference, and it
    costs nothing to remove.
    """
    with open(path, "w", newline="\n") as fh:
        fh.write("THEORY\n   {}\nTHEORY_END\n\n".format(theory))
        fh.write(
            "SCF\n   max_iter 500\n   rms_density {}\n"
            "   damping_thresh 1.0\n   damping_weight 0.0\n"
            "   diis_num_error_vect 5\n   diis_start_error 0.1\n"
            "   diis_end_error 0.00000002\nSCF_END\n\n".format(rms_density)
        )
        fh.write("GEOMETRY\n")
        for sym, x, y, z in atoms:
            fh.write(" {:<3s} {:>16.8f} {:>16.8f} {:>16.8f}\n".format(sym, x, y, z))
        fh.write("GEOMETRY_END\n")


def run(atoms, theory="cndo/2", workdir=None, rms_density="0.000000001"):
    """Run a locally built MolDS on a geometry (Tier 2)."""
    if not MOLDS_EXE:
        raise RuntimeError(
            "MOLDS_EXE is not set to a MolDS binary; Tier 2 needs a local build "
            "(see tests/data/ORACLE_NOTES.md item 24)"
        )
    tmp = workdir or tempfile.mkdtemp(prefix="xndors_molds_")
    os.makedirs(tmp, exist_ok=True)
    inp = os.path.join(tmp, "job.in")
    out = os.path.join(tmp, "job.dat")
    write_input(inp, atoms, theory, rms_density)
    if MOLDS_WSL:
        subprocess.run(
            ["wsl", "-e", "bash", "-c",
             "'{}' < '{}' > '{}'".format(MOLDS_EXE, _wsl_path(inp), _wsl_path(out))],
            check=True, capture_output=True,
        )
    else:
        with open(inp) as fin, open(out, "w") as fout:
            subprocess.run([MOLDS_EXE], stdin=fin, stdout=fout, check=True, cwd=tmp)
    return parse_output(out)


def main():
    import argparse
    import json

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("path", nargs="?", help="a MolDS .dat output to parse")
    ap.add_argument("--list", action="store_true", help="list the bundled cases")
    args = ap.parse_args()
    if args.list or not args.path:
        for name, (theory, path) in bundled_cases().items():
            parsed = parse_output(path)
            print(
                "{:<16s} {:<7s} {:>2d} atoms  {:>2d} AOs  {:>2d} occ".format(
                    name,
                    theory,
                    len(parsed["atoms"]),
                    parsed["n_valence_ao"] or 0,
                    parsed["n_occupied"],
                )
            )
        return 0
    print(json.dumps(parse_output(args.path), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
