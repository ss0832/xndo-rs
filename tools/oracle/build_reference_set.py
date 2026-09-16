#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Generate the MOPAC-backed oracle reference sets consumed by tests/oracle_matrix.rs.

Every value column is what MOPAC printed. Every geometry column is what was sent
to it. Nothing here is typed by a human except the geometry *parameters* and the
``why`` note on each row, both of which a reviewer can evaluate directly.

Usage::

    set MOPAC_EXE=...\\mopac.exe
    python tools/oracle/build_reference_set.py --method mndo
    python tools/oracle/build_reference_set.py --method all
    python tools/oracle/build_reference_set.py --method all --check   # fail if stale
    python tools/oracle/build_reference_set.py --method mndo --report # coverage only

``--check`` re-runs the oracle and byte-compares the rendered file, so it needs
the binary. Without it every molecule is skipped, and a file rendered from no
rows is reported as *unchecked* rather than stale -- the committed data is not
known to be wrong, it simply was not compared against anything. The hash
manifest in tests/data/MANIFEST.sha256 is the layer that detects a hand-edited
file without needing MOPAC; see tests/data/ORACLE_NOTES.md section (d).
"""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

import geometry_builders as gb  # noqa: E402
import run_mopac  # noqa: E402
import run_molds  # noqa: E402
import run_mopac7  # noqa: E402

OUT_DIR = os.path.join(ROOT, "tests", "data")


# --- the molecule sets -------------------------------------------------------
#
# Each entry: (name, charge, multiplicity, atoms, why)
#
# The `why` note is not decoration. It is the record of what the row is for, and
# the thing a reviewer checks when deciding whether the set still covers what it
# claims after an edit.


def mndo_molecules():
    m = []
    a = m.append

    # Diatomics: the smallest systems that isolate one bond each.
    a(("h2", 0, 1, gb.linear("H", ("H", 0.741)), "the smallest closed shell"))
    a(("n2", 0, 1, gb.linear("N", ("N", 1.098)), "triple bond, second row"))
    a(("o2", 0, 3, gb.linear("O", ("O", 1.208)), "ground-state triplet that is not a radical"))
    a(("f2", 0, 1, gb.linear("F", ("F", 1.412)), "the most electronegative pair"))
    a(("cl2", 0, 1, gb.linear("Cl", ("Cl", 1.988)), "third row"))
    a(("br2", 0, 1, gb.linear("Br", ("Br", 2.281)), "fourth row halogen"))
    a(("i2", 0, 1, gb.linear("I", ("I", 2.666)), "fifth row halogen"))
    a(("co", 0, 1, gb.linear("C", ("O", 1.128)), "polar triple bond"))
    a(("no", 0, 2, gb.linear("N", ("O", 1.151)), "doublet radical"))
    a(("hf", 0, 1, gb.linear("H", ("F", 0.917)), "largest dipole per bond in the second row"))
    a(("hcl", 0, 1, gb.linear("H", ("Cl", 1.275)), "H on a d-bearing partner in MNDO/d"))
    a(("hbr", 0, 1, gb.linear("H", ("Br", 1.414)), "fourth-row hydride"))
    a(("hi", 0, 1, gb.linear("H", ("I", 1.609)), "fifth-row hydride"))
    a(("lih", 0, 1, gb.linear("Li", ("H", 1.595)), "electropositive metal hydride"))
    a(("nah", 0, 1, gb.linear("Na", ("H", 1.887)), "third-row alkali hydride"))
    a(("bh", 0, 1, gb.linear("B", ("H", 1.232)), "boron, an electron-poor centre"))

    # Hydrides: one shape per row of the p block.
    a(("water", 0, 1, gb.bent("O", "H", "H", 0.9584, 0.9584, 104.45), "the reference small polar molecule"))
    a(("h2s", 0, 1, gb.bent("S", "H", "H", 1.336, 1.336, 92.1), "sulfur analogue, much smaller angle"))
    a(("ammonia", 0, 1, gb.pyramidal("N", "H", 1.012, 106.7), "pyramidal, lone pair"))
    a(("ph3", 0, 1, gb.pyramidal("P", "H", 1.420, 93.3), "phosphorus, nearly orthogonal bonds"))
    a(("methane", 0, 1, gb.tetrahedral("C", "H", 1.087), "the reference saturated carbon"))
    a(("silane", 0, 1, gb.tetrahedral("Si", "H", 1.480), "silicon, second-row d-bearing in MNDO/d"))
    a(("germane", 0, 1, gb.tetrahedral("Ge", "H", 1.525), "germanium"))
    a(("beh2", 0, 1, gb.linear("H", ("Be", 1.334), ("H", 1.334)), "linear, electron-deficient"))

    # Polyatomics: shapes and bonding motifs.
    a(("co2", 0, 1, gb.linear("O", ("C", 1.160), ("O", 1.160)), "linear, two polar double bonds"))
    a(("cs2", 0, 1, gb.linear("S", ("C", 1.553), ("S", 1.553)), "sulfur analogue of CO2"))
    a(("so2", 0, 1, gb.bent("S", "O", "O", 1.431, 1.431, 119.5), "hypervalent-looking bent sulfur"))
    a(("o3", 0, 1, gb.bent("O", "O", "O", 1.278, 1.278, 116.8), "strained homonuclear bent"))
    a(("hcn", 0, 1, gb.linear("H", ("C", 1.065), ("N", 1.153)), "linear with a triple bond"))
    a(("formaldehyde", 0, 1, [gb.atom("C"), gb.atom("O", 0, 0, 1.203),
                              gb.atom("H", 0.943, 0, -0.545), gb.atom("H", -0.943, 0, -0.545)],
       "carbonyl; also the ZINDO/S n->pi* reference"))
    a(("ethene", 0, 1, [gb.atom("C", 0, 0, 0.667), gb.atom("C", 0, 0, -0.667),
                        gb.atom("H", 0.929, 0, 1.237), gb.atom("H", -0.929, 0, 1.237),
                        gb.atom("H", 0.929, 0, -1.237), gb.atom("H", -0.929, 0, -1.237)],
       "C=C double bond, planar"))
    a(("ethyne", 0, 1, gb.linear("H", ("C", 1.061), ("C", 1.203), ("H", 1.061)),
       "C#C triple bond, linear"))
    a(("benzene", 0, 1, gb.ring(["C"] * 6, 1.397, "H", 1.084), "aromatic ring, delocalised"))
    a(("bf3", 0, 1, gb.trigonal_planar("B", ["F"] * 3, [1.313] * 3), "trigonal planar, no lone pair"))
    a(("nh4+", 1, 1, gb.tetrahedral("N", "H", 1.021), "closed-shell cation"))
    a(("bh4-", -1, 1, gb.tetrahedral("B", "H", 1.236), "closed-shell anion"))
    a(("ch3f", 0, 1, gb.tetrahedral("C", "H", 1.087, [(0, "F", 1.383)]),
       "substituted carbon; the substituent carries its own bond length"))
    a(("ch3cl", 0, 1, gb.tetrahedral("C", "H", 1.087, [(0, "Cl", 1.785)]), "heavier halide"))
    a(("ccl4", 0, 1, gb.tetrahedral("C", "Cl", 1.767), "four heavy substituents"))
    a(("pf5", 0, 1, gb.bipyramidal("P", "F", 1.577, "F", 1.522), "trigonal bipyramidal, hypervalent"))
    a(("sf6", 0, 1, gb.octahedral("S", "F", 1.564), "octahedral, the hypervalent stress case"))

    # Open shells: the UHF path.
    a(("methyl_radical", 0, 2, [gb.atom("C")] + gb.methyl_group((0.0, 0.0, 0.0), (0, 0, 1), 1.079, 120.0),
       "planar doublet radical"))
    a(("oh_radical", 0, 2, gb.linear("O", ("H", 0.970)), "diatomic doublet, degenerate surface"))
    a(("nh2_radical", 0, 2, gb.bent("N", "H", "H", 1.024, 1.024, 103.4), "bent doublet"))
    a(("ch2_triplet", 0, 3, gb.bent("C", "H", "H", 1.078, 1.078, 133.9), "triplet carbene"))
    a(("no2", 0, 2, gb.bent("N", "O", "O", 1.193, 1.193, 134.1), "doublet with two equivalent bonds"))

    # d elements. The eight MNDO elements that carry an independent `poc` --
    # Sc, V, Cr, Fe, Mo, Pd, Ag, Pt -- deliberately have no row here: MOPAC
    # 23.2.5 refuses to run them under MNDO ("Parameters for some elements are
    # missing"), even though its own parameters_for_mndo_C.F90 carries values
    # that xndo-rs's table reproduces. See ORACLE_NOTES item 16. Zn does run,
    # and is a d element whose poc equals the derived po(1), so it checks the
    # d path without depending on the poc override.
    a(("zncl2", 0, 1, gb.linear("Cl", ("Zn", 2.072), ("Cl", 2.072)),
       "zinc, a d element whose poc equals the derived po(1)"))

    # Alkali and alkaline earth, where MNDO is often weakest.
    a(("licl", 0, 1, gb.linear("Li", ("Cl", 2.021)), "ionic bond"))
    a(("nacl", 0, 1, gb.linear("Na", ("Cl", 2.361)), "the classic ionic diatomic"))
    a(("mgh2", 0, 1, gb.linear("H", ("Mg", 1.703), ("H", 1.703)), "alkaline earth hydride"))
    a(("alh3", 0, 1, gb.trigonal_planar("Al", ["H"] * 3, [1.584] * 3), "aluminium, trigonal planar"))
    return m


def mndod_molecules():
    """MNDO/d covers 22 elements; the set is built around the ones with d shells."""
    keep = {
        "h2", "water", "h2s", "ammonia", "ph3", "methane", "silane", "hcl", "hbr", "hi",
        "cl2", "br2", "i2", "co2", "cs2", "so2", "formaldehyde", "ethene", "benzene",
        "bf3", "ch3cl", "ccl4", "pf5", "sf6", "nh4+", "bh4-", "lih", "nah", "beh2",
        "mgh2", "alh3", "methyl_radical", "oh_radical", "nh2_radical", "no2", "f2", "hf",
    }
    return [row for row in mndo_molecules() if row[0] in keep]


def zindo_molecules():
    """ZINDO/S ground state.

    The s/p branch covers H, Li, Be, B, C, N, O, F, Na, Mg, Al, Si, P, S, Cl, K,
    Zn, Br and I (`zindo_s_parameters.csv` rows with `norb` in {1, 4}); the 9-AO
    transition-metal branch is rejected by the engine, so nothing here needs one.

    ZINDO/S is a spectroscopic parameterisation with no atomic heat terms, so
    what MOPAC prints as a heat of formation is a transformed total energy (see
    ORACLE_NOTES item 10). This set therefore compares the **orbital spectrum**,
    the charges and the dipole, which are convention-free, and leaves energies to
    the excited-state suite.
    """
    keep = {
        "h2", "n2", "f2", "cl2", "br2", "i2", "co", "hf", "hcl", "hbr", "hi", "lih",
        "nah", "bh", "water", "h2s", "ammonia", "ph3", "methane", "silane", "beh2",
        "co2", "cs2", "so2", "o3", "hcn", "formaldehyde", "ethene", "ethyne",
        "benzene", "bf3", "nh4+", "bh4-", "ch3f", "ch3cl", "ccl4", "licl", "nacl",
        "mgh2", "alh3", "zncl2",
    }
    rows = [row for row in mndo_molecules() if row[0] in keep]
    # ZINDO/S is parameterised for spectroscopy of closed-shell chromophores;
    # the open-shell rows live in the UCIS suite, not here.
    return [row for row in rows if row[2] == 1]


def zindo_cis_molecules():
    """ZINDO/S singlet CIS.

    Closed-shell chromophores with enough occupied and virtual orbitals for a
    window worth comparing. The elements are the same s/p set the ground-state
    suite covers; what is new here is the excited-state machinery -- the CI
    matrix, its eigenvectors, the transition dipoles and the relaxed state
    densities -- none of which the ground-state suite touches at all.
    """
    keep = {
        # carbonyl and azine n->pi*, the classic ZINDO/S test cases
        "formaldehyde", "co2", "cs2", "so2", "o3", "hcn",
        # unsaturated carbon: pi->pi*
        "ethene", "ethyne", "benzene",
        # diatomics and small hydrides with a real virtual space
        "n2", "co", "f2", "cl2", "br2", "i2", "hcl", "hbr", "hi",
        # heteroatom lone pairs
        "h2s", "ph3", "ammonia", "water",
        # halomethanes: sigma* virtuals
        "ch3f", "ch3cl", "ccl4",
        # ionic, to exercise a charged reference
        "licl", "nacl", "zncl2", "bf3",
        # closed-shell ions: the CI is built on a charged SCF reference, which
        # nothing else in this set does
        "nh4+", "bh4-",
    }
    rows = [row for row in mndo_molecules() if row[0] in keep]
    return [row for row in rows if row[2] == 1]


# MolDS ships parameters for a small element set, and the two theories do not
# cover the same one. `third_party/molds/NOTICE` records that the sulfur row is
# not a complete INDO one-centre set, so INDO stops at oxygen.
CNDO2_ELEMENTS = {"H", "Li", "C", "N", "O", "S"}
INDO_ELEMENTS = {"H", "Li", "C", "N", "O"}


def _molds_molecules(allowed):
    """Closed-shell neutral and ionic molecules inside an element set.

    Drawn from the same table the MOPAC sets use, so a molecule that appears in
    both suites is the same geometry in both, and a disagreement between the two
    oracles is about the model rather than about the coordinates.
    """
    out = []
    for name, charge, mult, atoms, why in mndo_molecules():
        if mult != 1:
            continue
        if not {sym for sym, _, _, _ in atoms} <= allowed:
            continue
        out.append((name, charge, mult, atoms, why))
    return out


def cndo2_molecules():
    """CNDO/2 against a locally built MolDS (Tier 2)."""
    return _molds_molecules(CNDO2_ELEMENTS)


def indo_molecules():
    """INDO against a locally built MolDS (Tier 2)."""
    return _molds_molecules(INDO_ELEMENTS)


# MINDO/3's published element set, and what `src/data/mindo3_parameters.csv`
# carries: H, B, C, N, O, F, Si, P, S, Cl. There are no d functions and no
# fourth row, so the MNDO set's Br/I/Ge/Zn/alkali rows drop out.
MINDO3_ELEMENTS = {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl"}
SYMBOL_TO_Z = {"H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
               "Si": 14, "P": 15, "S": 16, "Cl": 17}


def _mindo3_pairs():
    """The atom pairs `mindo3_pair_parameters.csv` actually carries.

    Read from the shipped table rather than listed here, so the set cannot claim
    a molecule the parameters do not support and cannot silently stop covering
    one that they do.
    """
    path = os.path.join(ROOT, "src", "data", "mindo3_pair_parameters.csv")
    pairs = set()
    with open(path, encoding="utf-8") as fh:
        rows = [line for line in fh if not line.startswith("#")]
    reader = csv.DictReader(rows)
    for row in reader:
        z1, z2 = int(row["z1"]), int(row["z2"])
        pairs.add((min(z1, z2), max(z1, z2)))
    return pairs


def mindo3_molecules():
    """MINDO/3 against MOPAC7, the program its parameter tables come from.

    Same table as the MOPAC sets, filtered twice: to the elements MINDO/3 has,
    and then to the molecules whose every atom pair has resonance and core
    parameters. MINDO/3's pair table covers 40 of the 55 pairs its ten elements
    could form -- F-S is one of the fifteen it does not -- so SF6 is outside the
    published model rather than a case this engine fails. Leaving it in the set
    to be skipped at run time would report it as a failure, which is a different
    claim.

    A molecule appearing in two suites is the same geometry in both, so a
    disagreement between two oracles is about the model and not the coordinates.
    """
    pairs = _mindo3_pairs()
    out = []
    for row in mndo_molecules():
        symbols = {sym for sym, _, _, _ in row[3]}
        if not symbols <= MINDO3_ELEMENTS:
            continue
        zs = sorted({SYMBOL_TO_Z[s] for s in symbols})
        if all((a, b) in pairs for i, a in enumerate(zs) for b in zs[i:]):
            out.append(row)
    return out


SETS = {
    "mndo": ("MNDO", "mopac_mndo_reference.tsv", mndo_molecules),
    # MOPAC 23.2.5 does not implement MINDO/3 -- `molkst_C.F90:331-367`
    # dispatches twenty methods and MINDO/3 is not one of them -- so the modern
    # binary cannot be the oracle. MOPAC7 1.15 can, and is the declared upstream
    # of the MINDO/3 tables rather than a proxy for it.
    "mindo3": ("MINDO3", "mopac7_mindo3_reference.tsv", mindo3_molecules),
    "mndod": ("MNDOD", "mopac_mndod_reference.tsv", mndod_molecules),
    "zindo": ("INDO", "mopac_zindo_reference.tsv", zindo_molecules),
    "zindo_cis": ("INDO", "mopac_zindo_cis_reference.tsv", zindo_cis_molecules),
    # MolDS is the program the CNDO/2 and INDO engines were ported from, so it
    # is the right oracle for them. These two sets are mined from the regression
    # outputs MolDS ships rather than from a local build; see `run_molds`.
    "cndo2": ("cndo/2", "molds_cndo2_reference.tsv", cndo2_molecules),
    "indo": ("indo", "molds_indo_reference.tsv", indo_molecules),
}

# The CI window. MOPAC takes `C.I.=(N, M)`: N orbitals in the window, M of them
# occupied. Both sides are pinned to the same window per molecule and the window
# is written into the reference file, because comparing CIS energies across
# different windows compares nothing.
#
# Five and five, clipped to what the molecule actually has. Five is the largest
# window every molecule in the set can supply on the virtual side, and 25
# configurations is enough that the roots are genuinely mixed rather than being
# relabelled configurations.
CIS_MAX_OCCUPIED = 5
CIS_MAX_VIRTUAL = 5
# Roots compared per molecule. MOPAC prints every root in the window; the low
# ones are the ones with any spectroscopic meaning, and the high ones are
# dominated by the window edge.
CIS_ROOTS = 8

COLUMNS = [
    "name", "charge", "multiplicity", "n_atoms", "geometry",
    "hof_kcal", "ionization_potential_ev", "homo_ev", "lumo_ev",
    "dipole_x_debye", "dipole_y_debye", "dipole_z_debye",
    "charges", "eigenvalues_ev", "why",
]

MOLDS_COLUMNS = [
    "name", "charge", "multiplicity", "n_atoms", "geometry",
    "total_energy_ev", "core_repulsion_ev",
    "charges", "eigenvalues_ev",
    "dipole_x_debye", "dipole_y_debye", "dipole_z_debye",
    "why",
]

CIS_COLUMNS = [
    "name", "charge", "multiplicity", "n_atoms", "geometry",
    "active_occupied", "active_virtual",
    "excitation_energies_ev", "oscillator_strengths", "state_dipoles_debye",
    "why",
]


def compute(method_keyword, molecules):
    rows = []
    version = None
    for name, charge, mult, atoms, why in molecules:
        xyz = os.path.join(OUT_DIR, "_scratch.xyz")
        os.makedirs(OUT_DIR, exist_ok=True)
        with open(xyz, "w") as f:
            f.write(gb.to_xyz_text(atoms, name))
        try:
            result = run_mopac.run(xyz, method=method_keyword, charge=charge,
                                   mult=mult, mode="1scf")
        except Exception as exc:  # noqa: BLE001 - report and skip, do not fabricate
            print("  SKIP {}: {}".format(name, exc), file=sys.stderr)
            continue
        finally:
            if os.path.exists(xyz):
                os.remove(xyz)

        version = version or result.get("mopac_version")
        eig = result.get("eigenvalues_ev") or []
        homo = lumo = None
        if eig:
            # MOPAC orders eigenvalues ascending; the occupied count follows from
            # the electron count, which run_mopac reports via the AUX block.
            n_occ = result.get("n_occupied")
            if n_occ:
                homo = eig[n_occ - 1] if n_occ - 1 < len(eig) else None
                lumo = eig[n_occ] if n_occ < len(eig) else None
        dip = result.get("dipole_debye") or [None, None, None]
        rows.append({
            "name": name,
            "charge": charge,
            "multiplicity": mult,
            "n_atoms": len(atoms),
            "geometry": gb.to_xyz_field(atoms),
            "hof_kcal": result.get("heat_of_formation_kcal"),
            "ionization_potential_ev": result.get("ionization_potential_ev"),
            "homo_ev": homo,
            "lumo_ev": lumo,
            "dipole_x_debye": dip[0], "dipole_y_debye": dip[1], "dipole_z_debye": dip[2],
            "charges": result.get("charges"),
            "eigenvalues_ev": eig,
            "why": why,
        })
    return rows, version


def cis_window_candidates(n_occ, n_virtual):
    """Windows to try, largest first.

    Not a single fixed window, because MOPAC will not always diagonalise one.
    `CIS_MAX_*` bounds the search; everything below that is tried in descending
    order of configuration count and the first window MOPAC treats *fully* wins.
    """
    occ_max = min(n_occ, CIS_MAX_OCCUPIED)
    vir_max = min(n_virtual, CIS_MAX_VIRTUAL)
    windows = [
        (o, v)
        for o in range(2, occ_max + 1)
        for v in range(2, vir_max + 1)
    ]
    # Biggest first; break ties towards the squarer window, which mixes more.
    windows.sort(key=lambda ov: (ov[0] * ov[1], -abs(ov[0] - ov[1])), reverse=True)
    return windows


def compute_cis(method_keyword, molecules):
    """One MOPAC run per molecule to size the window, then a search for a window
    MOPAC will diagonalise in full.

    Two things make this more than a single call.

    First, the window cannot be chosen before the orbital count is known, and
    asking for more occupied or virtual orbitals than exist makes MOPAC quietly
    pick a different one -- so the window is read back out of the output rather
    than assumed.

    Second, and less obviously, **MOPAC's INDO CI does not always use the whole
    window**. Its configuration list is screened, and in some molecules most of
    it is thrown away: CCl4 asked for 4x4 reports `CI excitations= 1`, and
    benzene asked for 5x5 keeps 21 of 26. A CIS in this crate has no such
    screening, so a screened CI is a different calculation wearing the same
    name. The search takes the largest window whose configuration count is
    exactly `occupied * virtual + 1`, which is the largest window on which the
    two sides are provably computing the same thing.
    """
    rows = []
    version = None
    for name, charge, mult, atoms, why in molecules:
        xyz = os.path.join(OUT_DIR, "_scratch.xyz")
        os.makedirs(OUT_DIR, exist_ok=True)
        with open(xyz, "w") as f:
            f.write(gb.to_xyz_text(atoms, name))
        chosen = None
        try:
            ground = run_mopac.run(xyz, method=method_keyword, charge=charge,
                                   mult=mult, mode="1scf")
            n_occ = ground.get("n_occupied")
            n_mo = len(ground.get("eigenvalues_ev") or [])
            if not n_occ or not n_mo:
                print("  SKIP {}: no orbital count".format(name), file=sys.stderr)
                continue
            tried = []
            for occ, vir in cis_window_candidates(n_occ, n_mo - n_occ):
                result = run_mopac.run(
                    xyz, method=method_keyword, charge=charge, mult=mult, mode="1scf",
                    extra_keywords=("CIS", "C.I.=({},{})".format(occ + vir, occ),
                                    "WRTCI={}".format(occ * vir + 1)),
                )
                window = result.get("ci_window")
                count = result.get("ci_count")
                if not window or not result.get("ci_states"):
                    tried.append("{}x{}:no table".format(occ, vir))
                    continue
                used_occ = window[1] - window[0] + 1
                used_vir = window[3] - window[2] + 1
                if (used_occ, used_vir) != (occ, vir):
                    tried.append("{}x{}:became {}x{}".format(occ, vir, used_occ, used_vir))
                    continue
                if count != occ * vir + 1:
                    tried.append("{}x{}:{} of {}".format(occ, vir, count, occ * vir + 1))
                    continue
                chosen = (occ, vir, result)
                break
        except Exception as exc:  # noqa: BLE001 - report and skip, do not fabricate
            print("  SKIP {}: {}".format(name, exc), file=sys.stderr)
            continue
        finally:
            if os.path.exists(xyz):
                os.remove(xyz)

        if chosen is None:
            print("  SKIP {}: no window MOPAC diagonalises in full ({})".format(
                name, "; ".join(tried[:6])), file=sys.stderr)
            continue
        occ, vir, result = chosen
        version = version or result.get("mopac_version")
        keep = result["ci_states"][:CIS_ROOTS]
        rows.append({
            "name": name,
            "charge": charge,
            "multiplicity": mult,
            "n_atoms": len(atoms),
            "geometry": gb.to_xyz_field(atoms),
            "active_occupied": occ,
            "active_virtual": vir,
            "excitation_energies_ev": [st["energy_ev"] for st in keep],
            "oscillator_strengths": [st["oscillator_strength"] for st in keep],
            "state_dipoles_debye": [st["dipole_magnitude_debye"] for st in keep],
            "why": why,
        })
    return rows, version


def compute_molds(theory, molecules):
    """Run a locally built MolDS over a chosen molecule set (Tier 2).

    Tier 1 -- mining the regression outputs MolDS ships -- came first and is
    still what validates the build: `third_party/molds/tests/` is reproduced
    byte for byte by the binary these numbers come from. What Tier 1 could not
    do is choose the molecules, and a suite resting on the three cases upstream
    happened to test is narrow however many scalars each one yields. This set is
    drawn from the same table the MOPAC suites use, so a molecule in both is the
    same geometry in both.
    """
    rows = []
    for name, charge, mult, atoms, why in molecules:
        if charge != 0:
            # MolDS's input deck has no total-charge keyword; every calculation
            # is neutral. Saying so is better than silently dropping the row.
            print("  SKIP {}: MolDS takes no molecular charge".format(name), file=sys.stderr)
            continue
        try:
            parsed = run_molds.run(atoms, theory)
        except Exception as exc:  # noqa: BLE001 - report and skip, do not fabricate
            print("  SKIP {}: {}".format(name, exc), file=sys.stderr)
            continue
        dip = parsed["dipole_debye"] or [None, None, None]
        rows.append({
            "name": name,
            "charge": charge,
            "multiplicity": mult,
            "n_atoms": len(atoms),
            # The geometry MolDS echoed back, at the precision it echoed it, so
            # the two programs are given identical coordinates (item 23).
            "geometry": gb.to_xyz_field(parsed["atoms"]),
            "total_energy_ev": parsed["total_energy_ev"],
            "core_repulsion_ev": parsed["core_repulsion_ev"],
            "charges": parsed["charges"],
            "eigenvalues_ev": parsed["eigenvalues_ev"],
            "dipole_x_debye": dip[0], "dipole_y_debye": dip[1], "dipole_z_debye": dip[2],
            "why": why,
        })
    return rows, "0.3.1"


MINDO3_COLUMNS = [
    "name", "charge", "multiplicity", "n_atoms", "geometry",
    "hof_kcal", "total_energy_ev", "electronic_energy_ev", "core_repulsion_ev",
    "ionization_potential_ev", "homo_ev", "lumo_ev",
    "dipole_magnitude_debye",
    "charges", "eigenvalues_ev", "why",
]


def uncollinear_first_three(atoms):
    """Reorder so the first three atoms are not in a straight line, or `None`.

    MOPAC7 rejects Cartesian input whose first three atoms are collinear
    (`getgeo.f:373`). The molecule is not the problem and neither is the model:
    it is which three atoms happen to be written first. Moving one non-collinear
    atom to the third position is enough, and changes nothing physical -- the
    same coordinates in a different order, and the reordered list is what goes
    into the reference file, so both programs are given the same atom order and
    the per-atom charges line up.

    `None` when every atom is collinear with the first two, which is a genuinely
    linear molecule (ethyne) that no ordering can help. The caller declares
    those rather than working around them.
    """
    if len(atoms) < 4 or not run_mopac7.first_three_are_collinear(atoms):
        return atoms
    for k in range(3, len(atoms)):
        candidate = list(atoms)
        candidate[2], candidate[k] = candidate[k], candidate[2]
        if not run_mopac7.first_three_are_collinear(candidate):
            return candidate
    return None


def compute_mopac7(molecules):
    """Run MOPAC7 1.15 over the MINDO/3 set.

    MOPAC7 prints the total energy, the electronic energy and the core-core
    repulsion as three separate numbers, which MOPAC 23 does not. That is worth
    having: it says whether a disagreement lives in the density or in the
    geometry term, which is the distinction that took the CNDO/2 and INDO suites
    apart (item 24).

    `n_occupied` is read from MOPAC7's own "NO. OF FILLED LEVELS", so the HOMO
    and LUMO are indexed by the oracle's count rather than by one recomputed
    here from the electron count -- the two agreeing is part of what is being
    checked, so deriving one from the other would check nothing.
    """
    rows = []
    for name, charge, mult, atoms, why in molecules:
        ordered = uncollinear_first_three(atoms)
        if ordered is None:
            print("  SKIP {}: every atom is collinear, and MOPAC7 rejects "
                  "Cartesian input whose first three atoms are "
                  "(getgeo.f:373)".format(name), file=sys.stderr)
            continue
        atoms = ordered
        try:
            parsed = run_mopac7.run(
                atoms, run_mopac7.keywords_for(charge, mult, atoms=atoms))
        except Exception as exc:  # noqa: BLE001 - report and skip, do not fabricate
            print("  SKIP {}: {}".format(name, str(exc)[:120]), file=sys.stderr)
            continue
        eig = parsed.get("eigenvalues_ev") or []
        homo = lumo = None
        n_occ = parsed.get("n_occupied")
        if eig and n_occ:
            homo = eig[n_occ - 1] if 0 < n_occ <= len(eig) else None
            lumo = eig[n_occ] if n_occ < len(eig) else None
        rows.append({
            "name": name,
            "charge": charge,
            "multiplicity": mult,
            "n_atoms": len(atoms),
            "geometry": gb.to_xyz_field(atoms),
            "hof_kcal": parsed.get("heat_of_formation_kcal"),
            "total_energy_ev": parsed.get("total_energy_ev"),
            "electronic_energy_ev": parsed.get("electronic_energy_ev"),
            "core_repulsion_ev": parsed.get("core_repulsion_ev"),
            "ionization_potential_ev": parsed.get("ionization_potential_ev"),
            "homo_ev": homo,
            "lumo_ev": lumo,
            # The magnitude, not the components: a molecule of three atoms or
            # fewer goes to MOPAC7 as a Z-matrix, so the components are in
            # MOPAC7's reconstructed frame. See `run_mopac7.write_input`.
            "dipole_magnitude_debye": parsed.get("dipole_magnitude_debye"),
            "charges": parsed.get("charges"),
            "eigenvalues_ev": eig,
            "why": why,
        })
    return rows, "1.15"


def cell(value):
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ",".join("" if v is None else repr(float(v)) for v in value)
    if isinstance(value, float):
        return repr(value)
    return str(value)


def render(rows, method_keyword, oracle_version, generator_sha, columns=None, keywords=None,
           program="OpenMOPAC"):
    columns = columns or COLUMNS
    out = []
    out.append("# PROVENANCE: values are what OpenMOPAC printed; geometries are what was sent to it.")
    out.append("# oracle_program\t{}".format(program))
    out.append("# oracle_version\t{}".format(oracle_version))
    out.append("# oracle_keywords\t{}".format(
        keywords or "{} PRECISE {} 1SCF AUX(PRECISION=9)".format(
            method_keyword, run_mopac.SCF_TIGHTENING)))
    out.append("# generator\ttools/oracle/build_reference_set.py")
    out.append("# generator_sha256\t{}".format(generator_sha))
    out.append("# generated_utc\t{}".format(
        datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")))
    out.append("# geometry: atoms separated by '|', fields within an atom by spaces, Angstrom.")
    out.append("# An empty cell means the oracle did not report that quantity; the test skips it.")
    out.append("\t".join(columns))
    for row in rows:
        out.append("\t".join(cell(row[c]) for c in columns))
    return "\n".join(out) + "\n"


def report_coverage(rows):
    elements = {}
    charges = {}
    mults = {}
    for row in rows:
        for atom_field in row["geometry"].split("|"):
            elements[atom_field.split()[0]] = elements.get(atom_field.split()[0], 0) + 1
        charges[row["charge"]] = charges.get(row["charge"], 0) + 1
        mults[row["multiplicity"]] = mults.get(row["multiplicity"], 0) + 1
    print("  molecules: {}".format(len(rows)))
    print("  elements ({}): {}".format(
        len(elements), " ".join("{}:{}".format(k, v) for k, v in sorted(elements.items()))))
    print("  charges: {}".format(dict(sorted(charges.items()))))
    print("  multiplicities: {}".format(dict(sorted(mults.items()))))


STALE, UP_TO_DATE, UNCHECKED, INCOMPLETE = "stale", "ok", "unchecked", "incomplete"


def _strip_stamp(text):
    # The generated_utc line legitimately differs between runs.
    return "\n".join(l for l in text.splitlines() if not l.startswith("# generated_utc"))


def data_row_count(text):
    """Molecule rows in a rendered reference file: no comments, no column header."""
    body = [line for line in text.splitlines() if line and not line.startswith("#")]
    return max(len(body) - 1, 0)


def emit(path, text, rows, check, program):
    """Write the rendered file, or compare it against the committed one.

    ``--check`` re-runs the oracle, so it says something only when the oracle
    ran. None of the three binaries is bundled, so the first thing a reader of
    the distributed source meets is every molecule skipped and a file rendered
    from no rows, which differs from the committed one in every line. Calling
    that STALE is a false alarm that reads as "the shipped data is wrong": it
    is not known to be wrong, it is unchecked. Report the two separately, and
    fail on the second only.

    The staleness check that works without the oracle is the hash manifest in
    tests/data/MANIFEST.sha256, which pins the reference files *and* this
    generator, and which tests/oracle_matrix.rs runs in every checkout.
    """
    name = os.path.basename(path)
    if not check:
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
        print("  wrote {} ({} rows)".format(name, len(rows)))
        return UP_TO_DATE

    existing = open(path, encoding="utf-8").read() if os.path.exists(path) else ""
    committed = data_row_count(existing)
    if not rows:
        print("  CANNOT CHECK: {}: {} produced no rows, so there is nothing to "
              "compare the {} committed rows against".format(name, program, committed),
              file=sys.stderr)
        return UNCHECKED
    if len(rows) < committed:
        print("  INCOMPLETE: {}: {} produced {} of the {} committed rows; repair the "
              "oracle environment rather than trusting this answer".format(
                  name, program, len(rows), committed), file=sys.stderr)
        return INCOMPLETE
    if _strip_stamp(existing) != _strip_stamp(text):
        print("  STALE: {}".format(name), file=sys.stderr)
        return STALE
    print("  up to date: {}".format(name))
    return UP_TO_DATE


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--method", default="all", choices=["all"] + sorted(SETS))
    ap.add_argument("--check", action="store_true",
                    help="regenerate and fail if the committed file differs")
    ap.add_argument("--report", action="store_true",
                    help="print the coverage histogram without running the oracle")
    args = ap.parse_args()

    with open(__file__, "rb") as f:
        generator_sha = hashlib.sha256(f.read()).hexdigest()

    names = sorted(SETS) if args.method == "all" else [args.method]
    verdicts = {}
    for key in names:
        keyword, filename, builder = SETS[key]
        molecules = builder()
        print("{}: {} molecules".format(key, len(molecules)))
        if args.report:
            rows = [{"geometry": gb.to_xyz_field(a), "charge": c, "multiplicity": m}
                    for _, c, m, a, _ in molecules]
            report_coverage(rows)
            continue
        if key == "mindo3":
            rows, version = compute_mopac7(molecules)
            report_coverage(rows)
            text = render(
                rows, keyword, version, generator_sha,
                columns=MINDO3_COLUMNS,
                keywords="{} (plus UHF <multiplicity> / CHARGE=q as needed)".format(
                    run_mopac7.keywords_for()),
                program="MOPAC7",
            )
            verdicts[filename] = emit(
                os.path.join(OUT_DIR, filename), text, rows, args.check, "MOPAC7")
            continue
        molds = key in ("cndo2", "indo")
        if molds:
            rows, version = compute_molds(keyword, molecules)
            report_coverage(rows)
            text = render(
                rows, keyword, version, generator_sha,
                columns=MOLDS_COLUMNS,
                keywords="{} rms_density=1e-9 (locally built MolDS 0.3.1)".format(keyword),
                program="MolDS",
            )
            verdicts[filename] = emit(
                os.path.join(OUT_DIR, filename), text, rows, args.check, "MolDS")
            continue
        cis = key.endswith("_cis")
        if cis:
            rows, version = compute_cis(keyword, molecules)
        else:
            rows, version = compute(keyword, molecules)
        report_coverage(rows)
        if version is None:
            # Under --check with no rows, MOPAC simply is not installed. emit()
            # reports that as unchecked; aborting here instead would turn a
            # missing oracle into a hard failure of the whole run.
            if not (args.check and not rows):
                raise SystemExit("no MOPAC version reported; refusing to stamp the file")
            version = "unknown"
        if cis:
            text = render(
                rows, keyword, version, generator_sha,
                columns=CIS_COLUMNS,
                keywords="{} PRECISE {} 1SCF CIS C.I.=(occ+vir,occ) WRTCI".format(
                    keyword, run_mopac.SCF_TIGHTENING),
            )
        else:
            text = render(rows, keyword, version, generator_sha)
        verdicts[filename] = emit(
            os.path.join(OUT_DIR, filename), text, rows, args.check, "OpenMOPAC")

    if not args.check:
        return 0
    unchecked = sorted(n for n, v in verdicts.items() if v == UNCHECKED)
    if unchecked:
        print("{} of {} reference files were not checked because their oracle is not "
              "installed: {}".format(len(unchecked), len(verdicts), ", ".join(unchecked)),
              file=sys.stderr)
        print("That is not a staleness verdict. tests/data/MANIFEST.sha256, checked by "
              "tests/oracle_matrix.rs, detects edits to these files without any oracle.",
              file=sys.stderr)
    if any(v in (STALE, INCOMPLETE) for v in verdicts.values()):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
