#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Parameterised molecular geometry builders for the oracle reference sets.

Why builders rather than typed coordinate blocks
------------------------------------------------

A reference set of a hundred molecules is a hundred coordinate blocks, and a
hundred hand-typed coordinate blocks cannot be reviewed: a transposed digit in
one of them is invisible. These builders express each molecule as a shape plus
the bond lengths and angles that define it, so what a reviewer checks is
"pyramidal, N-H 1.01 A, 107 degrees", which is a claim they can evaluate.

What is being tested is agreement between two programs *at the same geometry*.
The geometry therefore has to be reasonable and reproducible; it does not have
to be a minimum on either program's surface. Bond lengths here are close to
experiment.

Every builder returns ``[(symbol, x, y, z), ...]`` in Angstrom.

One trap, worth stating because it is easy to hit when extending these sets:
when you substitute a heavier atom into an existing molecule, substitute its
bond length too. Copying a row and changing only the symbol produces a geometry
that is not merely unrelaxed but physically wrong, and the resulting comparison
still "passes" because both programs are evaluated at the same wrong geometry.
That is why every builder takes its lengths explicitly instead of defaulting
them per element.
"""

from __future__ import annotations

import math

# Unit vectors towards the corners of a regular tetrahedron and octahedron,
# used by the builders below so that the shapes are defined once.
TETRAHEDRON = [
    (1.0, 1.0, 1.0),
    (1.0, -1.0, -1.0),
    (-1.0, 1.0, -1.0),
    (-1.0, -1.0, 1.0),
]
OCTAHEDRON = [
    (1.0, 0.0, 0.0),
    (-1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, -1.0, 0.0),
    (0.0, 0.0, 1.0),
    (0.0, 0.0, -1.0),
]


def _unit(v):
    n = math.sqrt(sum(c * c for c in v))
    return tuple(c / n for c in v)


def atom(symbol, x=0.0, y=0.0, z=0.0):
    return (symbol, float(x), float(y), float(z))


def linear(first, *rest):
    """``linear("O", ("C", 1.16), ("O", 1.16))`` - atoms strung along z.

    Each ``(symbol, length)`` is placed that far beyond the previous atom, so
    the arguments read as a chain of bond lengths rather than as absolute
    coordinates.
    """
    atoms = [atom(first)]
    z = 0.0
    for symbol, length in rest:
        z += length
        atoms.append(atom(symbol, 0.0, 0.0, z))
    return atoms


def bent(centre, a, b, r1, r2, angle_deg):
    """A-centre-B with the centre at the origin and the bisector along +z."""
    half = math.radians(angle_deg) / 2.0
    return [
        atom(centre),
        atom(a, r1 * math.sin(half), 0.0, r1 * math.cos(half)),
        atom(b, -r2 * math.sin(half), 0.0, r2 * math.cos(half)),
    ]


def pyramidal(centre, ligand, r, angle_deg, n=3):
    """C3v (n=3) or equivalent: `n` ligands at `angle_deg` from the axis.

    `angle_deg` is the ligand-centre-ligand angle, which is what tables report.
    """
    # Convert the ligand-centre-ligand angle to the angle from the C3 axis.
    cos_ll = math.cos(math.radians(angle_deg))
    cos_axis = math.sqrt(max(0.0, (1.0 + 2.0 * cos_ll) / 3.0)) if n == 3 else 0.0
    theta = math.acos(max(-1.0, min(1.0, cos_axis)))
    atoms = [atom(centre)]
    for i in range(n):
        phi = 2.0 * math.pi * i / n
        atoms.append(
            atom(
                ligand,
                r * math.sin(theta) * math.cos(phi),
                r * math.sin(theta) * math.sin(phi),
                r * math.cos(theta),
            )
        )
    return atoms


def trigonal_planar(centre, ligands, lengths):
    """Three ligands at 120 degrees in the xy plane."""
    atoms = [atom(centre)]
    for i, (symbol, r) in enumerate(zip(ligands, lengths)):
        phi = 2.0 * math.pi * i / 3.0
        atoms.append(atom(symbol, r * math.cos(phi), r * math.sin(phi), 0.0))
    return atoms


def tetrahedral(centre, ligand, r, substitutions=()):
    """Td centre with four ligands; ``substitutions`` replaces corners.

    ``substitutions`` is ``[(index, symbol, length), ...]`` so that a
    substituted corner always carries its own bond length - see the trap in the
    module docstring.
    """
    atoms = [atom(centre)]
    lengths = [r] * 4
    symbols = [ligand] * 4
    for index, symbol, length in substitutions:
        symbols[index] = symbol
        lengths[index] = length
    for corner, symbol, length in zip(TETRAHEDRON, symbols, lengths):
        u = _unit(corner)
        atoms.append(atom(symbol, length * u[0], length * u[1], length * u[2]))
    return atoms


def octahedral(centre, ligand, r, substitutions=()):
    """Oh centre with six ligands; ``substitutions`` as in `tetrahedral`."""
    atoms = [atom(centre)]
    lengths = [r] * 6
    symbols = [ligand] * 6
    for index, symbol, length in substitutions:
        symbols[index] = symbol
        lengths[index] = length
    for corner, symbol, length in zip(OCTAHEDRON, symbols, lengths):
        atoms.append(atom(symbol, length * corner[0], length * corner[1], length * corner[2]))
    return atoms


def bipyramidal(centre, equatorial, r_eq, axial, r_ax):
    """D3h: three equatorial ligands in the xy plane, two axial along z."""
    atoms = [atom(centre)]
    for i in range(3):
        phi = 2.0 * math.pi * i / 3.0
        atoms.append(atom(equatorial, r_eq * math.cos(phi), r_eq * math.sin(phi), 0.0))
    atoms.append(atom(axial, 0.0, 0.0, r_ax))
    atoms.append(atom(axial, 0.0, 0.0, -r_ax))
    return atoms


def ring(symbols, radius, substituent=None, r_sub=1.09):
    """A planar regular ring; ``substituent`` adds one radial atom per vertex.

    ``ring(["C"] * 6, 1.397, "H", 1.084)`` is benzene.
    """
    n = len(symbols)
    atoms = []
    for i, symbol in enumerate(symbols):
        phi = 2.0 * math.pi * i / n
        atoms.append(atom(symbol, radius * math.cos(phi), radius * math.sin(phi), 0.0))
    if substituent is not None:
        for i in range(n):
            phi = 2.0 * math.pi * i / n
            rr = radius + r_sub
            atoms.append(atom(substituent, rr * math.cos(phi), rr * math.sin(phi), 0.0))
    return atoms


def methyl_group(carbon_at, axis, r_ch=1.09, angle_deg=109.5):
    """Three hydrogens completing a methyl carbon whose fourth bond is `axis`.

    `axis` points from the carbon towards its heavy-atom partner.
    """
    u = _unit(axis)
    # Any vector not parallel to u, to seed the perpendicular frame.
    seed = (0.0, 0.0, 1.0) if abs(u[2]) < 0.9 else (1.0, 0.0, 0.0)
    p = _unit(
        (
            seed[1] * u[2] - seed[2] * u[1],
            seed[2] * u[0] - seed[0] * u[2],
            seed[0] * u[1] - seed[1] * u[0],
        )
    )
    q = (
        u[1] * p[2] - u[2] * p[1],
        u[2] * p[0] - u[0] * p[2],
        u[0] * p[1] - u[1] * p[0],
    )
    theta = math.radians(180.0 - angle_deg)
    atoms = []
    for i in range(3):
        phi = 2.0 * math.pi * i / 3.0
        d = [
            -u[k] * math.cos(theta) + math.sin(theta) * (p[k] * math.cos(phi) + q[k] * math.sin(phi))
            for k in range(3)
        ]
        atoms.append(
            atom(
                "H",
                carbon_at[0] + r_ch * d[0],
                carbon_at[1] + r_ch * d[1],
                carbon_at[2] + r_ch * d[2],
            )
        )
    return atoms


def water_cluster(n, spacing=3.0):
    """`n` water molecules on a cubic grid - the performance ladder, not chemistry.

    Used by the benchmark suite, where what matters is a reproducible system of
    a given size rather than a sensible structure.
    """
    side = math.ceil(n ** (1.0 / 3.0))
    atoms = []
    placed = 0
    for i in range(side):
        for j in range(side):
            for k in range(side):
                if placed >= n:
                    break
                ox, oy, oz = i * spacing, j * spacing, k * spacing
                atoms.append(atom("O", ox, oy, oz))
                atoms.append(atom("H", ox + 0.9584, oy, oz))
                atoms.append(atom("H", ox - 0.2400, oy + 0.9278, oz))
                placed += 1
    return atoms


def to_xyz_field(atoms):
    """Encode a geometry into the single TSV `geometry` column.

    Atoms are separated by ``|`` and the four fields within an atom by spaces,
    matching the encoding every other data file in this repository uses.

    **Eight decimals, matching ``to_xyz_text`` and ``run_mopac.write_mop``
    exactly.** This used to be six, which meant the reference file recorded a
    geometry up to 5e-7 A away from the one the oracle was actually given, so
    every comparison carried a geometry-mismatch error on top of whatever it was
    meant to be measuring. On methane that was the *entire* reported
    disagreement: at matched geometries MOPAC and xndo-rs agree to 5e-9
    kcal/mol, against the 3.2e-5 the suite had been reporting.
    """
    return "|".join(
        "{} {:.8f} {:.8f} {:.8f}".format(sym, x, y, z) for sym, x, y, z in atoms
    )


def to_xyz_text(atoms, comment=""):
    """Standard XYZ text, for handing a geometry to an external program."""
    lines = [str(len(atoms)), comment]
    lines += ["{:<2} {:>14.8f} {:>14.8f} {:>14.8f}".format(s, x, y, z) for s, x, y, z in atoms]
    return "\n".join(lines) + "\n"


def elements(atoms):
    return sorted({sym for sym, _, _, _ in atoms})
