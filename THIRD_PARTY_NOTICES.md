# Third-party notices

`xndo-rs` is distributed under GPL-3.0-or-later. This file is the file-by-file
record of what it takes from other projects.

The retained upstream copyright notices are in [`NOTICE`](NOTICE). Licence texts
are in [`LICENSES/`](LICENSES). Rust crate licences, copied verbatim from what
each crate ships, are in
[`THIRD_PARTY_LICENSES_RUST.md`](THIRD_PARTY_LICENSES_RUST.md).

> **Correction.** Through v0.2.4 this file said it recorded "parameter-data
> sources" and that no upstream source was redistributed. That was wrong, and
> materially so: roughly a dozen Rust modules are **line-level ports of OpenMOPAC
> Fortran**, which their own doc comments have always said. They are derivative
> works of Apache-2.0 code and are listed as such below. The corresponding
> attribution and modification notices, which Apache-2.0 sections 4(b) and 4(c)
> and GPL-3.0 section 5(a) require, are now carried in `NOTICE` and in a
> `PROVENANCE:` header on each file.

## OpenMOPAC 23.2.5 (Apache-2.0)

Upstream: https://github.com/openmopac/mopac, tag `v23.2.5`, commit
`052691223d19935a89f0fe18cd12301bd83e4201`.
Licence: `third_party/mopac/LICENSE` (byte-identical to upstream) and
`LICENSES/Apache-2.0.txt`.
Provenance detail, including the evidence that upstream ships no NOTICE file:
`third_party/mopac/PROVENANCE.md`.
Retained attribution: `third_party/mopac/AUTHORS.rst`,
`third_party/mopac/CITATION.cff`.

### Ported code

Each file carries a `PROVENANCE:` header naming the upstream source and the
modification.

| xndo-rs file | derived from |
| --- | --- |
| `src/constants.rs` | `src/conref_C.F90` |
| `src/data_tables.rs` | `src/models/parameters_C.F90` and the per-method parameter modules |
| `src/integrals.rs` | `src/integrals/mndod.F90` (`reppd`, `spcore`) |
| `src/integrals_d.rs` | `src/integrals/mndod.F90`, `src/integrals/mndod_C.F90` |
| `src/onecenter.rs` | `src/integrals/mndod.F90` (`rsc`, `scprm`, `inighd`, `wstore`, `eiscor`), `src/integrals/mndod_C.F90` |
| `src/overlap.rs` | `src/integrals/diat.F90`, `src/integrals/set.F90` |
| `src/params.rs` | `src/models/calpar.F90`, `src/integrals/mndod.F90` (`inid`, `ddpo`, `poij`, `aijm`), `src/input/readmo.F90` |
| `src/repulsion.rs` | `src/integrals/ccrep.F90` |
| `src/rotations.rs` | `src/integrals/mndod.F90` (`rotmat`) |
| `src/scf.rs` | `src/SCF`, `src/models/parameters_C.F90` (shell occupancies) |
| `src/zindo.rs` | `src/INDO/scf.F90`, `integral.F90`, `ovlp.F90`, `ci.F90`, `reimers_C.F90`, `parameters_for_INDO_C.F90` |

The `src/INDO/` lineage is recorded in `NOTICE`: that code is by Rebecca
Gieseking, derived from Jeff Reimers' CNDO/INDO program (DOI:10.4231/D3R49G96G),
itself derived from the CNDO/S program of Del Bene, Jaffe, Ellis and Kuehnlenz
(QCPE 174).

### Extracted parameter data

All are embedded with `include_str!` (`src/data_tables.rs`) and hashed in
`src/data/legacy_parameter_manifest.sha256`.

- `src/data/element_data.csv` (from `src/models/parameters_C.F90`: `ios`, `iop`,
  `iod`, `npq`, `main_group`, `ndelec`, `eheat`, `ams`)
- `src/data/mndo_parameters.csv`
- `src/data/mndo_pair_parameters.csv`
- `src/data/mndod_parameters.csv`
- `src/data/mndod_pair_parameters.csv`
- `src/data/zindo_s_parameters.csv`

> `element_data.csv` was embedded but absent from this list through v0.2.4.
> `tests/attribution.rs` now asserts that every `include_str!`-ed data file
> appears here, so the set cannot silently narrow again.

### Reference values used in tests

`src/gradient.rs` and `src/hessian.rs` contain gradients and frequencies
computed by OpenMOPAC and used as test oracles, with the keywords and geometry
recorded beside them. `tests/data/*.tsv` holds the same kind of thing at larger
scale, described in `tests/data/ORACLE_NOTES.md`. These are measurements
produced by running the program, not copied code.

## MolDS 0.3.1 (GPL-3.0-or-later)

Upstream: https://sources.debian.org/src/molds/0.3.1-3/
Original project: https://osdn.net/projects/molds/
Copyright: (C) 2011-2012 Mikiya Fujii; (C) 2012-2013 Katsuhiko Nishimra.
Licence: `LICENSES/GPL-3.0-or-later.txt`. Detail: `third_party/molds/NOTICE`.

### Ported code

| xndo-rs file | derived from |
| --- | --- |
| `src/cndo_indo.rs` | `src/cndo/Cndo2.cpp`, `src/indo/Indo.cpp` |

### Extracted parameter data

- `src/data/molds_cndo2_indo_parameters.csv`
- `src/data/molds_zindo_s_parameters.csv`

### Retained upstream test data

`third_party/molds/tests/` holds five of the regression `.in`/`.dat` pairs MolDS
ships, verbatim. They are used to verify that a locally built MolDS reproduces
its author's own numbers before it is trusted as an oracle; all five come back
byte-identical. See `tests/data/ORACLE_NOTES.md` item 24.

## MOPAC7 1.15 (public domain)

Upstream: https://sources.debian.org/src/mopac7/1.15-4/fortran/
Detail: `third_party/mopac7/NOTICE`.

### Ported code

| xndo-rs file | derived from |
| --- | --- |
| `src/mindo3.rs` | `fortran/block.f` (MINDO/3 tables), `analyt.f`, `delri.f`, `calpar.f`, `compfg.f` |

### Extracted parameter data

- `src/data/mindo3_parameters.csv` (no deviation from upstream; the one this
  file used to declare, in `hp2_ev`, was fixed in v0.3.0 -- see
  `tests/data/ORACLE_NOTES.md` item 4)
- `src/data/mindo3_pair_parameters.csv`

## STO-nG expansions (`src/sto.rs`)

Added in v0.3.0 for the Molden export, which needs a Gaussian description of a
Slater basis.

### Numerical tables

The expansion coefficients are R. F. Stewart, "Small Gaussian Expansions of
Slater-Type Orbitals", *J. Chem. Phys.* **52**, 431-438 (1970): published
least-squares fits, reproduced as **cited numerical facts**. A table of fitted
constants in a subscription journal carries no licence to propagate, and
declaring one would pass an invented obligation to every downstream recipient.

### Ported code

| file | derived from | licence |
| --- | --- | --- |
| `src/sto.rs` | gfn2-rs `src/sto.rs` (same authors) | GPL-3.0-or-later |
| `src/sto.rs` (table layout, index arithmetic, `slater_to_gauss` normalisation) | xtb `slaterToGauss` | LGPL-3.0-or-later |

The second row is the one that is easy to leave out. The *numbers* are published
facts, but their arrangement -- the fifteen-row layout, the `n`/`l` index
arithmetic, the separate STO-6G tables for 6s and 6p, the `zeta^2` scaling and
Cartesian normalisation -- follows xtb, which is LGPL-3.0-or-later. That is
compatible with this crate's GPL-3.0-or-later, and `LICENSES/LGPL-3.0-or-later.txt`
is carried on that basis.

`src/gto.rs` and `src/molden.rs`, which use it, are this project's own work.

## Files this project generates about the above

Embedded with `include_str!` and hashed like the rest, but xndo-rs's own work
rather than anything taken from upstream:

- `src/data/legacy_parameter_catalog.csv` -- which bundled dataset serves which
  method, and under what licence
- `src/data/legacy_parameter_manifest.sha256` -- the hashes, checked by
  `tools/validate_legacy_parameters.py` and by `tests/attribution.rs`

## Rust crates

Roughly ninety crates are linked into the library and binaries. Their licences
are reproduced in `THIRD_PARTY_LICENSES_RUST.md`, generated by
`tools/collect_rust_licenses.py` from the files each crate ships rather than
from its declared SPDX expression.

`faer` is why that distinction is not academic: it declares `license = "MIT"`
but ships `COPYING.EIGEN.MPL2`, `COPYING.LAPACK.BSD`,
`COPYING.SUITE_SPARSE.AMD.BSD` and `COPYING.SUITE_SPARSE.COLAMD.BSD`, because
parts of it are ported from Eigen, LAPACK and SuiteSparse. It is a direct
dependency present in every artifact, so the MPL-2.0 obligation is real. The
text is in `LICENSES/MPL-2.0.txt`.

## Removed material

Version 0.2.4 removed D3, H4, HX, PM3, Sparkle and MSINDO-sCIS data. Notices
that applied only to removed files were removed with them.
