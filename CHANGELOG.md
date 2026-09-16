# Changelog

## 0.3.0

`docs/v0.3.0-delivery.md` is the full record, including what was planned and
dropped, and the known problems in this release.

### Licensing

- Brought the attribution record up to what Apache-2.0 4(b)/4(c) and GPL-3.0
  5(a) require. v0.2.4 retained no upstream copyright notice, stated no
  modification or date, described a dozen line-level Fortran ports as
  "parameter data", and embedded `src/data/element_data.csv` without
  attributing it anywhere.
- Added `NOTICE`, `LICENSES/` (ten texts), `THIRD_PARTY_LICENSES_RUST.md`,
  `src/licenses.rs`, `xndo_rs_cli licenses`, and `tests/attribution.rs`.
- Rust crate licences are collected from the files each crate ships rather than
  from its SPDX field, which is how `faer`'s MPL-2.0 obligation is caught: it
  declares `MIT` and ships `COPYING.EIGEN.MPL2`.
- Added `tools/privacy_scan.py`, which reads inside archives and searches
  binaries in UTF-16LE as well as UTF-8.

### Verification

- All six methods now have an enforced oracle suite against the program their
  parameters came from: 2415 independent points over 220 molecules, with no
  known disagreements. MINDO/3's, against MOPAC7 1.15, is new in this release.

### Fixed

- **`hp2_ev` was stored rounded** for N, Si, S and Cl. MOPAC7 forms that
  integral as `(gpp - gp2)/2` and stores nothing. Worth 2.2e-2 eV per chlorine
  atom; CCl4 was 8.6e-2 eV out.
- **MINDO/3's Slater exponents were not restated in MOPAC7's Bohr**, a 1.93e-5
  relative error in every overlap exponential and a few times 1e-4 eV per atom
  pair.
- **The SCF spin guess placed the excess spin in proportion to the density**,
  which is backwards and put it on the *most* occupied orbital. Planar methyl
  under MINDO/3 converged to a real SCF solution 3.07 eV above the right one.
  An open shell now runs from two starts and keeps the lower.
- **The CPHF solvers reported nothing on running out of iterations**, so every
  Hessian in every release up to v0.2.4 was built from a response that had not
  met its own tolerance -- water needs more than the 100 cycles it was given.
- **The ZINDO/S UHF path had no SCF accelerator.**
- **`NddoResult` discarded its beta orbitals**, so MNDO and MNDO/d reported the
  lowest unoccupied *alpha* orbital as the LUMO for open shells.

### Added

- Orbital energies, occupations and the HOMO/LUMO pair for every method, through
  Rust, Python, ASE and a CLI `orbitals` subcommand.
- Molden wavefunction output, with coefficients back-transformed by `S^(-1/2)`
  so that they are orthonormal over the Gaussians in the file. Through Rust,
  Python, ASE, a `molden` subcommand and a `--molden` flag.
- A photoexcitation study (`studies/photoexcitation/`) run as a test, measuring
  ZINDO/S-CIS against exactly known limits. It recovers the model's
  charge-transfer coefficient as 1.1988 against its `TOMK = 1.2`.
- `third_party_licenses()` in the Python API.

### Not in this release

No correlated method of any kind, and no implicit solvation. See
`docs/scope.md`.

## 0.2.4

- Replaced the MNDO-specific geometry optimizer with one unified L-BFGS
  `optimize` API for CNDO/2, INDO, MNDO, MNDO/d, MINDO/3, and ZINDO/S across
  Rust, Python, and CLI interfaces.
- Exposed adjustable `max_iter`, maximum-gradient `gtol`, and L-BFGS `history`
  controls in Python and CLI while retaining `OptOptions` in Rust.
- Added all-method optimization execution, capability-discovery, and invalid-
  control regression tests.
- Resolved all 61 strict Clippy findings without changing published scientific
  parameters, and resolved all seven private-item rustdoc links.
- Updated API guides, release records, package metadata, and examples to
  version 0.2.4.

## 0.2.3

- Made `XNDO` the documented generic ASE calculator name, selected with an
  explicit `method`; the method-specific `MNDO` convenience class remains
  available and exported, while `XNDO` is intentionally omitted from the ASE
  module's `__all__` because XNDO is the library name rather than a method.
- Expanded the Python/ASE API guide with complete signatures, accepted inputs,
  return dictionaries and units, calculator properties, and failure behavior.
- Expanded the Rust API guide to cover every crate-root function and type,
  result fields and units, method discovery, configuration, errors, and the
  advanced public module layout.
- Updated Rust, Python, notices, validation records, and discovery messages to
  version 0.2.3.

## 0.2.2

- Added analytic ground-state gradients and Hessians for CNDO/2, INDO, and
  MINDO/3 with both RHF and UHF references.
- Added analytic ZINDO/S RHF/UHF ground-state gradients and Hessians.
- Added state-specific analytic gradients and Hessians for spin-adapted
  singlet and triplet ZINDO/S CIS roots, including excitation-only and absolute
  state derivatives.
- Added spin-orbital ZINDO/S UCIS spectra, UV-visible properties, analytic
  state gradients, and analytic state Hessians on UHF references.
- Extended the MNDO/d analytic UHF Hessian through d AOs using spd automatic
  differentiation and spin-coupled UCPHF, removing the numerical fallback.
- Exposed the new derivatives and derived frequencies through Rust,
  Python-native, ASE, and CLI APIs.
- Added full-coordinate finite-difference energy/gradient checks for every new
  derivative path, including CH3, NH2, HO2, OH, and an open-shell d-AO SH case.
- Retained OH UHF gradients and vertical UCIS spectra while explicitly
  rejecting non-unique state-specific derivatives on its degenerate surface.
- Added a 0.0001 eV root-isolation guard so near-degenerate CIS/UCIS spectra
  remain available without reporting basis-dependent state derivatives.
- Reused one LU factorization for all Cartesian CPHF right-hand sides.
- Updated PyPI metadata, documentation, validation records, and API discovery.
- Resolved all compiler warnings for all targets and features.

## 0.2.1

- Removed the inherited D3, H4, and HX correction subsystem and all associated
  numerical data and notices.
- Removed PM3, Sparkle, capped-bond, point-charge, and PM3 oracle assets.
- Replaced PM3-named public NDDO types with `NddoCalculator`,
  `NddoParameters`, `NddoOptions`, `NddoResult`, and `XndoError`.
- Limited native NDDO dispatch to MNDO and MNDO/d.
- Added MolDS-derived Li parameters to CNDO/2 and ground-state INDO.
- Removed incomplete MSINDO literature transcriptions.
- Added per-process caches for embedded MNDO, MNDO/d, and ZINDO/S tables.
- Added explicit Rust, Python, ASE, and CLI UV-visible spectrum access for
  ZINDO/S, including wavelengths, oscillator/dipole strengths, state dipoles
  and charges, configurations, and hole/electron descriptors.
- Updated independent OpenMOPAC MNDO regression coverage, including UHF,
  gradient, optimization, and optimized-geometry frequency checks.
- Reworked documentation, provenance, API tests, and privacy checks.
- Updated package versions to 0.2.1 and resolved all Rust compiler warnings.

## 0.2.0

Initial generalized xndo-rs prototype.
