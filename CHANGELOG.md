# Changelog

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
