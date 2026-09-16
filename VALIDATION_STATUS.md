# Validation status for xndo-rs 0.3.0

This record separates executable checks, finite-difference validation,
independent-oracle evidence, and scientific scope boundaries.

## Build and API validation

The release was checked on Windows with one Cargo build job:

```text
cargo check --all-targets --all-features
cargo clippy --all-targets --all-features -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --all-features --no-deps
cargo test --all-targets --all-features
cargo test --test derivative_fd_matrix --all-features -- --ignored --nocapture --test-threads=1
maturin build --release --out dist
maturin sdist --out dist
python -m pytest -q tests/test_python.py tests/test_python_api.py
python -m twine check dist/*
python tools/validate_legacy_parameters.py
python tools/validate_cndo_indo_wiring.py
python tools/validate_method_aliases.py
```

Since 0.3.0 the release also runs:

```text
python tools/privacy_scan.py --tree .
python tools/collect_rust_licenses.py --check
python tools/oracle/build_reference_set.py --check    # needs the oracle binaries
python tools/study/render_report.py --check
```

Results:

- 195 ordinary Rust tests passed; 0 failed;
- both ignored full finite-difference release tests passed explicitly;
- 30 tests passed against the installed CPython 3.9+ ABI3 wheel, including
  Python-native and ASE RHF, UHF, RHF-CIS, UCIS, and six-method optimization
  workflows;
- strict Clippy and rustdoc both passed with `-D warnings`; this includes the
  61 Clippy findings and seven private-item documentation links inherited from
  the input source;
- the Rust, Python-native, and ASE surfaces executed at runtime; CLI
  optimization converged for CNDO/2, INDO, MNDO, MNDO/d, MINDO/3, and ZINDO/S
  while accepting custom iteration, gradient-tolerance, and history controls;
- all parameter SHA-256, sentinel, row-count, wiring, and alias audits passed;
- wheel and single-root source distribution passed `twine check`.

ASE 3.29.0 emitted six `DeprecationWarning` messages from its own array-shape
assignment when used with NumPy 2.5.1. They are not Rust compiler or xndo-rs
runtime warnings.

## Analytic derivative validation

All production gradients were compared component by component with central
differences of fully reconverged energies. All production Hessians were
compared element by element with central differences of analytic gradients;
matrix symmetry and translational behavior were checked separately.

Representative maximum absolute differences:

| Path | Maximum difference |
| --- | ---: |
| CNDO/2, INDO, MINDO/3, ZINDO/S RHF/UHF ground gradients | 1.69e-6 eV/Bohr |
| CNDO/2, INDO, MINDO/3, ZINDO/S RHF/UHF ground Hessians | 3.12e-6 eV/Bohr^2 |
| MNDO/MNDO-d RHF/UHF release-matrix gradients | 4.05e-6 eV/Bohr |
| MNDO/MNDO-d RHF/UHF release-matrix Hessians | 2.96e-5 eV/Bohr^2 |
| MNDO/d UHF SH d-AO gradient | 7.50e-7 eV/Bohr |
| MNDO/d UHF SH d-AO analytic Hessian | 4.20e-6 eV/Bohr^2 |
| ZINDO/S singlet RHF-CIS state gradient | 4.53e-6 eV/Bohr |
| ZINDO/S triplet RHF-CIS state gradient | 1.23e-6 eV/Bohr |
| ZINDO/S singlet RHF-CIS state Hessian | 7.85e-5 eV/Bohr^2 |
| ZINDO/S triplet RHF-CIS state Hessian | 6.11e-6 eV/Bohr^2 |
| ZINDO/S UHF-UCIS state gradient | 1.62e-6 eV/Bohr |
| ZINDO/S UHF-UCIS state Hessian | 2.04e-5 eV/Bohr^2 |

The MNDO/MNDO-d matrix covers H2O, NH3, H2CO, CH3, OH, representative
transition-metal d-shell cases, every explicit main-group MNDO/d d-bearing
element, and an open-shell SH d-AO case. The ZINDO/S UHF and UCIS matrix covers
distorted CH3, NH2, and HO2 radicals.

OH is retained. Its ZINDO/S UHF ground gradient agrees with reconverged energy
differences to 1.01e-7 eV/Bohr, and its vertical UCIS spectrum is tested.
State-specific UHF Hessians and UCIS derivatives for an open-shell diatomic are
not unique on the degenerate electronic surface without a state-averaged or
diabatic definition, so those APIs return an explicit error. State-specific
CIS/UCIS derivatives likewise require root separation of at least 0.0001 eV;
vertical spectra remain available for closer roots.

No empirical rescaling or relaxed tolerance is used to conceal an SCF branch
change, orbital degeneracy, or root crossing.

These numbers are unchanged from 0.2.4, and that is itself worth recording.
0.3.0 fixed a defect in which all three coupled-perturbed solvers returned their
last iterate after running out of iterations without saying so, which means
every Hessian in every release up to 0.2.4 was built from a response that had
not met its own 1e-9 tolerance -- water reaches only 1.88e-8 in the 100 cycles
it was allowed. The table above did not move, because 2e-8 is well under the
1e-5 scale these comparisons resolve. The agreement was real; the *guarantee*
was not, and now is.

## Independent oracle evidence

This is the part of this document that changed most in 0.3.0. Through 0.2.4 the
record said, correctly, that no oracle binary had been run: the evidence was
frozen values and finite differences. It now rests on three oracle programs,
each one the program the method's own parameters came from.

| method | oracle | molecules | elements | independent points |
| --- | --- | ---: | ---: | ---: |
| MNDO | OpenMOPAC 23.2.5 | 51 | 19 | 309 |
| MNDO/d | OpenMOPAC 23.2.5 | 37 | 17 | 230 |
| ZINDO/S | OpenMOPAC 23.2.5 `INDO` | 41 | 10 | 465 |
| ZINDO/S CIS | OpenMOPAC 23.2.5 `INDO CIS` | 24 | 10 | 322 |
| CNDO/2 | MolDS 0.3.1 | 17 | 6 | 270 |
| INDO | MolDS 0.3.1 | 14 | 5 | 208 |
| MINDO/3 | MOPAC7 1.15 | 36 | 10 | 611 |

2415 independent scalar comparisons. An "independent point" is one scalar the
oracle printed and this crate produced separately, with derived quantities and
symmetry-forced zeros excluded and linear constraints (charges summing to the
molecular charge) discounted. The harness computes the count and a meta-test
asserts every method clears fifty.

Worst observed deviation per quantity, per suite, is printed on every run and
committed beside each tolerance floor. No floor may be raised without the
measurement that justifies it in the same commit.

**The list of known disagreements is empty.** Every deviation the suites found
was either fixed at the root or shown to be the oracle's own arithmetic:

- MNDO germane (1.209e-3 kcal/mol) and MNDO/d H2S (9.014e-4) sit at the branch
  point of MOPAC's truncated B-integral series, confirmed by locating the
  predicted crossing and watching the disagreement step across it;
- CNDO/2 and INDO energies floor at about 4e-4 eV on MolDS's seven-figure print
  format;
- MINDO/3 floors at MOPAC7's four-decimal charges and eigenvalues.

Two things are deliberately **not** compared, and saying so is part of the
claim: CNDO/2 and INDO carry no ions, because MolDS's input deck has no
total-charge keyword; and their dipole is not compared at all, because MolDS and
this crate use two different definitions of it in a zero-differential-overlap
model, which differ by 0.23 D on CH4.

`tests/data/ORACLE_NOTES.md` records all of this at length, including every
convention and trap found along the way.

## Parameters and licenses

Bundled reusable parameter data are:

- OpenMOPAC 23.2.5 under Apache-2.0 for MNDO, MNDO/d, and ZINDO/S;
- MolDS 0.3.1 under GPL-3.0-or-later for CNDO/2, INDO, and a ZINDO/S
  cross-check subset;
- public-domain MOPAC7 for the consolidated MINDO/3 tables and sentinels.

The missing-parameter search was repeated for every registered but gated
method. No additional candidate provided a complete numerical table, an
unambiguous matching Hamiltonian convention, and a compatible redistribution
license together. No nearby-method values were substituted, and rejected
source-code names are intentionally omitted.

## Performance and resource behavior

Embedded data tables are parsed once per process. RHF and UHF ZDO response
engines assemble all Cartesian right-hand sides and solve them from one LU
factorization. Second-order UHF response uses paired alpha/beta DIIS. The
MNDO-family Hessian streams response chunks under the configured memory budget;
the open-shell d-AO path now uses analytic spd automatic differentiation and
spin-coupled UCPHF instead of a numerical-Hessian fallback.
