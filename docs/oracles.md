# Independent oracle validation

Every method is checked against the program **its own parameter tables came
from**, rather than against a proxy or against another method. That is the whole
policy, and it is what decides which oracle is used where.

| method | oracle | why that one |
| --- | --- | --- |
| MNDO, MNDO/d, ZINDO/S | OpenMOPAC 23.2.5 | the `src/data/*.csv` provenance headers name it |
| CNDO/2, INDO | MolDS 0.3.1 | those two engines are ports of it |
| MINDO/3 | MOPAC7 1.15 | MOPAC 23.2.5 does not implement MINDO/3 at all; MOPAC7 is the declared upstream of the tables |

None of the three binaries is bundled. `tools/oracle/run_mopac.py`,
`run_molds.py` and `run_mopac7.py` drive them,
`tools/oracle/build_reference_set.py` writes the reference files, and
`tests/oracle_matrix.rs` compares against them without needing any of the
binaries present.

## What is compared

2415 independent scalar comparisons over 220 molecules. Per method:
heat of formation or total energy, core-core repulsion where the oracle reports
it separately, Mulliken charges, the frontier orbital pair or the whole orbital
spectrum, dipole components, and for ZINDO/S CIS the excitation energy,
oscillator strength and excited-state dipole magnitude of each root.

An "independent point" is one scalar the oracle printed and this crate produced
separately. Derived quantities are not double-counted, values that are zero by
symmetry in both programs are compared but not counted, and linear constraints
are discounted -- N charges that must sum to the molecular charge count as N-1.
The harness computes the count, and a meta-test asserts every method clears
fifty from at least twelve molecules.

`VALIDATION_STATUS.md` has the per-method table.
`tests/data/ORACLE_NOTES.md` has twenty-six numbered items covering every
convention, oracle limitation and defect the comparison turned up.

## How disagreements are handled

Three categories, and each has one permitted response:

1. **a convention difference** (units, sign, frame, print precision, active-space
   definition) -- fix the comparison. The tolerance is never touched.
2. **a different SCF solution** -- report which is variationally lower and record
   it with an observable that distinguishes the two.
3. **a model difference** -- it is a bug or a scope decision, and it is resolved
   as one or the other.

Tolerance floors are committed with the measured worst case beside them, and a
floor may only be raised in the same commit that adds the measurement showing
why. Empirical scale factors are never applied, and a different method's values
are never substituted to force agreement.

The list of known disagreements is currently **empty**. Two residuals remain
above 1e-5 and both were traced to the oracle's own truncated series rather than
to this crate: MNDO germane at 1.209e-3 kcal/mol and MNDO/d H2S at 9.014e-4,
which sit at the branch point of MOPAC's B-integral expansion.

## What the oracles cannot do

Stated because a limitation that is not written down looks like coverage:

- MolDS's input deck has no total-charge keyword, so the CNDO/2 and INDO suites
  contain no ions.
- MolDS and this crate use different definitions of the dipole in a
  zero-differential-overlap model -- 0.23 D apart on CH4 -- so those two suites
  compare everything except the dipole.
- MOPAC7 assumes internal coordinates for three atoms or fewer, rebuilds its own
  Cartesian frame from a Z-matrix, rejects Cartesian input whose first three
  atoms are collinear, and stops on any bond under 0.8 A. Each is worked around
  where it can be and declared where it cannot; ethyne is the one molecule no
  atom ordering rescues.

## Finite differences

Analytic gradients and Hessians are compared with central differences of fully
reconverged energies and gradients respectively, in `tests/derivative_fd_matrix.rs`
and `tests/cis_derivatives.rs`. That is a different kind of evidence from an
oracle -- it checks the derivative against the energy this crate itself produces,
not against another program -- and both are kept.
