# Independent oracle validation

The principal independent oracle is the official OpenMOPAC 23.2.5 Windows
release. `tools/oracle/run_mopac.py` creates a fresh input, executes the binary
selected by `MOPAC_EXE`, and parses AUX results. The binary is not included.

MNDO regression coverage includes:

- water fixed-geometry heat of formation, dipole, and charges;
- water and distorted-water Cartesian gradients;
- an open-shell methyl radical UHF energy;
- distorted-water geometry optimization;
- optimized-water frequencies.

For MOPAC `FORCE` comparisons, use an optimized geometry and `NOREOR` so the
Cartesian frame and mode selection are comparable. Numerical derivatives are
also checked against central finite differences within the Rust test suite.

CNDO/2 and INDO tests use MolDS equations and source constants as an independent
implementation reference. MINDO/3 tables are checked against public-domain
MOPAC7 source sentinels. ZINDO/S is checked against both OpenMOPAC extraction
conventions and the independent MolDS subset.

Oracle agreement is method-specific. Tests never apply empirical scale factors
or substitute a different method to force agreement.
