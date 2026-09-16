# Missing-parameter search

The audit, repeated for 0.3.0, retained the search for complete, redistributable parameter and
Hamiltonian definitions for CNDO/1, INDO/1, INDO/2, ZINDO/1, ZINDO/2,
MINDO/1, MINDO/2, SINDO1, and MSINDO.

No additional set satisfied all of these requirements:

- complete numerical parameters for the named method;
- an unambiguous Hamiltonian and parameter convention;
- a license compatible with redistribution in this GPL library;
- enough provenance to cite the exact version and source.

Therefore no additional parameter set is included and no nearby-method values
are substituted. Names of rejected candidate code bases are intentionally
omitted. The three adopted sources are documented in
`parameter-provenance.md`: OpenMOPAC 23.2.5, MolDS 0.3.1, and public-domain
MOPAC7.

The registered methods remain gated until both a complete numerical definition
and compatible provenance are available.
