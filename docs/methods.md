# Methods

| Method | Status | Elements or branch | API |
| --- | --- | --- | --- |
| CNDO/2 | native | H, Li, C, N, O, S | RHF/UHF energy, gradient, optimization, Hessian, frequencies |
| INDO | native | H, Li, C, N, O | RHF/UHF energy, gradient, optimization, Hessian, frequencies |
| MNDO | native | OpenMOPAC table | RHF/UHF energy, derivatives, optimization, frequencies |
| MNDO/d | native | OpenMOPAC table | RHF/UHF energy, derivatives, optimization, frequencies |
| MINDO/3 | native | H, B, C, N, O, F, Si, P, S, Cl | RHF/UHF energy, gradient, optimization, Hessian, frequencies |
| ZINDO/S | native | 1-AO/4-AO s/p branch | RHF/UHF ground optimization/derivatives, RHF-CIS/UCIS UV-vis properties and state derivatives |

Every native method additionally reports orbital energies, occupations and the
HOMO/LUMO pair, and can write a Molden wavefunction file. ZINDO/S rejects
elements that would need its unimplemented d branch for Molden output, rather
than writing a file missing their d coefficients.

No correlated method is implemented. See `docs/scope.md`.

`INDO` means the unversioned MolDS-compatible ground-state model. It is not an
alias for `INDO/S`; the latter parses as ZINDO/S.

Registered but non-executable methods are CNDO/1, INDO/1, INDO/2, ZINDO/1,
ZINDO/2, MINDO/1, MINDO/2, SINDO1, and MSINDO. The registry preserves names
and future API slots but returns an explicit error on execution.

CNDO/3, INDO/3, and ZINDO/3 remain parseable only to diagnose legacy input;
they are marked non-canonical and cannot run.

Accepted aliases and capabilities are discoverable through Rust
`Method::accepted_strings`, Python `available_methods()` and `api_methods()`,
or `xndo_rs_cli methods`.
