# Parameter provenance

Only complete tables with a compatible redistribution basis are bundled.

| Dataset | Runtime role | Upstream | License |
| --- | --- | --- | --- |
| MNDO elements and pairs | primary | OpenMOPAC 23.2.5 | Apache-2.0 |
| MNDO/d elements and pairs | primary | OpenMOPAC 23.2.5 | Apache-2.0 |
| ZINDO/S | primary | OpenMOPAC 23.2.5 | Apache-2.0 |
| CNDO/2 and INDO | primary | MolDS 0.3.1 | GPL-3.0-or-later |
| ZINDO/S subset | cross-check | MolDS 0.3.1 | GPL-3.0-or-later |
| MINDO/3 | runtime mirror and audit | MOPAC7 | public domain |

MOPAC7's program source explicitly places the program in the public domain;
the retained notice links both that source statement and OpenMOPAC's historical
distribution description.

The authoritative machine-readable catalog is
`src/data/legacy_parameter_catalog.csv`; hashes are stored in
`src/data/legacy_parameter_manifest.sha256`. Extraction and validation tools
are in `tools/`.

MolDS provides CNDO/2 for H, Li, C, N, O, and S. Its ground-state INDO path is
enabled only for H, Li, C, N, and O because the sulfur row is not treated as a
complete INDO one-center set.

OpenMOPAC's ZINDO/S table contains transition-metal entries, but 0.3.0 still gates
the irregular 9-AO d-shell path because its additional Hamiltonian machinery
is not implemented. A table row alone is not considered implementation
evidence.

See `parameter-search.md` for the methods whose parameters remain absent.
