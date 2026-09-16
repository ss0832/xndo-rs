# OpenMOPAC provenance

What was taken from OpenMOPAC, from which version, and what the licence
requires in return. The retained notices themselves are in the root `NOTICE`;
this file records the evidence behind them.

## Version pinned

| | |
| --- | --- |
| upstream | https://github.com/openmopac/mopac |
| tag | `v23.2.5` |
| commit | `052691223d19935a89f0fe18cd12301bd83e4201` |
| licence | Apache-2.0 |
| release date | 2026-05-03 (from `CITATION.cff`) |
| DOI | 10.5281/zenodo.6511958 |
| JOSS | 10.21105/joss.08025 |

## The licence copy is verbatim

`third_party/mopac/LICENSE` is a byte-for-byte copy of upstream's:

```
sha256  6d1d968fb225eca367cb7f0b8831ab012a35d92b547e945e17ef8e7b05c3e5cc
```

`tests/attribution.rs` pins that hash, so a copy that drifts from upstream fails
the test rather than passing quietly.

Note that upstream's `LICENSE` ends at `END OF TERMS AND CONDITIONS` and omits
the Apache-2.0 APPENDIX ("How to apply the Apache License to your work"). That
is upstream's copy as published; the omission is theirs and is reproduced rather
than corrected. `LICENSES/Apache-2.0.txt` holds the full text including the
appendix, taken from a crate that ships it.

## Upstream ships no NOTICE file

Apache-2.0 section 4(d) obliges a redistributor to propagate a `NOTICE` file
**if the original work includes one**. OpenMOPAC v23.2.5 does not. The complete
top level of the tagged tree is:

```
.gitignore        CODE_OF_CONDUCT.md   Dockerfile
AUTHORS.rst       CONTRIBUTING.rst     LICENSE
CITATION.cff      CMakeLists.txt       README.md
```

There is no `NOTICE` at any level of that tree, so 4(d) imposes nothing here.

**No file named `third_party/mopac/NOTICE` exists in this repository, and none
should be created.** Placing one at that path would assert that upstream
published a NOTICE, and would pass an obligation to propagate a document
upstream never wrote on to every downstream recipient. Recording the absence is
the accurate thing to do. `tests/attribution.rs` asserts that the path stays
empty, so a well-meaning future commit cannot invent one.

## What section 4 does require, and where it is satisfied

| clause | requirement | where |
| --- | --- | --- |
| 4(a) | give recipients a copy of the licence | `third_party/mopac/LICENSE`, `LICENSES/Apache-2.0.txt` |
| 4(b) | carry prominent notices stating that files were changed | `PROVENANCE:` / `MODIFIED:` header on each derived file; the table in `THIRD_PARTY_NOTICES.md` |
| 4(c) | retain the attribution notices of the source form | the verbatim header block in the root `NOTICE` |
| 4(d) | propagate a NOTICE file if one exists | not applicable; see above |

Since xndo-rs is distributed under GPL-3.0-or-later, GPL-3.0 section 5(a) also
requires prominent notices of modification with dates. The same `MODIFIED:`
headers carry those.

## Retained attribution documents

- `AUTHORS.rst` -- the upstream contributor history, verbatim. It names
  individuals; those names are attribution required by the licence and must not
  be stripped by any privacy scan.
- `CITATION.cff` -- the upstream citation metadata, verbatim.

## What was taken

Two kinds of thing, and the distinction matters because the earlier version of
`THIRD_PARTY_NOTICES.md` claimed only the first:

1. **Parameter data** -- element and pair tables, extracted into CSV.
2. **Line-level ports of Fortran into Rust** -- the one-centre integrals, the
   rotation matrices, the d-orbital integrals, the parameter derivations, the
   core-core repulsion, and the SCF and INDO/S drivers. These are derivative
   works of the Fortran, not merely users of its numbers.

`THIRD_PARTY_NOTICES.md` lists both, file by file, with the upstream `.F90` each
one derives from.

## Elements of the model that are upstream behaviour, not ours

Recorded here because they look like arbitrary constants in the Rust and are
not:

- **MNDO/d is evaluated with pre-2019 physical constants.** `readmo.F90:411`
  selects them for the `MNDOD` keyword as well as for `OLDFPC`. See
  `tests/data/ORACLE_NOTES.md` item 13.
- **MOPAC's INDO module carries its own constants** in `src/INDO/reimers_C.F90`,
  which are neither of the two sets `conref_C.F90` offers.
- **The `src/INDO/` lineage.** Per upstream's `src/INDO/README.md` that code is
  by Rebecca Gieseking, derived from Jeff Reimers' CNDO/INDO program
  (DOI:10.4231/D3R49G96G), itself derived from the CNDO/S program of Del Bene,
  Jaffe, Ellis and Kuehnlenz (QCPE 174). xndo-rs's ZINDO/S engine descends from
  that line and the root `NOTICE` records it.
