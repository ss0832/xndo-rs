<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Oracle notes — xndo-rs reference data

Companion document for `tests/data/*.tsv` and `tests/oracle_matrix.rs`.

Nothing in `tests/data/*.tsv` is typed by a human. Every value column is what an
oracle program printed; every geometry column is what was sent to it. The
generators live in `tools/oracle/`. See §(d) for how staleness is detected.

**Status.** MNDO and MNDO/d are enforced suites; see §(e) for what they
measure. ZINDO/S, MINDO/3, CNDO/2 and INDO are not built yet. The numbered items
in §(c) are the conventions, traps and defects found so far, and are cited by
the tests and by commit messages; four of them (3, 3b, 14 and 16) are resolved
and kept for the record, and item 13 is open.

---

## Oracle programs

| Program | Version | Licence | How obtained | Redistributed? | Produces |
| --- | --- | --- | --- | --- | --- |
| OpenMOPAC | 23.2.5 | Apache-2.0 | official Windows release binary, not bundled; selected by `MOPAC_EXE` | no | MNDO, MNDO/d, ZINDO/S |
| MOPAC7 | 1.15 | public domain | Debian source, built locally with gfortran | no | MINDO/3 |
| MolDS | 0.3.1 | GPL-3.0-or-later | Debian source | no | CNDO/2, INDO |

MOPAC 23.2.5 is the declared upstream of the MNDO, MNDO/d and ZINDO/S parameter
tables (`src/data/*.csv` provenance headers), and MOPAC7 is the declared upstream
of the MINDO/3 tables. Each method is therefore checked against the program its
own parameters came from, not against a proxy.

---

## (c) Ambiguities, traps and deliberate omissions

Numbered so that tests and commit messages can cite them.

### 1. MOPAC 23.2.5 implements neither MINDO/3 nor CNDO

`molkst_C.F90:331-367` enumerates the 20 methods the binary dispatches: MNDO,
AM1, PM3, RM1, MNDOD, PM6 and its corrected variants, PM7 variants, PM6-ORG,
PM8, INDO. There is no `MINDO3` and no `CNDO`. `MINDO` appears in the upstream
tree only in `AUTHORS.rst` and two README files, never in the keyword parser
`src/input/readmo.F90`.

Consequence: `run_mopac.py --method MINDO3` does **not** produce MINDO/3. MINDO/3
is routed to MOPAC7 instead, which is also its declared parameter upstream.

### 2. MOPAC's `INDO` is the code ZINDO/S was ported from

`src/data/zindo_s_parameters.csv:1-4` records its own provenance as
`src/models/parameters_for_INDO_C.F90` of OpenMOPAC v23.2.5. The engine lives in
the upstream `src/INDO/` directory. Per its `README.md`, that code was written by
Rebecca Gieseking, derived from Jeff Reimers' CNDO/INDO program
[DOI:10.4231/D3R49G96G], itself derived from the CNDO/S program [QCPE 174] by
J. Del Bene, H. H. Jaffe, R. L. Ellis and G. Kuehnlenz.

MOPAC `INDO` is therefore a same-lineage oracle for ZINDO/S, ground state and
excited states both. This lineage is recorded in the attribution documents.

### 3. xndo-rs omitted the Zerner sigma/pi resonance weighting — RESOLVED (see 3b)

MOPAC's INDO/S multiplies the sigma and pi components of the overlap that feeds
the resonance integral by Zerner's weighting factors:

* `src/INDO/reimers_C.F90:120` — `data fintfa / 1.D0,1.267D0,0.585D0,3*1.D0/`
* `:133-135` — `equivalence (fintfa(1),fssig),(fintfa(2),fpsig),(fintfa(3),fppi), ...`
* `src/INDO/ovlp.F90:134` — `! fspdf(i, j) is weighting factor.`
* `:139-147` — defaults are all 1, then `fspdf(1,1)=fssig`, `fspdf(2,1)=fpsig`,
  `fspdf(2,2)=fppi`. Applied to same-`l` pairs only.

Net effect: p-p sigma overlap x 1.267, p-p pi overlap x 0.585; s-s and s-p x 1.0.

`src/zindo.rs:490-511` stores the raw Slater block and `:631` builds the
resonance element as `0.5 * (beta_a + beta_b) * overlap[(mu, nu)]`, with no
sigma/pi factor anywhere. Grepping the whole tree for `1.267|0.585` finds only an
unrelated MNDO pair parameter at `src/data/mndo_pair_parameters.csv:7`.

**Measured consequence: the sigma/pi weighting is real, but it is NOT the main
problem.** A first comparison on formaldehyde (C at origin, O at z=1.203, H at
(+-0.943, 0, -0.545)) shows a disagreement far too large for a weighting factor:

| quantity | MOPAC `INDO` | xndo-rs ZINDO/S | ratio |
| --- | --- | --- | --- |
| electronic energy | -793.3576 eV | -796.0340 eV | 1.003 |
| S1 (HOMO->LUMO, 6->7) | **3.195 eV** | **13.837 eV** | **4.33** |
| S2 | 8.301 eV (6->8) | 17.140 eV (5->7) | — |
| S3 | 8.851 eV (4->7) | 17.707 eV (6->8) | — |

Both runs used the same geometry, and MOPAC was re-run with the full CI window
`C.I.=(10,6)` so the active spaces match (25 spin-adapted configurations); the
state list is unchanged from the `C.I.=(6,3)` run for the roots they share.

Two conclusions:

* **The SCF is approximately right** — the electronic energies agree to 0.34%,
  and the dominant configuration of the lowest root (6->7) is the same in both.
* **The CIS construction is wrong, not just weighted differently.** Overlap
  weights of 1.267 and 0.585 cannot move an excitation energy by a factor of
  4.33. MOPAC's HOMO/LUMO are -10.779 / +0.148 eV, a gap of 10.93 eV, and its
  S1 lies 7.7 eV *below* that gap — the expected behaviour for a singlet, where
  `A_ia,ia = (eps_a - eps_i) + 2(ia|ia) - (ii|aa)` is dominated by the large
  negative Coulomb term `-(ii|aa)`. xndo-rs places S1 *above* a comparable gap,
  which is the signature of that Coulomb term having the wrong sign, the wrong
  magnitude, or being swapped with the exchange term — not of a scaled overlap.

So item 3 (sigma/pi) and this CIS defect are **two separate problems**, and the
sigma/pi factor must not be presented as the explanation for the CIS numbers.

**The CIS defect is now found and fixed.** `build_jk` was storing `gamma_AB` in
`k[mu][nu]` for AO pairs on different atoms. That is what the *Fock* builder
needs — the two-centre element is `F_mu_nu = beta_AB S_mu_nu - 0.5 P_mu_nu
gamma_AB` — but `mo_eri` reads `k` as the exchange **integral** `(mu nu|mu nu)`,
which ZDO sets to zero across atoms. The `-0.5 P gamma` term is the Coulomb
integral `(mu mu|nu nu)` appearing in the exchange contraction at
`lambda = mu, sigma = nu`, not a two-centre exchange integral. Reading it as one
injected a spurious ~10 eV integral into every CIS matrix element.

`k` is now unambiguously the exchange integral: zero off-atom, `G1`/`3*F2`
on-atom. The Fock builders read `j` in the two-centre exchange position, where
`j` holds exactly the `gamma_AB` they used to read from `k`, so **the SCF is
bit-identical** (formaldehyde total energy `-516.378679113497` eV before and
after). The `kx`/`kxy` derivative matrices are now identically zero, which is
correct: one-centre integrals do not depend on geometry.

Effect on formaldehyde S1: **13.837 eV -> 4.651 eV** against MOPAC's 3.195 eV.
That removes 9.2 eV of the 10.6 eV discrepancy. The residual 1.46 eV is the
right order for the missing sigma/pi weighting of item 3, which is therefore
still expected to matter — as a second-order effect, not as the explanation.

### 3b. The sigma/pi weighting, applied — RESOLVED

For an s/p basis the weighting collapses to two factors on the local-frame
quantities that `slater_locals_numeric` already returns. MOPAC weights same-`l`
pairs only and averages the two centres' factors, so:

| local quantity | pair | weight |
| --- | --- | --- |
| `s111` (s-s sigma) | same `l` | `(fssig + fssig)/2 = 1` |
| `s211`, `s121` (s-p sigma) | cross `l` | 1 |
| `s221` (p-p sigma) | same `l` | `fpsig = 1.267` |
| `s222` (p-p pi) | same `l` | `fppi = 0.585` |

Applied in `sp_overlap_block_g` in the diatomic local frame, before
`build_di_g` rotates — the same place MOPAC applies them, and inside the
`Dual`/`Dual2` path so the analytic derivatives follow automatically.

Formaldehyde, matched active space, after both fixes:

| dominant config | xndo-rs (eV) | MOPAC (eV) | diff | xndo-rs f | MOPAC f |
| --- | ---: | ---: | ---: | ---: | ---: |
| 6->7 | 3.1866 | 3.195 | -0.008 | 0.0000 | 0.0000 |
| 6->8 | 8.2139 | 8.301 | -0.087 | 0.0271 | 0.0416 |
| 4->7 | 8.8194 | 8.851 | -0.032 | 0.0164 | 0.0204 |
| 5->7 | 8.8798 | 9.885 | **-1.005** | 0.3391 | 0.5483 |
| 5->8 | 10.4682 | 10.505 | -0.037 | 0.0330 | 0.0397 |
| 6->9 | 11.6954 | 11.776 | -0.081 | 0.1425 | 0.1420 |
| 3->7 | 11.8784 | 11.883 | -0.005 | 0.0000 | 0.0000 |

Every root now matches by dominant configuration and in the same order. Six of
seven agree to within 0.09 eV, against an initial error of 10.6 eV.

### 12. Both residuals were the harness reading the wrong table -- RESOLVED

Two residuals were recorded here as open model defects:

**(a) The 5->7 root is 1.0 eV low.** Reported as a missing off-diagonal CIS
coupling, on the grounds that xndo-rs put the 4->7 and 5->7 roots 0.06 eV apart
where MOPAC separated them by 1.03 eV.

**(b) Oscillator strengths run 20-35% low on the bright states** -- 0.0271 vs
0.0416, 0.3391 vs 0.5483 and so on -- while the symmetry-forbidden ones were
exactly zero in both.

**Neither existed.** Both numbers on the MOPAC side came from the *spin-adapted
configuration* table, which is the CI matrix diagonal before diagonalisation,
not from the CI eigenstate table below it (item 7). The "1.03 eV separation" is
the separation of two configurations; the roots they mix into are 0.014 eV
apart, which is what xndo-rs reports. The "0.5483" is a configuration's
oscillator strength; the root's is 0.331, which is what xndo-rs reports.

Re-measured on formaldehyde against the CI table, same geometry, full 6x4
window, all 24 roots:

| quantity | worst over 24 roots |
| --- | ---: |
| excitation energy | 1.191e-3 eV |
| oscillator strength | 8.060e-5 |
| state dipole magnitude | 5.237e-3 D |

and over the whole suite of 24 molecules and 192 roots, 4.625e-4 eV /
1.666e-5 / 2.233e-3 D. This is deviation bucket 1: a convention difference,
fixed on the comparison side, tolerance untouched.

Worth stating plainly, because it cuts the other way too: the two real ZINDO/S
defects found in this work -- the two-centre exchange integral and the missing
sigma/pi resonance weighting (items 3 and 3b) -- were found by the *same* kind of
comparison. The lesson is not to distrust large disagreements, it is to confirm
what the oracle is printing before naming a cause.

### 4. `mindo3_parameters.csv` declared a modified value — fixed, not tolerated

Through v0.2.4, `src/data/mindo3_parameters.csv:4` said `hp2_ev` followed a
"rounded xndo-rs v0.2.1 internal convention (nominally (gpp-gp2)/2)" rather than
MOPAC7's value. The prediction here was a systematic MINDO/3 offset against the
MOPAC7 oracle, to be fixed or registered in `KNOWN` and never absorbed into a
tolerance. It was fixed.

The first MINDO/3 suite run found it immediately, and localised it without
ambiguity: CCl4 was **8.6e-2 eV** out on the total energy while its core-core
repulsion agreed to **7.7e-6 eV**. The geometry term is exact and the density is
wrong. Across the set the residual was **2.2e-2 eV per chlorine atom** — Cl2
4.1e-2, HCl 2.3e-2, CH3Cl 2.2e-2, CCl4 8.6e-2 — against about 8e-4 eV for every
molecule without a chlorine in it.

`hp2` is the one-centre p-p' exchange integral, and **MOPAC7 does not store
one**: `fock1.f` writes `GPP(NI) - GP2(NI)` and `0.5*(GPP(NI) - GP2(NI))`
straight into the Fock matrix. The relation is structural, so "nominally" was
the wrong word for it. Four of the nine values had been rounded to two decimals
— N 0.695→0.70, Si 0.385→0.38, S 0.535→0.54, Cl 0.665→0.67 — each 0.005 eV, and
each a change to the model rather than to how it is printed. Chlorine shows it
worst because `p^5` gives the error the most p-p' pairs to act on.

All four now carry the exact half, and
`mindo3::tests::hp2_is_exactly_half_the_gpp_gp2_gap` asserts the relation for
every element by **exact equality** — these are two- and three-decimal numbers,
both sides are representable, and a tolerance is what let the rounding in.

### 5. Read scalars from AUX, not from the printed output

`AUX(PRECISION=9)` carries far more digits than the `.out` file's 3-4 decimals.
For the NDDO methods the AUX block provides `HEAT_OF_FORMATION`,
`ATOM_CHARGES`, `DIP_VEC[3]`, `EIGENVALUES`, `EIGENVECTORS`,
`MOLECULAR_ORBITAL_OCCUPANCIES`, `IONIZATION_POTENTIAL`, `AO_ZETA`, `ATOM_PQN`
(observed on `MNDO PRECISE 1SCF`, water).

### 6. The INDO CI results are printed, not written to AUX

Adding `CIS C.I.=(N,M)` changes what reaches AUX. Observed on formaldehyde with
`INDO PRECISE 1SCF CIS C.I.=(6,3) WRTCI=8 AUX(PRECISION=9)`: the AUX block
contains `HEAT_OF_FORMATION`, `DIPOLE`, `DIP_VEC[3]`, `IONIZATION_POTENTIAL`,
`SET_OF_MOS` — but **no `EIGENVALUES` and no `ATOM_CHARGES`**. The CI state table
goes to the output unit only.

So the two ZINDO/S suites are generated differently:

* ground state — plain `INDO PRECISE 1SCF AUX(PRECISION=9)`, read from AUX;
* excited states — `INDO ... CIS C.I.=(N,M) WRTCI=k`, read from the `.out` state
  table with a dedicated parser.

MOPAC also prints `WARNING: Specified print options may not behave as expected
with INDO` when `TDIP` is requested; `TDIP` produced no extra columns in the
observed output and is not used.

### 7. MOPAC's INDO CI prints two tables, and only the second is comparable

**This is the single most expensive misreading in this file's history.** It is
what ORACLE_NOTES item 12 recorded as two model defects that did not exist.

The first table is headed

```
     sym   eV   cm**-1 -dets- dipole oscilator X FRAG ....Excitations named from first reference determinate
                       tot  #  Debye  strength
```

and lists the **spin-adapted configurations**: the diagonal of the CI matrix
*before it is diagonalised*, each labelled by the single excitation it came
from. The second is headed

```
  CI trans.  energy frequency wavelength oscillator------ polarization------  dipole  --- components---
 st.  symm.    eV      cm - 1       nm      strength       x       y      z    moment   x     y     z
```

and lists the **CI eigenstates**. For formaldehyde in a 6x4 window:

| | configuration table | CI table |
| --- | ---: | ---: |
| `(5)->(7)` / root 4 energy | 9.887 eV | 8.8828 eV |
| its oscillator strength | 0.547281 | 0.017745 |

A CIS calculation produces eigenstates. Comparing it against the configuration
table produces exactly the two symptoms item 12 recorded -- one root a full eV
"low", and bright-state oscillator strengths "20-35% low" -- from a correct
implementation.

**What the CI table gives per root**, and the precision of each:

| quantity | printed as | usable to |
| --- | --- | --- |
| excitation energy | `8.8684904` | ~1e-6 eV |
| same, in cm^-1 | `71529.` | 1.24e-4 eV |
| oscillator strength | `0.331024` | 1e-6 |
| state dipole magnitude | `1.429693` | 1e-6 D |
| state dipole components | `0.000 -0.000 -1.430` | 1e-3 D |
| polarization | `-0.000 0.000 -1.000` | unit vector |

So **read the eV column here**, not cm^-1 -- the opposite of the configuration
table, where the 3-decimal eV column is coarser than the integer cm^-1 one. The
suite compares energy, oscillator strength and dipole magnitude: 3 scalars per
root. The components and the polarization vector add a factor of 1000 less
precision and no independent information.

Row 1 of the *configuration* table is the ground state at 0.000 eV and is not a
root. The CI table has no such row; it starts at state 2.

### 7b. MOPAC's INDO CI screens the configuration list; a CIS does not

MOPAC does not always diagonalise the window it is asked for. Asked for a 4x4
window it reports:

| molecule | point group | configurations kept | a full 4x4 CIS has |
| --- | --- | ---: | ---: |
| formaldehyde | C2v | 17 | 17 |
| benzene | D6h | 17 | 17 |
| CCl4 | Td | **1** | 17 |

and asked for 5x5, benzene keeps 21 of 26. A CIS in this crate has no screening,
so a screened CI is a different calculation wearing the same name -- comparing
against one would not be a loose comparison but a meaningless one.

The generator therefore reads `CI excitations=` back out of the output and takes
the **largest window whose count is exactly `occupied * virtual + 1`**, recording
that window in the reference file so both sides provably diagonalise the same
matrix. CCl4 lands on 3x4, benzene on 5x4, and most molecules on the full 5x5.

### 7c. Degenerate roots are manifolds, not roots

Inside a degenerate set the eigenvectors are fixed only up to a unitary rotation,
so a per-root property that is not invariant under that rotation is not an
observable of the root, and two correct programs will disagree on it.

The excited-state permanent dipole is exactly such a property. CCl4's first
excited state is triply degenerate; its three components carry dipoles of equal
magnitude pointing in whatever directions the diagonaliser happened to pick, and
comparing them root by root gives a **1.78 D** "disagreement" that means nothing.

The suite therefore compares the excitation energy per root (degenerate energies
are equal), the oscillator strength per **manifold sum** (invariant, being the
trace of a quadratic form over the block), and the permanent dipole only for
roots that are alone at their energy. Roots are grouped at 1e-5 eV.

### 8. Active-space correspondence: `C.I.=(N,M)`

MOPAC echoes `C.I.=(6,3)` as "3 DOUBLY FILLED LEVELS USED IN A C.I. INVOLVING 6
M.O.'S", and reports `CI excitations= 10` for a 6-occupied-orbital reference:
3 occupied x 3 virtual singles, plus the ground configuration.

So `N` is the total number of MOs in the CI window and `M` is how many of them
are occupied. In terms of `ZindoOptions` (`src/zindo.rs:3638-3645`):

```
C.I.=(active_occupied + active_virtual, active_occupied)
```

Comparing CIS energies across a mismatched window compares nothing. One window is
pinned for the whole suite and recorded in section (a).

### 9. MOPAC needs a short working directory

Running the binary with its input under a long path aborts at cleanup with
`forrtl: severe (9): permission to access file denied, unit 14` and a truncated
path in the message. The `.out` and `.aux` files are written before the abort, so
the failure is easy to miss. `run_mopac.py` already uses `tempfile.mkdtemp`,
which lands in the system temp directory and is short enough; generators must not
pass a long explicit working directory.

### 10. ZINDO/S has no heat of formation

`INDO 1SCF` on water reports `HEAT_OF_FORMATION = -8362.4 kcal/mol`. That is the
total electronic energy in disguise, not a heat of formation: INDO/S is a
spectroscopic parameterisation with no atomic heat terms. Compare total energies
for ZINDO/S, never heats of formation.

### 11. Triplet oscillator strengths are hard zero in xndo-rs

`src/zindo.rs:3718-3722` sets them to zero because there is no spin-orbit
intensity borrowing. Never compared.

### 16. MOPAC's MNDO refuses the transition metals whose parameters xndo-rs ships

`mndo_parameters.csv` carries usable MNDO parameters (non-zero `uss` and
`betas`) for Sc, Ti, V, Cr, Fe, Co, Ni, Cu, Zr, Mo, Pd, Ag and Pt, with the
provenance header pointing at `parameters_for_mndo_C.F90` of OpenMOPAC v23.2.5.
MOPAC 23.2.5 itself will not run them under `MNDO`: ScF3, VF5, CrF6, FeCl2,
MoF6, PdCl2, AgCl and PtCl2 all end with

```
* Parameters for some elements are missing
```

Zn does run, so the d path itself is reachable; it is these specific elements
that MOPAC gates off.

Two consequences:

* the MNDO suite has no transition-metal rows, and cannot have any while the
  oracle refuses them;
* the eight elements that carry an **independent `poc`** under MNDO — Sc, V,
  Cr, Fe, Mo, Pd, Ag, Pt — are exactly the ones with no MOPAC MNDO reference.
  The `poc` override (item 14) is therefore validated where it matters and
  measurable, under MNDO/d for Na and Mg, and is unvalidatable for MNDO
  transition metals. That is recorded here rather than presented as verified.

Whether xndo-rs should also refuse those elements under MNDO is a scope
question for `docs/scope.md`, not an oracle question.

### 13. MNDO/d is evaluated with the pre-2019 physical constants -- RESOLVED

Under `MNDOD`, xndo-rs disagreed with MOPAC on **every** molecule by a small
amount that grew with size: H2 2.25e-3, water 9.4e-3, methane 1.3e-2, CO2
2.8e-2, SF6 4.8e-2, benzene 7.4e-2 kcal/mol.

The decisive observation was H2. Hydrogen has **identical parameters in the MNDO
and MNDO/d tables** (`ussd/betasd/zsd/alpd/gssd(1)` equal `ussm/betasm/zsm/alpm/
gssm(1)` to the last digit), and xndo-rs's MNDO and MNDO/d agreed with each
other exactly on H2 -- yet MOPAC's `MNDOD` differed from its own `MNDO` by
2.25e-3 kcal/mol. So MOPAC applied something under `MNDOD` that does not depend
on d parameters at all.

The first guess recorded here was routing: that MOPAC's `MNDOD` sends every pair
through `mndod.F90` while xndo-rs keeps sp pairs on the sp kernel. **That guess
was wrong, and wrong in an instructive way** -- `rotate.F90:105` calls `rotatd`
for *every* method, so MOPAC has only one integral path and there is no routing
difference to find.

**The actual cause is `readmo.F90:411`:**

```fortran
if (index(keywrd,' OLDFPC') + index(keywrd, ' MNDOD') > 0) then
  fpc(:) = fpcref(2,:)   ! pre-2019 constants
else
  fpc(:) = fpcref(1,:)   ! 2018 CODATA
end if
```

`MNDOD` selects the same historical constant set as `OLDFPC`. From
`conref_C.F90`:

| constant | CODATA `fpcref(1,*)` | historical `fpcref(2,*)` | relative |
| --- | --- | --- | ---: |
| Hartree / eV | 27.211386245988 | 27.21 | -5.09e-5 |
| a0 / A | 0.529177210903 | 0.529167 | -1.93e-5 |
| eV / (kcal/mol) | 23.060547830619029 | 23.061 | +1.96e-5 |

That is not a rounding artefact but part of the model: Thiel and Voityuk fitted
MNDO/d against those constants. The offsets enter every two-centre integral
(through `rho = 0.5 ev / G` and through `r/a0`) and accumulate over pairs, which
is exactly the observed growth with size.

**Measurement that confirmed it.** MOPAC on H2 at 0.74 A:

| | MNDO | MNDOD | difference |
| --- | ---: | ---: | ---: |
| heat of formation / kcal/mol | 2.82589461639645 | 2.82813547927174 | 2.2409e-3 |
| ionization potential / eV | 15.2044976479027 | 15.2043140046896 | 1.8364e-4 |
| HOMO + LUMO / eV | -10.96455 | -10.96455 | 0 |

The sum of the eigenvalues is unchanged, which pins `F_AA = U_ss + G_ss/2` and
rules out anything in the one-centre terms or the core-electron attraction; the
whole shift sits in `F_AB = beta S - (ss|ss)/2`. Working back from it gives a
relative change in `rho` of about 4.5e-5, against 5.09e-5 predicted from the
Hartree alone. A single H atom is identical under both, which rules out `eisol`
and `eheat`.

**The fix** is [`ModelConstants`] in `src/constants.rs`. The crate's internal
length unit is the CODATA Bohr everywhere, because geometry is converted before
a method is chosen, so a parameterisation fitted against a different Bohr radius
is folded into the parameters instead of the coordinates: every Slater exponent
is scaled by `k = a0_CODATA / a0_model`, every energy prefactor becomes
`ev_model / k`, and `pocord` -- a length in the model's own Bohr -- is divided by
`k`. Both substitutions are exact, and `k == 1` for CODATA, so MNDO is
bit-identical (its worst deviations are unchanged to the digit).

The truncated solvers are the one subtlety. `additive_rho1`, `additive_rho2`
(5 secant steps) and `poij` (golden section, bracket [0.1, 5.0], absolute 1e-8)
are **not** scale-invariant, so each converts to the model's Bohr, runs MOPAC's
iteration verbatim, and converts the root back.

Result, over the same 37 molecules:

| quantity | before | after |
| --- | ---: | ---: |
| heat of formation | 7.421e-2 kcal/mol (benzene) | **8.863e-4** (h2s) |
| atomic charge | 3.674e-4 e | **3.283e-5** |
| frontier orbital | 9.102e-4 eV | **1.343e-4** |
| dipole component | 9.489e-4 D | **3.810e-4** |

MNDO/d now agrees with MOPAC as closely as MNDO does (worst 1.227e-3 kcal/mol),
and the suite runs on MNDO's floors rather than one two orders of magnitude
looser.

`src/integrals_d.rs:9-15` still records that the spd path bypasses MOPAC's
`ind2`/`isym`/`rep(491)` symmetry compression. That remains true and is a
plausible source of the residual 1e-4 kcal/mol, but it is no longer the dominant
term and is well inside the floor.

### 14. MNDO/d: Na, Mg and Al were catastrophically wrong — RESOLVED

On top of item 13, three elements fail by orders of magnitude:

| molecule | xndo-rs | MOPAC MNDOD | (MOPAC MNDO) |
| --- | ---: | ---: | ---: |
| nah | 55.62 | 34.11 | 37.25 |
| mgh2 | 124.92 | 65.23 | 42.77 |
| alh3 | 167.54 | 35.74 | 40.18 |

Both programs agree that `MNDOD` differs from `MNDO` for these elements, and the
parameter tables genuinely differ (Na `zs` 0.821 -> 0.988, Mg 0.939 -> 1.449,
Al 1.444 -> 1.794), so this is not a table-selection mistake — it is how the
MNDO/d parameters are consumed.

What Na, Mg and Al have in common in `mndod_parameters.csv` is a **large `poc`**
(1.531, 1.351, 1.584) together with non-zero internal exponents `zsn`/`zpn`. The
elements that come out close have smaller `poc` (Si 1.267 -> 1.6e-2 kcal/mol,
P 1.185, S 1.116, Cl 1.030). Na and Mg additionally have `zd = 0`: they are
MNDO/d elements that keep an **sp basis** but use the MNDO/d additive one-centre
and core treatment.

**What it actually was - two causes, both in the core-core repulsion.**

The decisive measurement was AlH3. The MO energies agree with MOPAC to 2e-4 eV
(printing precision) while the heat of formation is 131.8 kcal/mol out, so
whatever differs cannot be in the Fock matrix. `ENPART` gives MOPAC's
NUCLEAR-NUCLEAR REPULSION as 83.8227 eV against xndo-rs's 89.5405 eV - a
5.7178 eV gap, which is 131.86 kcal/mol: the whole discrepancy.

1. **The pair-fitted core-core form.** `ccrep.F90` puts the familiar
   `1 + 2 fff exp(-alpb R)` behind `.not. method_mndod` and gives MNDO/d its own
   branch at :134-149:

   ```fortran
   if (ni == nj) then
     scale = 1 + 2 exp(-alpb R)
   else
     select case (nj)
     case(11, 12, 13)   ! Na, Mg, Al
       scale = 1 + exp(-alpb R) + exp(-alp(ni) R)
     case default
       scale = 1 + exp(-alpb R) + exp(-alp(nj) R)
     end select
   end if
   ```

   `alpb` is the alpha of *one* element only; the other exponential uses the
   partner's per-element `alp`. xndo-rs used the PM6 form for both methods. The
   MNDO/d pair table has exactly eight entries, and Na-H, Mg-H and Al-H are
   three of them - which is why those three molecules and no others blew up.

2. **The `poc` core-monopole override.** `inid` sets `po(9) = pocord` where an
   element defines one. For almost every element carrying a `poc` this is a
   no-op, because `poc` *is* the derived `po(1) = 13.6057 / g_ss` to five
   figures. It is an independent number for exactly two MNDO/d elements:
   Na (`poc / po(1) = 0.517`) and Mg (`0.732`).

Fixing (1) alone took alh3 from 131.80 to 0.0139 kcal/mol and left nah at 6.44
and mgh2 at 9.24. Adding (2) took those to 0.0036 and 0.0106.

The previous code discarded `poc` with a comment predicting it would worsen
Sc/Fe/Ni. Those elements have no MOPAC MNDO reference at all (item 16), so that
prediction was made by eye rather than by measurement.

### 17. The AUX eigenvalue vector is incomplete; read the printed output

MOPAC's AUX block drops eigenvalues, and drops *occupied* ones, so it cannot be
used for an orbital comparison. On BF3 it declares `EIGENVALUES[014]` and lists
14 values where the printed output carries the full 16 -- the missing two are
the lowest orbital (-47.35454) and one member of a degenerate pair. Benzene
loses 10 of 30 and CCl4 6 of 20.

`MOLECULAR_ORBITAL_OCCUPANCIES` is truncated in step, so counting occupied
orbitals from it against the complete printed list puts HOMO and LUMO two
orbitals too low. `run_mopac.py` therefore reads eigenvalues from the printed
`EIGENVALUES` block, and uses the AUX occupancies only when their length matches
the eigenvalue list; otherwise the electron count decides.

This was worth catching: before the fix the MNDO suite reported a 17.5 eV LUMO
error on BF3 and the ZINDO/S suite a 25 eV eigenvalue error, both entirely
artefacts of the truncation. xndo-rs matches the printed values.

### 18. ZINDO/S did not converge on benzene -- RESOLVED

400 iterations left a 1.3e-1 density residual. The ZINDO/S SCF had **only linear
damping**: unlike `scf.rs`, which carried A-DIIS/CDIIS, `zindo.rs` had no
convergence accelerator at all, and benzene's degenerate frontier pair is exactly
what an unaccelerated Roothaan iteration oscillates on.

Fixed together with item 19 by `src/scf_accel.rs`, which lifts the accelerator
out of `scf.rs` and gives it to every engine. Benzene now converges and agrees
with MOPAC across its whole occupied spectrum. The `KNOWN` entry is gone.

### 19. ZINDO/S converged to a symmetry-broken solution on BH -- RESOLVED

BH's two pi orbitals are pure boron px and py, non-bonding, and degenerate by
symmetry. MOPAC returns them so:

```
    ROOT NO.    1           2           3           4           5
           -18.30538   -10.00468    -0.21296    -0.21296     4.09982
  S  B   1   0.74550    -0.61867     0.00000     0.00000     0.24793
 PX  B   1   0.00000     0.00000     1.00000     0.00000     0.00000
 PY  B   1   0.00000     0.00000     0.00000     1.00000     0.00000
 PZ  B   1   0.30362     0.64639     0.00000     0.00000     0.70000
  S  H   2   0.59333     0.44658     0.00000     0.00000    -0.66972
```

xndo-rs split the pair by 4.7 eV -- but **only when the bond lay along +z**:

| orientation | eigenvalues |
| --- | --- |
| along +x | -18.3054 -10.0047 -0.2130 -0.2130 4.0998 (matches MOPAC) |
| along +y | identical to +x |
| along (1,1,1) | identical to +x |
| **along +z** | **-17.7235 -6.4730 -4.0426 +0.6754 +4.9107** |

The rotation is exact, so an orientation-dependent answer could only mean the SCF
reached a different fixed point. It did, and the +z one was worse:

| case | iterations | electronic energy |
| --- | ---: | ---: |
| along +z | 20 | **-94.3866 eV** |
| along +x | 21 | **-101.7216 eV** |
| **+z, tilted by 1e-6** | 59 | **-101.7216 eV** |

**Root cause, corrected.** The first diagnosis written here -- "the ZINDO/S SCF
has no convergence accelerator" -- named a real defect but not this one. Adding
the accelerator alone did not fix BH: it converged in 15 iterations instead of
20, to the same wrong answer. The actual cause was in the initial guess.

`zindo.rs` built its starting density as
`initial_coefficients.leading_columns_gram(n_occ, 2.0)`, where
`initial_coefficients` came from diagonalising the SAD *density* matrix.
`symmetric_eigen` returns eigenvalues in **ascending** order, so the leading
columns of that eigenvector matrix are the **least**-populated atomic orbitals,
not the most. The guess was therefore an inverted-population density. It is
orientation-dependent because the AO basis is fixed to the global axes while the
molecule is not: along +z the inverted guess put its weight on boron's pz, which
is the orbital that actually bonds, and the iteration never left that basin.

The fix is one line -- use the SAD density directly, as every other engine does:

```rust
let mut p = guess.clone();   // was: initial_coefficients.leading_columns_gram(n_occ, 2.0)
```

All four orientations now give -101.721614273415 eV, and the pi pair comes back
degenerate at -0.21296 eV as MOPAC has it.

Ruled out along the way, and worth not re-checking: boron (BF3 agreed to
3e-4 eV), the one-centre diagonal (`uform_sp` branches on `orb == 0` only), the
sigma/pi weighting (needs p on both atoms; hydrogen has none), boron's parameters
(in line with Be, C, N, O, F), and the overlap rotation (`rotation_to_x_g(0,0,1)`
gives the correct `R` with `rot_ov[a][0] = dir[a] = (0,0,1)`, so px and py have
exactly zero coupling to hydrogen's s in both orientations).

### 20. SAD and A-DIIS/CDIIS are now the default everywhere

Following items 18 and 19, all four engines were put on one footing:

* the initial density is a **superposition of atomic densities** built by
  `scf_accel::sad_density_sp`, which fills s before p rather than spreading the
  valence electrons uniformly over every AO. MINDO/3 and CNDO/2/INDO both used
  the uniform form; CNDO/2/INDO additionally diagonalised a guess Fock, which is
  now only a test fixture;
* the convergence accelerator is **A-DIIS until the commutator norm falls below
  `adiis_switch`, then CDIIS**, for RHF and UHF alike;
* the unrestricted path stacks the two spin channels and extrapolates them with
  **one** set of coefficients. A comment in `scf.rs` had claimed A-DIIS could not
  be reused for UHF "without a spin-resolved energy functional", and so left UHF
  unaccelerated unless the caller asked for plain CDIIS by name. That is not so:
  the Frobenius product of the stacked matrices is `Tr[P_a F_a] + Tr[P_b F_b]`,
  which is exactly the pairing the UHF energy uses.

Every oracle point in this directory is unchanged by that work -- same worst
deviations to the digit -- so it is a path change, not a model change. The
helium exception in `scf.rs` (MOPAC's MNDO He has formal p AOs but an s-only core
attraction, and an extrapolated path can land on a different charge-transfer
fixed point) is preserved, and now applies to the UHF path as well.

### 21. The last heat-of-formation residuals are MOPAC's B-integral branch

After items 13, 22 and 23 the suite was left with exactly two disagreements
above 1e-5: germane under MNDO (1.209e-3 kcal/mol) and H2S under MNDO/d
(9.014e-4). Every other quantity on every other molecule had come down to about
1e-6. Both turned out to be **MOPAC's** truncation error, not this crate's.

**How it was found.** Three facts pointed away from a model difference:

* the isolated atoms agree exactly -- Ge, Si, S and C give identical heats of
  formation under both codes -- so it is not `eisol` or `eheat`;
* the charges, orbital energies and dipoles agree to about 1e-6, so the
  converged density is the same, and the residual has to sit in a
  density-independent term;
* silane, methane and **H2S under plain MNDO** agree to about 2e-9 kcal/mol at
  matched geometries. Whatever it is, it is not general.

Scanning the bond length made it unmistakable. Germane, Ge-H from 1.5250 to
1.5450 A at matched twelve-decimal geometries:

| R(Ge-H) / A | MOPAC | xndo-rs | delta |
| --- | ---: | ---: | ---: |
| 1.5250 | 14.884575406 | 14.883366230 | -1.209e-3 |
| 1.5300 | 15.222309683 | 15.221073650 | -1.236e-3 |
| 1.5350 | 15.593688188 | 15.592424850 | -1.263e-3 |
| **1.5375** | **15.790543794** | **15.790543800** | **+6.5e-9** |
| 1.5450 | 16.433852863 | 16.433852870 | +7.1e-9 |

A step, not a slope. First differences on a 2.5e-4 A grid identify which side
carries it:

```
xndo-rs   0.019606010  0.019688520  0.019770980  0.019853360  0.019935680  0.020017920  0.020100090
MOPAC     0.019607387  0.019689907  0.019772360  0.019854743  0.019937058  0.018744911  0.020100097
                                                                           ^^^^^^^^^^^
```

**xndo-rs is smooth to 1e-8; MOPAC has a single 1.27e-3 kcal/mol drop.** MOPAC's
own charges and eigenvalues are smooth straight through it, so the step is not
a different SCF solution -- it is in a term evaluated once, from the geometry.

**The cause is `bintgs` (`src/integrals/set.F90:97`)**, which forms the `B`
auxiliary integrals for the overlap. It chooses between an upward recursion and
a truncated power series on the magnitude of `x = (zeta_a - zeta_b) R / 2`:

| range of `|x|` | recursion when | otherwise series truncated at |
| --- | --- | ---: |
| `> 3` | always | -- |
| `2 .. 3` | `k <= 10` | 15 terms |
| `1 .. 2` | `k <= 7` | 12 terms |
| `0.5 .. 1` | `k <= 5` | **7 terms** |
| `<= 0.5` | -- | 6 terms |

The seven-term series is at its worst just below `|x| = 1`, where the first
neglected term is of order `x^8 / 8! = 2.5e-5`.

For the Ge(2p)-H(1s) pair, `zeta_p(Ge) = 2.0205640` and `zeta_s(H) = 1.3319670`,
so `|x| = 1` at

```
R = 2 / (2.0205640 - 1.3319670) = 2.90435 bohr = 1.536896 A
```

The step was measured between 1.53675 and 1.53700 A. Running it backwards, the
exponent that would put the branch exactly at the measured step is 2.020596 --
against the tabulated 2.0205640, agreeing to 3e-5, which is the precision the
step was located to.

H2S under MNDO/d is the same thing. There `zeta_p(S) = 2.0997056`, giving
`|x| = 0.969` at the reference geometry -- just inside the truncated-series
side -- and a predicted crossing at 1.3785 A. Scanning S-H:

| R(S-H) / A | `x(zeta_p)` | delta / kcal/mol |
| ---: | ---: | ---: |
| 1.250 | 0.9068 | -5.826e-4 |
| 1.330 | 0.9648 | -8.756e-4 |
| 1.370 | 0.9938 | -1.058e-3 |
| **1.390** | **1.0083** | **+2.0e-7** |
| 1.450 | 1.0518 | +1.8e-7 |

The disagreement grows monotonically as `x` approaches 1 from below and vanishes
the moment it crosses. Under plain MNDO the same molecule agrees to 2.4e-9,
because MNDO's `zeta_p(S) = 2.0091460` puts `x` at 0.855, far enough from the
branch that the truncation does not show.

**Consequences for the suite.** The heat-of-formation floors stay at 2.0e-3
kcal/mol while the other three quantities come down to 5e-5 or tighter. That
asymmetry is deliberate and is documented at the floors themselves: it bounds a
known artefact on the oracle's side, at two geometries out of eighty-eight, and
is not a statement about this engine's accuracy. Away from the branch the two
codes agree to about 5e-9 kcal/mol.

This is the one place in the file where the oracle is the less accurate of the
two, which is worth stating plainly: an oracle is a second implementation, not a
definition.

### 22. MOPAC's SCF is the limit unless you tighten it past `PRECISE`

CS2 under INDO disagreed with xndo-rs by 3.4e-4 eV on every orbital, spread
evenly across the spectrum rather than concentrated anywhere -- the shape of a
slightly different density, not a wrong term.

It was MOPAC's own convergence. Tightening MOPAC's SCF moves **MOPAC's own**
eigenvalues by 3.8e-4 eV, the same size as the disagreement:

| MOPAC setting | worst `|dEps|` vs xndo-rs |
| --- | ---: |
| `PRECISE` | 3.429e-4 eV |
| `PRECISE RELSCF=0.0001` | 1.050e-5 eV |
| `PRECISE RELSCF=0.000001` | 1.050e-5 eV |

1e-6 gains nothing over 1e-4, so `RELSCF=0.0001` is where MOPAC's SCF stops
being the limit and its print format starts -- 1.05e-5 eV is two units in the
last of the five decimals it prints. Every reference file is now generated with
it (`run_mopac.SCF_TIGHTENING`).

The effect was not confined to INDO. Regenerating all four sets moved MNDO's
charges from 4.9e-5 to 3.6e-6 e, MNDO/d's frontier orbitals from 1.3e-4 to
1.1e-5 eV, and the CIS excitation energies from 4.6e-4 to 1.9e-5 eV.

### 23. The reference geometry must be what the oracle was given, to the digit

`to_xyz_field` wrote the geometry into the reference file at six decimals while
`run_mopac.write_mop` handed MOPAC eight. The two codes were therefore run on
geometries up to 5e-7 A apart, and every comparison carried a geometry-mismatch
error underneath whatever it was meant to measure.

On methane that mismatch was the *whole* of the reported disagreement: at
matched geometries MNDO and xndo-rs agree to 5e-9 kcal/mol, against the 3.2e-5
the suite had been reporting. Silane likewise, 1.2e-4 down to 3e-9.

The sharpest case was BF3's excited-state dipole. A D3h geometry rounded to six
decimals is not quite D3h, the broken degeneracy redistributed the dipole across
a near-degenerate triple, and the suite reported a 2.2e-3 D disagreement on a
state whose dipole vanishes by symmetry. At eight decimals it is 4.5e-6 D.

The column now carries eight decimals, matching `to_xyz_text` and `write_mop`
exactly, so the reference geometry is bit-for-bit the one MOPAC saw.

### 24. CNDO/2 and INDO against MolDS

xndo-rs's CNDO/2 and INDO engines are ports of MolDS 0.3.1, so MolDS is the
program those tables came from and the right oracle for them. Both suites are
enforced and complete:

| | molecules | elements | points | orbitals | charges | total E | core repulsion |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CNDO/2 | 17 | 6 | 270 | 1.34e-5 eV | 1.20e-7 e | 3.46e-4 eV | 4.35e-4 eV |
| INDO | 14 | 5 | 208 | 1.78e-5 eV | 1.22e-7 e | 4.22e-4 eV | 4.35e-4 eV |

**Building MolDS.** It is 2012 C++ and was the riskiest step in the whole oracle
plan. Three things needed handling that its shipped `Makefile_GNU` does not
anticipate:

* g++ 13 defaults to C++17 and MolDS predates C++11, while Boost 1.83 requires
  at least C++11. `-std=c++14 -fpermissive` is new enough for Boost and lenient
  enough for MolDS.
* Debian's `liblapacke.so` has undefined references to the LAPACK *testing*
  routines (`clagsy_`, `slatms_`, ...), which live in `libtmglib` and are no
  longer packaged separately. Linking the static `liblapacke.a` pulls in only
  the objects MolDS actually calls, so they never come up.
* `apt-get install` needs root, but `apt-get download` and `dpkg -x` do not, so
  LAPACKE can be unpacked into a private prefix without touching the system.

**The build earns its place before its numbers are used.** The five CNDO/2 and
INDO regression cases MolDS ships are vendored under `third_party/molds/tests/`
and re-run against this binary: all five come back **byte identical** in every
orbital energy, charge, energy and dipole. That is the same check the upstream
harness performs, and it is what separates "it compiled" from "it computes what
its author computed". Those files stay in the tree for exactly that purpose;
they are no longer the reference data.

**The reference sets are our own geometries.** Drawn from the same table the
MOPAC suites use, so a molecule appearing in both is the same geometry in both
and a disagreement between the two oracles is about the model rather than the
coordinates. `rms_density` is tightened to 1e-9, well past the 1e-6 the bundled
cases use -- the same lesson as item 22.

**MolDS has its own constants.** `src/base/Parameters.cpp:45` sets
`eV2AU = 0.03674903`, a Hartree of **27.211602592 eV**, 7.95e-6 above CODATA.
Invisible on an orbital energy, but it multiplies the core repulsion, which runs
to hundreds of eV:

| molecule | core repulsion | predicted offset | measured |
| --- | ---: | ---: | ---: |
| CH4 | 264.5165 eV | 2.103e-3 eV | -- |
| C2H6 | ~723 eV | 5.748e-3 eV | 5.748e-3 eV |

`ModelConstants::MOLDS` applies it the way MNDO/d's set is applied (item 13). It
took the core repulsion from 5.75e-3 to 2.62e-5 eV on the bundled cases.

**What sets the floor now.** MolDS prints seven significant figures, so
benzene's 2803.807 eV core repulsion carries 5e-4 eV of rounding; the observed
4.35e-4 is inside one unit of that. No keyword widens the format, so the energy
floors are the oracle's output format rather than either program's arithmetic.
The orbital and charge floors are still the model's, at about 1e-5 eV and 1e-7 e.

**Two things are deliberately not compared.**

*Ions.* MolDS's input deck has no total-charge keyword; every calculation it
runs is neutral. `requires_ions` is therefore false for these two suites, which
records the oracle's limitation rather than weakening the check everywhere.

*The dipole.* `Cndo2::CalcElectronicTransitionDipoleMoment`
(`src/cndo/Cndo2.cpp:1796`) contracts the density with a `cartesianMatrix` over
the real Slater AOs expanded in Gaussians -- two-centre elements included -- and
with the AO overlap matrix. xndo-rs uses the NDDO form: atom-centred point
charges plus the one-centre s-p hybridisation term
`<ns|r|np> = (2n+1)/(2 zeta sqrt(3))`, which is what MOPAC's `dipol` and the
ZINDO/S engine here use, and what every other suite in this file compares.

Those are two different definitions of the dipole in a zero-differential-overlap
model, not a right and a wrong one, and on CH4 they differ by 0.23 D. Matching
MolDS needs contracted-Gaussian position integrals, which is the same machinery
the Molden export needs, so the two belong together. Until then these suites
compare everything **except** the dipole, and `CndoIndoResult::dipole_debye` is
documented as the NDDO form.

### 25. An open shell's unpaired electron was being assigned by the guess

Found while checking the new HOMO/LUMO surface against MOPAC7. Planar methyl
under MINDO/3 came out **3.07 eV above** MOPAC7, with an ionisation potential of
6.609 eV against the oracle's 9.656 and a spin population on the carbon `p_z` of
exactly **0.0000** where MOPAC7 has **1.00000**.

It was not a Hamiltonian error. The ladder below localises it: the atoms and the
diatomic are exact, closed-shell methane is at the usual residual, and only the
planar open-shell case moves.

| | MOPAC7 | xndo-rs, before | after |
| --- | ---: | ---: | ---: |
| H atom | -12.50500 | 0 | 0 |
| C atom, triplet | -119.47000 | 0 | 0 |
| CH, doublet | -135.66445 | 1.9e-4 | 1.9e-4 |
| **methyl, doublet** | **-169.36017** | **3.07e+0** | **6.8e-4** |
| methane | -186.19578 | 8.7e-4 | 8.7e-4 |

Pyramidalising the methyl by 0.05 A -- which lets the out-of-plane p mix and
destroys the degeneracy the trap needs -- made the *unchanged* code agree to
6.7e-4 eV. That is what ruled out the Hamiltonian and pointed at the SCF.

**What was happening.** Which orbital an unpaired electron occupies is decided by
the aufbau ordering of the **first** Fock matrix, built from the guess -- the
worst density of the whole run -- and once the electron is placed,
self-consistency defends the choice. The SAD guess is the *free* atom, so carbon
arrives with 2.0 electrons in its 2s where a bonded carbon has about 1.34
(MOPAC7's own number for this molecule). With a one-centre `gsp` of 11.47 eV that
surplus lifts the carbon `p_z` several eV, and at iteration 1 it sat at -3.091 eV
against the C-H antibonding `a1*` at -3.774. The odd electron went into `a1*` and
stayed.

The result was a fully converged, aufbau-satisfying, perfectly real SCF solution
-- just not the lowest one. **Every accelerator and every damping factor reached
it**, A-DIIS, CDIIS, none, and damping from 0.0 to 0.9, all to the same 15
digits. It is not a convergence failure and cannot be fixed by converging harder.

**The fix, in `scf_accel`.** Two parts, both shared by all four engines:

* `split_by_spin` now puts the excess spin where an AO has *room* for it,
  weighting by `min(P, 2-P)`, instead of in proportion to `P` -- which was
  backwards, giving the largest share of the unpaired spin to the *most*
  occupied AO. The rule has no free parameters and is not fitted to this case:
  applied to a free atom it reproduces Hund's rule exactly (C `s(1,1) p(2/3,0)`,
  N `s(1,1) p(1,0)`, O `s(1,1) p(1,1/3)`, F `s(1,1) p(1,2/3)`).
* `open_shell_starts` + `lowest_solution`: an open shell is run from the SAD
  start **and** from the equal-valence-occupancy start at the other end of the
  same axis, and the variationally lower fixed point is kept. A closed shell
  still gets exactly one start, so nothing restricted changes.

Neither part alone is enough -- the improved spin split still cannot separate
`p_x`, `p_y` and `p_z`, which the geometry distinguishes and the guess does not.

**Every oracle point in this directory is unchanged**, to the digit, in all six
suites. A second start can only ever lower the energy, so this is monotone by
construction: it fixes cases that were wrong and cannot disturb cases that were
right.

While fixing it, the ZINDO/S **UHF** loop turned out to have no accelerator at
all -- the one place item 20's "A-DIIS/CDIIS everywhere" did not hold -- and its
own inline proportional spin split. Both now go through `scf_accel` like the
other three.

### 26. MINDO/3 against MOPAC7, and the two defects the first run found

MINDO/3 was the last method without a suite, and the only one MOPAC 23.2.5
cannot be the oracle for: `molkst_C.F90` dispatches twenty methods and MINDO/3
is not among them. MOPAC7 1.15 is the oracle instead, and is the declared
upstream of `src/data/mindo3_parameters.csv` rather than a proxy for it.

| | molecules | elements | points | HoF | charges | orbitals | total E | core |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| MINDO/3 | 36 | 10 | 611 | 6.0e-5 kcal | 5.9e-5 e | 9.2e-5 eV | 6.5e-6 eV | 7.7e-6 eV |

All five floors are MOPAC7's **printing**, not either program's arithmetic: it
prints charges and eigenvalues to four decimals and the energies to five, and no
keyword widens the format.

**Two real defects, both on this side, both fixed rather than tolerated.**

*`hp2` was stored rounded* for four of nine elements. Item 4 has the detail; the
short version is that MOPAC7 computes the p-p' exchange integral as
`(GPP - GP2)/2` inside `fock1.f` and stores nothing, so the relation is
structural. Cost: 2.2e-2 eV per chlorine atom, 8.6e-2 eV on CCl4.

*The Slater exponents were never restated in MOPAC7's Bohr.* This is item 13's
defect, in the one engine it had not been carried to. The overlap depends on
`zeta * R` and nothing else, so exponents fitted against a0 = 0.529167 A
(`analyt.f:47`, `delri.f:22`) and a distance in CODATA Bohr disagree by 1.93e-5
in every exponential in the integral. It multiplies the resonance energy, which
is several eV per bond, so it showed up as a few times 1e-4 eV **per atom pair**
-- benzene 4.5e-3 eV, methane 8.7e-4, and a constant floor under every molecule
in the set. `params.rs` and `cndo_indo.rs` had applied the scaling for years;
`mindo3.rs` had not.

What made the second one findable was comparing the core-core repulsion
separately from the total. MOPAC7 prints the total energy, the electronic energy
and the core repulsion as three numbers, which MOPAC 23 does not, and the core
term agreed to 7.7e-6 eV throughout. Since `gamma` enters the core repulsion
with a `Z_A Z_B` prefactor -- 49 for Cl2 -- that agreement bounds `gamma` to
about 1.6e-7 relative and leaves only the resonance term. The residual scaled
with atom pairs rather than with bonds or atoms, which is what the two-centre
overlap does.

CCl4 moved by five orders of magnitude between the two fixes.

**What the oracle cannot do, recorded rather than worked around.**

*MOPAC7 assumes internal coordinates for three atoms or fewer* (`getgeo.f:256`:
`IF(NATOMS.GT.3) THEN INT=(NA(4).NE.0) ELSE INT=.TRUE.`), and no keyword
overrides it -- `XYZ` is consulted later and does something else. Handing it
Cartesians there does not fail, which is the dangerous part: water came back as
a straight O-H-H chain with a 0.7575 A bond, converged, and reported an energy.
`run_mopac7.write_input` writes a Z-matrix below four atoms and Cartesians at or
above it. A Z-matrix that small never needs a dihedral, so the usual
linear-molecule difficulty does not arise.

*MOPAC7 rebuilds its own frame from a Z-matrix*, so for those molecules the
dipole components are in MOPAC7's orientation. This suite therefore compares the
dipole **magnitude**, which MOPAC7 already prints beside them and which is
invariant. Weaker than three components, and the honest thing to compare.

*MOPAC7 rejects Cartesian input whose first three atoms are collinear*
(`getgeo.f:373`). Which three atoms are written first is not a property of the
molecule, so the generator reorders one non-collinear atom into third place and
writes the reordered list into the reference file, leaving both programs the
same atom order. Ethyne is the one molecule no ordering helps, and it is
declared rather than worked around.

*MOPAC7 stops on any interatomic distance below 0.8 A* unless `GEO-OK` is given
(`moldat.f:638`). H2's bond is 0.741 A. `GEO-OK` is added only where the guard
would actually fire, so it stays armed for every other molecule; switching it on
unconditionally would also switch off the one check that catches a mistyped
coordinate.

*MINDO/3's pair table covers 40 of the 55 pairs its ten elements could form.*
F-S is one of the fifteen it does not, so SF6 is outside the published model
rather than a case this engine fails; `mindo3_molecules` reads the shipped table
and filters on it, rather than listing exclusions by hand.

### 15. Element coverage is currently partial, and is not yet gated

The MNDO reference set covers 19 elements (Al B Be Br C Cl F Ge H I Li Mg N Na O
P S Si Zn) of the 57 that have usable MNDO parameters -- though see item 16: eight
of the remaining ones are transition metals MOPAC will not run under MNDO at all.
MNDO/d covers 17 of 22. The
plan's CSV-driven coverage assertion — every parameterised element appears in the
set, or is listed with a reason — is **not** implemented yet, because asserting
it today would mean writing 39 exemptions, which documents nothing. The suite
asserts the weaker properties that are meaningful now: at least ten elements, at
least one open-shell row, at least one cation and one anion.

Expanding the sets to full element coverage is remaining work, and that is when
the CSV-driven gate should go in.

---

## (d) How the files are produced and checked

Three layers, so that a checkout without the oracle binaries still detects drift:

1. `python tools/oracle/build_reference_set.py --check` re-runs the oracle and
   byte-compares the rendered text. Requires the oracle binary.
2. `tests/data/MANIFEST.sha256` records the SHA-256 of each TSV **and of the
   generator script that produced it**, plus the oracle version string. Checked
   by a plain `#[test]`; needs no oracle.
3. Each TSV carries a `#` provenance header (`oracle_program`, `oracle_version`,
   `generator`, `generator_sha256`, `generated_utc`) which the harness asserts
   against the version it expects, so swapping MOPAC 23.2.5 for a later release
   fails loudly.

Layer 1 says something only when the oracle ran. With none of the three binaries
installed every molecule is skipped, and a file rendered from no rows differs
from the committed one in every line: reporting that as `STALE` would read as
"the shipped data is wrong" when the truth is that nothing was compared. It is
reported as `CANNOT CHECK` and does not fail. A run that produces *some* but not
all of the committed rows is `INCOMPLETE` and does fail, because a half-working
oracle environment is a thing to repair, not to draw conclusions from.

The `generator_sha256` in a TSV header names the generator revision that
produced those numbers, which is a historical fact and does not change when the
generator is later edited. It therefore need not equal the hash the manifest
records for the current `build_reference_set.py`, and as of v0.3.0 it does not:
the generator was edited after the last oracle run, in `main()` and in the new
`emit()` reporting helper only. `render()` and every `compute_*` are untouched,
so the bytes a given oracle run produces are unchanged. An edit that *can* change
those bytes obliges a regeneration, and the manifest is what makes it impossible
to make one quietly.

The rule: **no number in `tests/data/` is ever typed by a human.** The single
exception is literature-transcribed rows, which follow the two-script protocol —
one script extracts the literals, a second, independently written script
re-parses both sides and compares.

---

## (e) Validation outcome

### MNDO - enforced

51 molecules, **295 independent points**, against OpenMOPAC 23.2.5.

| quantity | comparisons | worst deviation | at | floor |
| --- | ---: | ---: | --- | ---: |
| heat of formation | 51 | **1.227e-3 kcal/mol** | germane | 2.0e-3 |
| atomic charge | 176 | **4.992e-5 e** | ch2_triplet | 5.0e-4 |
| frontier orbital | 88 | **4.350e-5 eV** | zncl2 LUMO | 1.5e-3 |
| dipole component | 153 | **6.969e-5 D** | licl | 1.5e-3 |

### MNDO/d - enforced

37 molecules, **222 independent points**.

| quantity | comparisons | worst deviation | at | floor |
| --- | ---: | ---: | --- | ---: |
| heat of formation | 37 | **7.421e-2 kcal/mol** | benzene | 1.0e-1 |
| atomic charge | 136 | **3.674e-4 e** | silane | 5.0e-4 |
| frontier orbital | 66 | **9.095e-4 eV** | sf6 LUMO | 1.5e-3 |
| dipole component | 111 | **9.489e-4 D** | so2 | 1.5e-3 |

The heat-of-formation floor is an order of magnitude looser than MNDO's, and
deliberately so: it is bounding open item 13, a **named and understood** model
difference whose residual is per-pair and grows with size (2.25e-3 kcal/mol for
H2, 7.42e-2 for benzene). It is not absorbing an unexplained disagreement. When
item 13 is resolved this floor comes down to MNDO's.

Every floor in both tables is set from the measurement, a little above it,
rather than from what would be nice. Raising one requires the measurement that
justifies it in the same commit, naming the deviation bucket.

### ZINDO/S - enforced

39 molecules compared (benzene and bh are in `KNOWN`, items 18 and 19), **418
independent points**.

| quantity | comparisons | worst deviation | at | floor |
| --- | ---: | ---: | --- | ---: |
| orbital spectrum | 332 | **3.412e-4 eV** | cs2 | 1.0e-3 |
| atomic charge | 125 | **1.336e-4 e** | cs2 | 5.0e-4 |

The orbital-spectrum check compares *every* eigenvalue the oracle printed, not
just the frontier pair, which is the strongest convention-free comparison
available for a method with no atomic heat terms. That this agrees to 3.4e-4 eV
across 39 molecules is what establishes that the two ZINDO/S defects found
earlier (items 3 and 3b) are fixed.

Neither `KNOWN` entry is an accepted disagreement. Both are open defects with a
named cause, recorded so the suite states them rather than hiding them, and both
are to be fixed at the root.

### Deviation buckets

1. **Convention difference** (units, sign, frame, printed precision, active-space
   definition). Fix the *comparison*, never the tolerance. Document here with
   both `file:line` references.
2. **Different SCF solution** — same functional, different fixed point. Report
   which is variationally lower, register in `KNOWN` with the observation that
   distinguishes it. Never loosen.
3. **Model difference** — a term one program has and the other does not. A bug
   or a scope decision, and resolved as one. Items 3, 13 and 14 are all bucket 3.
