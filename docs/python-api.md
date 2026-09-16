# Python and ASE API

This document describes every public Python function and ASE calculator in
xndo-rs 0.3.0. The thin wrappers live in `xndo_rs.native`; the same functions
are re-exported at package level, so `xndo_rs.gradient(...)` and
`xndo_rs.native.gradient(...)` are equivalent.

## Installation and imports

Build a wheel from the source tree with:

```text
maturin build --release --features python
```

Install ASE support with `xndo-rs-python[ase]`. The package requires Python
3.9 or newer and NumPy. The version is available as `xndo_rs.__version__`.

```python
import xndo_rs
from xndo_rs import native
```

## Common conventions

All `numbers` arguments are one-dimensional atomic-number sequences of length
`N`. All `positions` arguments have shape `(N, 3)` and are in Angstrom. A
non-integer `charge` is accepted by the API, although ordinary molecular use
normally supplies an integer. `multiplicity` is `2S + 1` and must be positive.

`reference` accepts:

- `"auto"`: RHF for a closed shell and UHF for an open shell;
- `"rhf"`: force a restricted reference; incompatible electron counts fail;
- `"uhf"`: force an unrestricted reference, including for a singlet.

Method parsing is case-insensitive. Spaces, underscores, and hyphens are
normalized where applicable. The exact enabled strings are:

| Method | Accepted strings | Ground derivatives | Optimization | Excited states |
|---|---|---:|---:|---:|
| CNDO/2 | `cndo2`, `cndo/2`, `cndo` | yes | yes | no |
| INDO | `indo` | yes | yes | no |
| MNDO | `mndo` | yes | yes | no |
| MNDO/d | `mndod`, `mndo/d`, `mndo-d` | yes | yes | no |
| MINDO/3 | `mindo3`, `mindo/3`, `mindo` | yes | yes | no |
| ZINDO/S | `zindo/s`, `zindos`, `zindo`, `indo/s`, `indos` | yes | yes | yes |

Bare `indo` selects the
ground-state INDO engine. Use `indo/s` for ZINDO/S. Registered but
non-executable historical labels are returned by
`available_methods()` and fail explicitly if passed to a computational API.
Use `api_methods()` when method availability must be discovered at runtime.

Output dictionaries use unit-bearing key names. Arrays returned by native
functions are Python lists; ASE adapters convert numerical results to NumPy
arrays. Cartesian gradients and forces have shape `(N, 3)`, and Hessians have
shape `(3N, 3N)`.

## Ground-state functions

### `single_point`

```python
single_point(
    numbers,
    positions,
    charge=0.0,
    multiplicity=1,
    reference="auto",
    method="mndo",
) -> dict
```

Runs a fixed-geometry SCF calculation for any native method. Every result has:

| Key | Meaning |
|---|---|
| `method` | Canonical method name. |
| `energy_hartree`, `energy_ev` | Total electronic-plus-core energy. |
| `electronic_ev`, `core_ev` | Electronic and core-repulsion contributions. |
| `charges` | Net atomic NDO charges, length `N`. |
| `unrestricted` | Whether the UHF path was used. |
| `spin_density` | AO spin-density matrix for UHF, otherwise `None`. |
| `iterations` | SCF iterations taken. |
| `converged` | SCF convergence flag. Non-convergence normally raises before a result is returned. |
| `mo_energies_ev`, `mo_energies_beta_ev` | Orbital energies in eV, ascending. Beta is `None` for RHF. |
| `n_occ`, `n_alpha`, `n_beta` | Occupied counts; `n_occ` is `n_alpha`. |
| `occupations`, `occupations_beta` | Per-orbital occupation numbers (2.0/0.0 restricted, 1.0/0.0 unrestricted). `occupations_beta` is `None` for RHF. |
| `homo_ev`, `lumo_ev`, `homo_lumo_gap_ev` | Frontier pair over **both** spin channels, and the gap. |
| `homo_alpha_ev`, `lumo_alpha_ev`, `homo_beta_ev`, `lumo_beta_ev` | The same per channel. |
| `dipole_debye` | Permanent dipole components. Not returned by MINDO/3. |

MNDO, MNDO/d and MINDO/3 additionally return `heat_of_formation_kcal`; the
other three have no atomic heat terms, so their total energy is the comparable
quantity.

Any frontier value that does not exist -- a full or an empty shell -- is `None`
rather than a substituted number, and so is the gap that would need it.

`homo_ev` and `lumo_ev` are read over both spin channels, so for an open-shell
doublet the LUMO is usually the beta partner of the singly occupied orbital
rather than the lowest unoccupied *alpha* orbital. Use `lumo_alpha_ev` when the
alpha diagram alone is what you want. Alpha and beta eigenvalues come from
different Fock operators and are not levels of one orbital diagram.

### `orbital_energies`

```python
orbital_energies(
    numbers,
    positions,
    charge=0.0,
    multiplicity=1,
    reference="auto",
    method="mndo",
) -> dict
```

The orbital block of `single_point` on its own, plus `method` and `energy_ev`,
for any method with a native engine. Use it when the orbitals are what you want
and the rest of the result is not.

```python
import xndo_rs
o = xndo_rs.orbital_energies([6, 1, 1, 1], methyl_xyz, multiplicity=2, method="mindo3")
print(o["homo_ev"], o["lumo_ev"], o["homo_lumo_gap_ev"])
```

### `molden`

```python
molden(
    numbers,
    positions,
    charge=0.0,
    multiplicity=1,
    reference="auto",
    method="mndo",
    coefficients="lowdin",
) -> str
```

A Molden wavefunction file, as a string. The Slater basis is expanded as STO-6G
(Stewart, *J. Chem. Phys.* **52**, 431 (1970)) and d shells are written in
Molden's `[5D]` order.

The MO coefficients are back-transformed with `S^(-1/2)`, and that is not
cosmetic. Every engine here assumes an orthonormal AO basis -- that is what zero
differential overlap means -- while a Molden file describes real Gaussians,
which are not orthonormal. A reader given the raw coefficients forms
`P = C n C^T` over a non-orthogonal basis and gets a density that does not
integrate to the electron count and orbitals that are not normalised.

`coefficients="raw"` writes the untransformed coefficients, for comparison with
programs that make that identification. Such a file says so in its own title.

ZINDO/S rejects elements that would need its unimplemented d branch, rather
than writing a file missing their d coefficients.

```python
import xndo_rs
text = xndo_rs.molden([8, 1, 1], water_positions, method="mndo")
open("water.molden", "w").write(text)
```

`XNDO.write_molden(path)` does the same from the ASE calculator.

### `gradient`

```python
gradient(
    numbers,
    positions,
    charge=0.0,
    multiplicity=1,
    reference="auto",
    method="mndo",
) -> dict
```

Returns an analytic ground-state Cartesian gradient for every native method.
Keys are `method`, `energy_hartree`, `energy_ev`, `iterations`,
`unrestricted`, `gradient_hartree_per_bohr`, and
`gradient_ev_per_angstrom`. MNDO, MNDO/d, and MINDO/3 also return
`heat_of_formation_kcal`.

### `forces`

```python
forces(
    numbers,
    positions,
    charge=0.0,
    multiplicity=1,
    reference="auto",
    method="mndo",
) -> dict
```

Uses the analytic gradient and returns its negative. Keys are `method`,
`energy_hartree`, `energy_ev`, `iterations`, `unrestricted`,
`forces_hartree_per_bohr`, and `forces_ev_per_angstrom`. MNDO, MNDO/d, and
MINDO/3 also return `heat_of_formation_kcal`.

### `hessian`

```python
hessian(
    numbers,
    positions,
    charge=0.0,
    multiplicity=1,
    reference="auto",
    method="mndo",
) -> dict
```

Returns the analytic Cartesian Hessian for every native method. The result has
`method`, `hessian_hartree_per_bohr2`, and `ndof`, where `ndof == 3 * N`.
The Hessian ordering is `x, y, z` for atom 0, then atom 1, and so on.

### `frequencies`

```python
frequencies(
    numbers,
    positions,
    charge=0.0,
    multiplicity=1,
    reference="auto",
    method="mndo",
) -> dict
```

Mass-weights and diagonalizes the analytic ground-state Hessian. The result
has `method`, `frequencies_cm`, and `eigenvalues`. Frequencies are in
`cm^-1`, sorted ascending; negative values represent imaginary modes.
`eigenvalues` are the mass-weighted Hessian eigenvalues in
`eV / (Angstrom^2 amu)`. The routine returns all `3N` Cartesian modes and does
not project translations or rotations.

### `optimize`

```python
optimize(
    numbers,
    positions,
    charge=0.0,
    multiplicity=1,
    reference="auto",
    method="mndo",
    max_iter=200,
    gtol=1.0e-3,
    history=8,
) -> dict
```

Runs the unified L-BFGS optimizer for CNDO/2, INDO, MNDO, MNDO/d, MINDO/3,
or ZINDO/S. `gtol` is convergence on the largest absolute Cartesian gradient
component in eV/Bohr. `max_iter` limits accepted L-BFGS iterations and
`history` selects the number of stored correction pairs. All three must be
positive; `gtol` must also be finite. The result keys are:

| Key | Meaning |
|---|---|
| `method` | Canonical native-method name. |
| `positions_angstrom` | Final `(N, 3)` geometry. |
| `energy_hartree`, `energy_ev` | Final total energy. |
| `heat_of_formation_kcal` | Final heat of formation for MNDO, MNDO/d, and MINDO/3; absent otherwise. |
| `converged` | Whether the gradient threshold was reached. |
| `iterations` | Number of L-BFGS iterations. |

If the line search cannot find an energy-lowering Armijo step, optimization
returns the last accepted geometry with `converged == False`. Numerical and
SCF failures still raise `ValueError`.

## ZINDO/S excited-state functions

The excited-state APIs execute only the ordinary one-/four-AO s/p ZINDO/S
branch. Transition-metal nine-AO ZINDO/S is rejected until its full
method-specific Hamiltonian is available. `active_occupied` retains the
requested number of occupied orbitals nearest the HOMO; `active_virtual`
retains virtual orbitals nearest the LUMO. `None` means the full corresponding
space.

For an RHF reference, `state_type` accepts `"singlet"`, `"triplet"`, or
`"both"`. For UHF, the API selects spin-orbital UCIS and reports
`"unrestricted"`, regardless of the supplied RHF `state_type` value.

### `excited_states` and `excited_properties`

```python
excited_states(
    numbers,
    positions,
    charge=0.0,
    n_states=10,
    active_occupied=None,
    active_virtual=None,
    method="zindo/s",
    state_type="singlet",
    multiplicity=1,
    reference="auto",
) -> dict

excited_properties(...) -> dict
```

`excited_properties` is an exact alias that emphasizes the property-rich
output. The outer result has `method`, `reference`, `state_type`,
`ground_energy_hartree`, `ground_energy_ev`, `ground_charges`,
`ground_dipole_au`, `ground_dipole_debye`, `mo_energies_ev`,
`mo_energies_beta_ev`, and `states`.

Every entry in `states` has:

| Key | Meaning |
|---|---|
| `state` | One-based root number within the returned spin sector. |
| `spin` | `singlet`, `triplet`, or `unrestricted`. |
| `spin_multiplicity`, `s2_expectation` | Spin labels (`1, 0` for singlet; `3, 2` for triplet). |
| `energy_ev`, `energy_cm1`, `wavelength_nm` | Vertical excitation in three representations. |
| `state_total_energy_ev`, `state_total_energy_hartree` | Ground energy plus vertical excitation. |
| `oscillator_strength` | Electric-dipole oscillator strength; spin-free triplets are zero. |
| `transition_dipole_au`, `transition_dipole_debye` | Transition-dipole vector. |
| `transition_dipole_magnitude_au`, `dipole_strength_au2` | Transition-dipole scalar measures. |
| `permanent_dipole_au`, `permanent_dipole_debye`, `permanent_dipole_magnitude_debye` | CIS-state permanent dipole. |
| `difference_dipole_au`, `difference_dipole_debye`, `difference_dipole_magnitude_debye` | Excited-minus-ground dipole. |
| `charges` | State-specific net atomic charges. |
| `hole_population`, `electron_population` | Per-atom transition populations, each summing to approximately one. |
| `hole_centroid_angstrom`, `electron_centroid_angstrom` | Population centroids. |
| `charge_transfer_distance_angstrom` | Distance from hole to electron centroid. |
| `dominant_configurations` | `(occupied, virtual, coefficient)` tuples with one-based MO indices. |

These are vertical state properties. No geometry relaxation, line broadening,
solvent shift, or molar-absorptivity convention is applied.

### `uv_vis_spectrum`

```python
uv_vis_spectrum(
    numbers,
    positions,
    charge=0.0,
    n_states=10,
    active_occupied=None,
    active_virtual=None,
    method="zindo/s",
    multiplicity=1,
    reference="auto",
) -> dict
```

Returns the same dictionary schema as `excited_states`. RHF uses the singlet
CIS sector; UHF uses UCIS. The output is a stick spectrum without an imposed
line shape.

### `excited_state_gradients`

```python
excited_state_gradients(
    numbers,
    positions,
    charge=0.0,
    n_states=10,
    active_occupied=None,
    active_virtual=None,
    state_type="singlet",
    multiplicity=1,
    reference="auto",
) -> dict
```

ZINDO/S is implicit. The outer result contains `method`, `reference`,
`ground_energy_ev`, and `states`. Each state contains `root`, `spin`,
`excitation_energy_ev`, `state_total_energy_ev`,
`excitation_gradient_hartree_per_bohr`,
`excitation_gradient_ev_per_angstrom`,
`state_gradient_hartree_per_bohr`, and
`state_gradient_ev_per_angstrom`. The excitation gradient differentiates the
vertical gap; the state gradient differentiates ground energy plus that gap.

### `excited_state_hessians`

```python
excited_state_hessians(
    numbers,
    positions,
    charge=0.0,
    n_states=10,
    active_occupied=None,
    active_virtual=None,
    state_type="singlet",
    multiplicity=1,
    reference="auto",
) -> dict
```

ZINDO/S is implicit. The outer result contains `method`, `reference`, and
`states`. Each state contains `root`, `spin`, `excitation_energy_ev`,
`state_total_energy_ev`, `excitation_hessian_hartree_per_bohr2`,
`excitation_hessian_ev_per_angstrom2`,
`state_hessian_hartree_per_bohr2`, and
`state_hessian_ev_per_angstrom2`.

## Discovery and parameter-data functions

### `available_methods`

```python
available_methods() -> list[dict]
```

Returns one row for each registered `Method`, including executable and
non-executable labels. Every row has `name`, `family`, `status` (`native`,
`registered`, or `noncanonical`), boolean `energy`, `gradient`, `hessian`,
`uhf`, `spectrum`, and `excited_properties` capabilities,
`accepted_strings`, enabled `apis`, and an implementation `note`.

### `api_methods`

```python
api_methods() -> dict[str, list[dict]]
```

Maps each computational API name to the methods it accepts. Each row contains
`method`, `accepted_strings`, `uhf`, `reference_argument`,
`reference_strings`, `fixed_reference`, and `state_type_strings`. This is the
preferred compatibility-safe way to populate a UI or validate a user method
string.

### `parameter_datasets`

```python
parameter_datasets() -> list[dict]
```

Parses and returns the bundled provenance catalog. Dictionary keys are the
columns of `legacy_parameter_catalog.csv`, including dataset identifier,
upstream source, license, and runtime role.

### `parameter_dataset`

```python
parameter_dataset(name: str) -> str
```

Returns one provenance-preserving embedded CSV or manifest as text. Canonical
dataset names and useful aliases include `mndo`, `mndo_pair`, `mndod`,
`mndod_pair`, `zindo_s`, `mindo3`, `mindo3_pair`, `molds_cndo2_indo`,
`molds_zindo_s`, `catalog`, and `manifest`. Unknown names raise `ValueError`;
call `parameter_datasets()` to discover canonical IDs.

### `third_party_licenses`

```python
third_party_licenses() -> list[dict]
```

The licence and attribution documents this build embeds, one dict per document
with `path`, `role` and the verbatim `text`. These are the notices Apache-2.0
4(c) and GPL-3.0 5(a) require to be carried, so they travel with the installed
wheel rather than only with the source tree. The same documents are on disk in
`.dist-info/licenses/` and printable from the CLI with `xndo_rs_cli licenses`.

## ASE calculators

### Exports and naming

The ASE module's star-import surface is:

```python
xndo_rs.ase.__all__ == [
    "CNDO2", "INDO", "MNDO", "MNDOD", "MINDO3", "ZINDOS"
]
```

`XNDO` is the generic calculator implementation and can be imported explicitly
with `from xndo_rs.ase import XNDO`. It is intentionally absent from
`__all__`: XNDO is the library name, not a semiempirical method name. `MNDO`
remains the exported method-specific convenience class and is equivalent to
`XNDO(method="mndo")`.

### `XNDO`

```python
XNDO(
    charge=0,
    multiplicity=1,
    reference="auto",
    method="mndo",
    **ase_calculator_kwargs,
)
```

This generic ground-state calculator accepts every native method in the table
above. It advertises `energy`, `forces`, `charges`, `dipole`,
`heat_of_formation_kcal`, and `hessian`. `dipole` and
`heat_of_formation_kcal` are populated only when the selected native
single-point result supplies them.

ASE result units are eV and Angstrom:

- `energy`: eV;
- `forces`: eV/Angstrom;
- `charges`: elementary-charge units;
- `dipole`: e Angstrom;
- `hessian`: eV/Angstrom^2;
- `heat_of_formation_kcal`: kcal/mol.

Additional methods are:

- `get_gradient(atoms=None) -> ndarray`: negative ASE forces, eV/Angstrom;
- `get_hessian(atoms=None) -> ndarray`: `(3N, 3N)`, eV/Angstrom^2;
- `get_frequencies(atoms=None) -> ndarray`: all `3N` modes, `cm^-1`.

The method-specific constructors `MNDO`, `MNDOD`, `CNDO2`, `INDO`, and
`MINDO3` accept the same keyword arguments except `method`, which they set to
`mndo`, `mndod`, `cndo2`, `indo`, and `mindo3`, respectively.

```python
from ase import Atoms
from xndo_rs.ase import XNDO

atoms = Atoms(
    "H2",
    positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]],
)
atoms.calc = XNDO(method="mndo")
energy_ev = atoms.get_potential_energy()
forces_ev_per_angstrom = atoms.get_forces()
```

### `ZINDOS`

```python
ZINDOS(
    charge=0,
    multiplicity=1,
    reference="auto",
    n_states=10,
    active_occupied=None,
    active_virtual=None,
    **ase_calculator_kwargs,
)
```

`ZINDOS` combines ZINDO/S ground-state ASE properties with RHF-CIS or UHF-UCIS
state properties. It advertises `energy`, `forces`, `charges`, `dipole`, and
`hessian`. Its additional methods are:

- `get_gradient(atoms=None)`: ground-state gradient in eV/Angstrom;
- `get_hessian(atoms=None)`: ground-state Hessian in eV/Angstrom^2;
- `get_frequencies(atoms=None)`: ground-state frequencies in `cm^-1`;
- `get_uv_vis_spectrum(atoms=None, n_states=None)`: stick spectrum and all
  state properties described above;
- `get_excited_properties(atoms=None, n_states=None, state_type="singlet")`:
  explicit spin-sector state properties;
- `get_excited_state_gradients(atoms=None, n_states=None,
  state_type="singlet")`: excitation and absolute-state gradients;
- `get_excited_state_hessians(atoms=None, n_states=None,
  state_type="singlet")`: excitation and absolute-state Hessians.

The spectrum calculated during an ordinary ASE `calculate()` call is cached in
`calculator.results["uv_vis_spectrum"]`. `get_uv_vis_spectrum()` reuses it only
when called for the current `Atoms` object and the configured `n_states`.

## Errors and state-uniqueness boundary

Native validation and numerical errors are exposed to Python as `ValueError`.
Typical causes are a malformed geometry, incompatible charge/multiplicity/RHF
combination, unsupported element, non-native method, SCF failure, memory-limit
guard, or singular linear-algebra problem. Missing ASE raises `ImportError`
when `xndo_rs.ase` is imported.

State-specific CIS/UCIS derivatives require each requested root to be isolated
from every other root by at least 0.0001 eV. Near-degenerate roots retain their
vertical spectra but state gradients and Hessians raise an explanatory
`ValueError`. Open-shell diatomic state derivatives are likewise non-unique
without a specified state-averaged or diabatic model. Their ground gradient
and vertical UCIS spectrum remain available; the ambiguous state-specific
derivative is rejected rather than assigned an arbitrary adiabatic basis.
