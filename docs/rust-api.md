# Rust API

This document describes every item re-exported from the `xndo_rs` crate root
in version 0.2.4. A final module map also identifies the lower-level public
modules. The crate has `#![forbid(unsafe_code)]`.

## Units and result convention

Rust molecular positions are stored in Bohr. `Molecule::from_xyz_str` and
`Molecule::from_xyz_file` read conventional Angstrom XYZ coordinates and
convert them to Bohr. Ground-state energies are eV. Unless a field name states
otherwise, Rust gradients are eV/Bohr and Hessians are eV/Bohr^2.

Every fallible API returns `xndo_rs::Result<T>`, an alias of
`std::result::Result<T, XndoError>`. Functions do not silently substitute one
semiempirical Hamiltonian for another: a registered but non-executable method
returns `XndoError::InvalidInput` with the method's implementation status.

## High-level example

```rust
use xndo_rs::{
    run_gradient, run_hessian, run_method, CalculationResult, Method,
    Molecule, NddoOptions, Reference,
};

let molecule = Molecule::from_xyz_str(
    "2\nhydrogen\nH 0 0 0\nH 0 0 0.74\n",
    0.0,
)?;
let options = NddoOptions {
    reference: Reference::Rhf,
    ..NddoOptions::default()
};

let single_point = run_method(&molecule, Method::Mndo, &options)?;
assert!(matches!(single_point, CalculationResult::Nddo(_)));

let gradient = run_gradient(&molecule, Method::Mndo, &options)?;
let hessian = run_hessian(&molecule, Method::Mndo, &options)?;
assert_eq!(gradient.gradient.len(), 2);
assert_eq!((hessian.hessian.rows, hessian.hessian.cols), (6, 6));

# Ok::<(), xndo_rs::XndoError>(())
```

## Unified dispatch

### `run_method`

```rust
pub fn run_method(
    molecule: &Molecule,
    method: Method,
    nddo_options: &NddoOptions,
) -> Result<CalculationResult>
```

Runs a native fixed-geometry ground/reference-state engine. It accepts CNDO/2,
INDO, MNDO, MNDO/d, MINDO/3, and ZINDO/S. `NddoOptions` is translated to the
method-specific options for non-NDDO engines. The result enum is:

```rust
pub enum CalculationResult {
    CndoIndo(CndoIndoResult),
    Nddo(NddoResult),
    Mindo3(Mindo3Result),
    ZindoS(ZindoResult),
}
```

ZINDO/S excited states are separate because they require a converged reference
plus CIS/UCIS configuration settings.

### `run_gradient`

```rust
pub fn run_gradient(
    molecule: &Molecule,
    method: Method,
    options: &NddoOptions,
) -> Result<MethodGradientResult>
```

Dispatches the analytic ground-state gradient for every native method.
`MethodGradientResult` contains:

| Field | Meaning |
|---|---|
| `method: Method` | Executed Hamiltonian. |
| `energy_ev: f64` | Total energy. |
| `gradient: Vec<Vec3>` | `dE/dR`, eV/Bohr. |
| `forces: Vec<Vec3>` | Negative gradient, eV/Bohr. |
| `iterations: usize` | SCF iterations. |
| `unrestricted: bool` | Whether UHF was used. |
| `heat_of_formation_kcal: Option<f64>` | Present for MNDO, MNDO/d, and MINDO/3. |

### `run_hessian`

```rust
pub fn run_hessian(
    molecule: &Molecule,
    method: Method,
    options: &NddoOptions,
) -> Result<MethodHessianResult>
```

Dispatches the analytic ground-state Cartesian Hessian for every native
method. `MethodHessianResult` has `method: Method` and a `(3N, 3N)`
`hessian: Matrix` in eV/Bohr^2.

### `run_nddo` and `run_nddo_with_parameters`

```rust
pub fn run_nddo(
    molecule: &Molecule,
    method: Method,
    options: &NddoOptions,
) -> Result<NddoResult>

pub fn run_nddo_with_parameters(
    molecule: &Molecule,
    params: &NddoParameters,
    options: &NddoOptions,
) -> Result<NddoResult>
```

`run_nddo` accepts only `Method::Mndo` and `Method::MndoD` and uses the cached
built-in table. `run_nddo_with_parameters` executes an explicitly supplied
MNDO-family parameter set. `NddoCalculator::new(params)` constructs a reusable
calculator with default options; `NddoCalculator::with_options(params,
options)` uses custom settings; `calculate(&self, molecule)` returns
`NddoResult`.

## Method discovery

### `Method` and `MethodStatus`

`Method` has 18 variants:

```text
Cndo1, Cndo2, Cndo3, Indo, Indo1, Indo2, Indo3,
Zindo1, Zindo2, Zindo3, ZindoS,
Mindo1, Mindo2, Mindo3, Sindo, Msindo, Mndo, MndoD
```

The default is `Method::Mndo`. `Method::ALL` contains all variants.
`Cndo2`, `Indo`, `Mndo`, `MndoD`, `Mindo3`, and `ZindoS` have
`MethodStatus::Native`. `Cndo3`, `Indo3`, and `Zindo3` are
`NonCanonical`; the remaining labels are `Registered` but not executable.

Public methods on `Method` are:

| Method | Purpose |
|---|---|
| `as_str()` | Canonical display name. |
| `family()` | Method-family label. |
| `status()` | `Native`, `Registered`, or `NonCanonical`. |
| `supports_energy()`, `supports_gradient()`, `supports_hessian()` | Ground-state capability flags. |
| `supports_uhf()` | UHF capability flag. |
| `supports_spectrum()`, `supports_excited_properties()` | Excited-state capability flags. |
| `accepted_strings()` | Explicit user-facing aliases. |
| `api_names()` | Computational APIs enabled for the method. |
| `implementation_note()` | Status and provenance summary. |
| `execution_error()` | Human-readable reason a non-native method cannot run. |
| `parse(s)` | Case-insensitive parser returning `Option<Method>`. |

`Method` also implements `Default`, `Display`, and `FromStr`. Hyphens, spaces,
and underscores are normalized. Bare `indo` is ground-state INDO; `indo/s` is
ZINDO/S.

`MethodStatus::as_str()` returns `native`, `registered`, or `noncanonical`.

### `methods_for_api`

```rust
pub fn methods_for_api(api: &str) -> Vec<Method>
```

Returns native methods whose `api_names()` contains the supplied public API
name, such as `single_point`, `gradient`, `hessian`, `excited_states`,
`excited_properties`, or `uv_vis_spectrum`. Unknown API names return an empty
vector.

## Molecular data and small linear-algebra types

### `Atom` and `Molecule`

```rust
pub struct Atom {
    pub z: u8,
    pub position: Vec3, // Bohr
}

pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub charge: f64,
    pub multiplicity: usize,
}
```

`Molecule::new(atoms)` creates a neutral singlet. Builder methods
`with_charge` and `with_multiplicity` update those values; multiplicity is
clamped to at least one. `len`, `is_empty`, and `total_nuclear_charge` inspect
the geometry. `from_xyz_str(text, charge)` and `from_xyz_file(path, charge)`
read Angstrom XYZ coordinates. `symbol_to_z` accepts symbols or atomic-number
text; `z_to_symbol` performs the reverse lookup.

### `Vec3` and `Mat3`

`Vec3` exposes public `x`, `y`, and `z` values plus `new`, `zero`, `dot`,
`cross`, `norm2`, `norm`, `normalized`, `to_array`, and indexed-component
`get`. It implements the ordinary vector addition, subtraction, scalar
multiplication/division, negation, and assignment operators used by the
library.

`Mat3` stores three public column vectors in `col`. `from_columns`, `zero`,
and `mul_vec` create and apply a 3-by-3 matrix.

### `Matrix`

`Matrix` is a row-major dense `f64` matrix with public `rows` and `cols` and
private storage. Its public methods are:

- `zeros(rows, cols)` and `identity(n)`;
- `from_row_major(rows, cols, data)`, which panics if the length is wrong;
- `as_slice()` and `as_mut_slice()`;
- `transpose()`;
- `matmul()` using global faer parallelism;
- `matmul_seq()`, `transpose_matmul_seq()`, and `matmul_transpose_seq()` for
  nested-parallel contexts;
- `leading_columns_gram(count, scale)` for `scale * C * C^T`;
- `frobenius_dot(other)` and `rms_difference(other)`.

It supports `(row, column)` indexing. The `linalg` module additionally exposes
`symmetric_eigen`, `solve_linear`, and `solve_linear_matrix`.

## SCF configuration and NDDO results

### `Reference` and `ScfAccelerator`

`Reference::{Auto, Rhf, Uhf}` chooses the restricted or unrestricted path.
`Auto` uses RHF for a closed shell and UHF otherwise. `Uhf` can be forced for
an even-electron singlet.

`ScfAccelerator::{None, Cdiis, AdiisCdiis}` selects plain iteration, Pulay
CDIIS, or the robust A-DIIS-to-CDIIS hybrid.

### `NddoOptions`

| Field | Default | Meaning |
|---|---:|---|
| `charge` | `0.0` | Molecular charge. |
| `multiplicity` | `1` | `2S + 1`. |
| `max_scf` | `200` | Maximum SCF iterations. |
| `e_tol` | `1e-8` | Energy convergence tolerance. |
| `p_tol` | `1e-7` | Density convergence tolerance. |
| `use_diis` | `true` | Legacy master switch; false forces no accelerator. |
| `accelerator` | `AdiisCdiis` | SCF extrapolation strategy. |
| `adiis_switch` | `0.1` | Commutator threshold for switching to CDIIS. |
| `reference` | `Auto` | RHF/UHF selection. |
| `level_shift_ev` | `0.0` | Virtual-space level shift. |
| `damping` | `0.0` | Old-density fraction in `[0, 1)`. |
| `d_penalty_start` | `0.0` | Initial transition-metal d-orbital Fock penalty, eV. |
| `d_penalty_iters` | `0` | Linear penalty annealing length. |
| `hessian_cutoff` | `None` | Optional CPHF pair cutoff radius, Bohr. |
| `hessian_memory_mb` | `1024` | Soft Hessian-workspace limit; zero disables. |
| `integral_memory_mb` | `0` | Pair-cache limit; zero uses the process default. |
| `scf_memory_mb` | `512` | SCF accelerator-history limit; zero disables. |

When `integral_memory_mb == 0`, `XNDO_MAX_PAIR_CACHE_MB` can override the
process default; otherwise the default pair-cache ceiling is 4 GiB.

### `NddoResult`

| Field | Meaning |
|---|---|
| `density` | Total AO density matrix. |
| `spin_density` | UHF `P_alpha - P_beta`, otherwise `None`. |
| `mo_energies`, `mo_coeff` | MO eigenvalues in eV and coefficient matrix. |
| `n_occ` | Occupied-orbital count; alpha count on UHF. |
| `electronic_ev`, `core_ev`, `total_ev` | Energy decomposition. |
| `heat_of_formation_kcal` | Heat of formation, kcal/mol. |
| `charges` | Net atomic NDO charges. |
| `dipole_debye`, `dipole_magnitude` | Molecular dipole vector and magnitude, Debye. |
| `homo_ev`, `lumo_ev` | Frontier levels when defined. |
| `iterations`, `converged`, `unrestricted` | SCF status. |

## MNDO-family parameters and derivatives

### `NddoParameters`, `NddoElement`, and `PairParams`

`NddoParameters` contains public `method`, `elements`, and symmetric `pair`
maps. Constructors `mndo()` and `mndod()` parse the bundled OpenMOPAC 23.2.5
tables. `for_method(method)` returns an owned table and
`cached_for_method(method)` returns a process-wide shared table. Both accept
only MNDO or MNDO/d. `element(z)` and `pair(zi, zj)` return checked parameter
references.

`PairParams` exposes pairwise core-core `alpha` and `x` values.

`NddoElement` exposes the complete runtime parameter record. The field groups
are:

- identity and basis: `z`, `n`, `n_s`, `n_p`, `n_d`, `n_orb`,
  `core_charge`, `main_group`;
- one-electron terms and Slater exponents: `u_ss`, `u_pp`, `u_dd`, `zeta_s`,
  `zeta_p`, `zeta_d`, `beta_s`, `beta_p`, `beta_d`, `zsn`, `zpn`, `zdn`;
- one-center interactions: `g_ss`, `g_sp`, `g_pp`, `g_p2`, `h_sp`, `f0sd`,
  `g2sd`, and optional `onecenter` spd integrals;
- core-core data: `alpha`, `poc`, and Gaussian triples in `gauss`;
- atomic bookkeeping: `occ_s`, `occ_p`, `occ_d`, `ndelec`, `eheat_ev`,
  `e_isol`, and `mass`;
- derived multipole data in Bohr: `dd`, `qq`, `rho0`, `rho1`, `rho2`, `ddp`,
  and `po`.

`has_p()` and `has_d()` test the AO layout.

### Ground-state gradient APIs

```rust
pub fn analytic_gradient(
    molecule: &Molecule,
    params: &NddoParameters,
    options: &NddoOptions,
    step: f64,
) -> Result<GradientResult>

pub fn closed_form_gradient(
    molecule: &Molecule,
    params: &NddoParameters,
    options: &NddoOptions,
) -> Result<GradientResult>

pub fn numerical_gradient(
    molecule: &Molecule,
    params: &NddoParameters,
    options: &NddoOptions,
    step: f64,
) -> Result<GradientResult>
```

`closed_form_gradient` is the production dual-number Hellmann-Feynman path.
`analytic_gradient` keeps the core term closed form but differences the
fixed-density electronic energy by `step` in Bohr. `numerical_gradient`
re-runs SCF at both displacements and is the independent finite-difference
reference.

`GradientResult` contains `scf`, `energy_ev`, `gradient`, `forces`, and
`max_gradient`; vector quantities use eV/Bohr.

### Hessian and vibrational APIs

```rust
pub fn analytic_hessian(
    molecule: &Molecule,
    params: &NddoParameters,
    options: &NddoOptions,
    step: f64,
) -> Result<Matrix>

pub fn numerical_hessian(
    molecule: &Molecule,
    params: &NddoParameters,
    options: &NddoOptions,
    step: f64,
) -> Result<Matrix>

pub fn vibrational_analysis(
    molecule: &Molecule,
    params: &NddoParameters,
    options: &NddoOptions,
    step: f64,
) -> Result<VibrationalModes>

pub fn vibrational_analysis_from_hessian(
    molecule: &Molecule,
    hessian: Matrix,
) -> Result<VibrationalModes>
```

`analytic_hessian` is the production CPHF/analytic-integral path, including
RHF/UHF and d AOs. `step` is retained for the guarded numerical fallback.
`numerical_hessian` differences analytic gradients. Both return eV/Bohr^2.

`vibrational_analysis` computes the analytic Hessian before mass weighting;
`vibrational_analysis_from_hessian` accepts a precomputed `(3N, 3N)` Hessian.
`VibrationalModes` contains the original `hessian`, ascending
`frequencies_cm`, and mass-weighted `eigenvalues` in
eV/(Angstrom^2 amu). Negative frequency values denote imaginary modes. The
conversion constant is available as `hessian::SQRT_EV_PER_ANG2_AMU_TO_CM`.

### Geometry optimization

```rust
pub fn optimize(
    molecule: &Molecule,
    method: Method,
    scf_options: &NddoOptions,
    opt: &OptOptions,
) -> Result<OptResult>
```

This is the single unified L-BFGS optimizer for CNDO/2, INDO, MNDO, MNDO/d,
MINDO/3, and ZINDO/S. Trial geometries use the unified single-point dispatcher;
accepted geometries use the matching analytic gradient dispatcher. Registered
methods without an engine return `XndoError::InvalidInput` and are never
substituted.

`OptOptions` fields are `max_iter` (`200`), `gtol` (`1e-3` eV/Bohr), and
`history` (`8`). Each must be positive and `gtol` must be finite. `OptResult`
contains `method`, final `molecule`, `energy_ev`, optional
`heat_of_formation_kcal`, `converged`, `iterations`, and a `trajectory` of
`OptStep`. Each `OptStep` records `energy_ev`, optional
`heat_of_formation_kcal`, `max_gradient`, and Bohr `positions`. A failed
Armijo line search returns the last accepted geometry with `converged == false`;
input, SCF, and linear-algebra failures return `Err`.

## CNDO/2 and INDO engine

```rust
pub fn run_cndo_indo(
    molecule: &Molecule,
    method: Method,
    options: &CndoIndoOptions,
) -> Result<CndoIndoResult>

pub fn cndo_indo_element(z: u8) -> Result<CndoIndoElement>
```

`run_cndo_indo` accepts only `Method::Cndo2` or `Method::Indo`. The bundled
CNDO/2 table supports H/Li/C/N/O/S; INDO supports H/Li/C/N/O. RHF and the
spin-resolved UHF extension are available.

`CndoIndoOptions` fields are `charge`, `multiplicity`, `reference`, `max_scf`,
`e_tol_ev`, `p_tol`, and `damping`; defaults are `0`, `1`, `Auto`, `300`,
`1e-8`, `1e-7`, and `0.20`.

`CndoIndoResult` fields are `method`, `density`, `fock`, optional `fock_beta`,
`mo_coeff`, optional `mo_coeff_beta`, `mo_energies_ev`, optional
`mo_energies_beta_ev`, `n_occ`, `n_alpha`, `n_beta`, optional `spin_density`,
`unrestricted`, `electronic_ev`, `core_ev`, `total_ev`, `charges`,
`iterations`, and `converged`.

`CndoIndoElement` exposes `z`, `symbol`, `core_charge`, `valence_electrons`,
`valence_shell`, `n_orb`, `has_d`, `bonding_parameter_ev`, `imu_s_ev`,
`imu_p_ev`, `imu_d_ev`, `zeta_s`, `zeta_p`, `zeta_d`, `indo_g1_ev`,
`indo_f2_ev`, and the six one-center coefficients
`indo_f0_coeff_{s,p}`, `indo_g1_coeff_{s,p}`, and
`indo_f2_coeff_{s,p}`.

## MINDO/3 engine

```rust
pub fn run_mindo3(
    molecule: &Molecule,
    options: &Mindo3Options,
) -> Result<Mindo3Result>

pub fn mindo3_element(z: u8) -> Result<Mindo3Element>
pub fn mindo3_pair(za: u8, zb: u8) -> Result<(f64, f64)>
```

`Mindo3Options` fields are `charge`, `multiplicity`, `reference`, `max_scf`,
`e_tol_ev`, `p_tol`, and `damping`; defaults are `0`, `1`, `Auto`, `250`,
`1e-8`, `1e-7`, and `0.20`.

`Mindo3Result` fields are `density`, `fock`, optional `fock_beta`, `mo_coeff`,
optional `mo_coeff_beta`, `mo_energies_ev`, optional `mo_energies_beta_ev`,
`n_occ`, `n_alpha`, `n_beta`, optional `spin_density`, `unrestricted`,
`electronic_ev`, `core_ev`, `total_ev`, `heat_of_formation_kcal`, `charges`,
`iterations`, and `converged`.

`Mindo3Element` exposes `z`, `n_s`, `n_p`, `n_orb`, `core_charge`, `uss`,
`upp`, `vs`, `vp`, `zeta_s`, `zeta_p`, `gss`, `gsp`, `gpp`, `gp2`, `hsp`,
`hp2`, `f03`, `e_isol_ev`, and `e_heat_kcal`. `mindo3_pair` returns the
method's pair coefficients for a supported atomic-number pair.

## ZINDO/S ground and excited states

### Parameters, options, and ground result

`ZindoParameters::standard()` parses a fresh embedded OpenMOPAC table;
`ZindoParameters::cached()` returns the process-wide table. `element(z)`
returns a checked `ZindoElement`, and `supported_sp_elements()` returns sorted
atomic numbers enabled by the one-/four-AO implementation.

`ZindoElement` fields are `z`, `n_orb`, `core_charge`, `zeta_sp`, two-component
`zeta_d`, `zeta_weight`, three-component `beta`, the 1-based `fg` array,
`n_s`, `n_p`, `n_d`, and `mass`.

`ZindoOptions` fields are `charge`, `multiplicity`, `reference`, `max_scf`,
`e_tol_ev`, `p_tol`, `damping`, `n_states`, `active_occupied`, and
`active_virtual`. Defaults are `0`, `1`, `Auto`, `250`, `1e-8`, `1e-7`,
`0.20`, `10`, `None`, and `None`.

```rust
pub fn run_zindo_s(
    molecule: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
) -> Result<ZindoResult>
```

`ZindoResult` fields are `density`, `fock`, optional `fock_beta`, `mo_coeff`,
optional `mo_coeff_beta`, `mo_energies_ev`, optional `mo_energies_beta_ev`,
`n_occ`, `n_alpha`, `n_beta`, optional `spin_density`, `unrestricted`,
`electronic_ev`, `core_ev`, `total_ev`, `charges`, `iterations`, and
`converged`.

### Spin labels and state data

`CisSpin::{Singlet, Triplet, Unrestricted}` identifies a CIS sector.
`as_str()` returns its lowercase name. `parse()` also accepts `s`/`1`, `t`/`3`,
and `ucis`/`u` aliases.

`CiContribution` stores one-based `occupied` and `virtual_orbital` indices plus
the CI `coefficient`.

`ExcitedState` exposes:

- `spin`, `spin_multiplicity`, and `s2_expectation`;
- vertical `energy_ev`, `energy_cm1`, and optional `wavelength_nm`;
- `state_total_energy_ev` and `state_total_energy_hartree`;
- `oscillator_strength`, transition dipoles in `transition_dipole_au` and
  `transition_dipole_debye`, `transition_dipole_magnitude_au`, and
  `dipole_strength_au2`;
- permanent and difference dipole vectors in atomic units and Debye, with
  corresponding Debye magnitudes;
- state `charges`, per-atom `hole_population` and `electron_population`, their
  Angstrom centroids, `charge_transfer_distance_angstrom`, and `dominant`
  configurations.

`ZindoSpectrum` contains `ground: ZindoResult`, `spin`, ground dipoles in
atomic units and Debye, and `states: Vec<ExcitedState>`.

### Spectrum functions

All spectrum functions have the same arguments:

```rust
fn(
    molecule: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
) -> Result<ZindoSpectrum>
```

- `zindo_s_cis` is the backward-compatible singlet RHF-CIS entry point.
- `uv_vis_spectrum` is the explicitly named singlet stick-spectrum entry
  point; it applies no broadening or solvent shift.
- `zindo_s_cis_spin(..., spin: CisSpin)` additionally selects singlet or
  triplet RHF-CIS. `CisSpin::Unrestricted` must use the UCIS API.
- `zindo_s_ucis` requires a UHF `ZindoResult` and performs spin-conserving
  spin-orbital CIS.

The historical constants re-exported at crate root are
`DEBYE_PER_E_BOHR` and `EV_TO_WAVENUMBER_CM1`. Other ZINDO-specific constants
are available in `xndo_rs::zindo`.

### State-gradient functions

```rust
pub fn zindo_s_cis_gradients(
    molecule: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
    spin: CisSpin,
) -> Result<ZindoCisGradientResult>

pub fn zindo_s_ucis_gradients(
    molecule: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
) -> Result<ZindoCisGradientResult>
```

`ExcitedStateGradient` contains one-based `root`, `spin`,
`excitation_energy_ev`, `state_total_energy_ev`,
`excitation_gradient`, and `state_gradient`. The vectors are eV/Bohr;
`state_gradient` differentiates ground energy plus the vertical excitation.

`ZindoCisGradientResult` contains the converged `ground`, its
`ground_gradient`, the `spin` sector, and `states`.

### State-Hessian functions

```rust
pub fn zindo_s_cis_hessians(
    molecule: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
    spin: CisSpin,
) -> Result<ZindoCisHessianResult>

pub fn zindo_s_ucis_hessians(
    molecule: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
) -> Result<ZindoCisHessianResult>
```

`ExcitedStateHessian` adds `excitation_hessian` and `state_hessian` matrices
in eV/Bohr^2 to the root, spin, energy, and gradient fields.
`ZindoCisHessianResult` contains `ground`, `ground_gradient`,
`ground_hessian`, `spin`, and `states`.

State-specific CIS/UCIS derivatives require each requested root to be isolated
from all other roots by at least 0.0001 eV. Near-degenerate spectra remain
available, but derivative APIs return a linear-algebra error rather than a
basis-dependent state derivative. Open-shell diatomic state derivatives are
also rejected as non-unique without a state-averaged or diabatic definition.

## Error type

`XndoError` variants are:

| Variant | Meaning |
|---|---|
| `Io(std::io::Error)` | File I/O failure. |
| `Parse { line, message }` | Text parsing failure with line number. |
| `InvalidInput(String)` | Invalid geometry, method, spin/reference, or option. |
| `MissingElement(u8)` | Selected method lacks an element block. |
| `UnsupportedDOrbitals(u8)` | An intentionally s/p-only low-level path received a d shell. |
| `MissingParameter(String)` | Required named or pair parameter is absent. |
| `LinearAlgebra(String)` | Eigensolver/response failure or non-unique state derivative. |
| `ResourceLimit { operation, required_mb, limit_mb }` | A guarded allocation exceeds its soft limit. |
| `ScfNotConverged { iterations, error }` | SCF exhausted `max_scf`. |

The type implements `Display`, `Error`, and conversion from `std::io::Error`.

## Public module map

The crate-root exports above are the supported application-facing surface.
The crate also keeps numerical building blocks public for validation and
advanced use:

| Module | Public role |
|---|---|
| `basis` | MNDO AO metadata and basis construction. |
| `cndo_indo` | CNDO/2 and INDO SCF, analytic derivatives, and donor parameters. |
| `constants` | Energy, length, dipole, and thermochemical conversions. |
| `data_tables` | Embedded CSV text, parsed tables, and element data. |
| `dual`, `dual2` | First- and second-order forward-mode scalar arithmetic. |
| `error` | `Result` and `XndoError`. |
| `fock` | NDDO RHF/UHF Fock construction. |
| `frame` | Generic-frame rotation and back-transformation helpers. |
| `gradient` | MNDO fixed-density, analytic, and numerical gradient internals. |
| `hamiltonian` | NDDO core-Hamiltonian and bounded pair-cache construction. |
| `hessian` | CPHF response, Hessians, and vibrational analysis. |
| `integrals`, `integrals_d` | s/p and explicit-d NDDO two-electron integrals. |
| `linalg` | Dense matrix type, eigensolver, and linear solves. |
| `math` | `Vec3` and `Mat3`. |
| `method` | Method registry, parsing, and capability discovery. |
| `mindo3` | MINDO/3 SCF, parameters, and analytic derivatives. |
| `onecenter` | One-center spd integral construction. |
| `optimizer` | Unified native-method L-BFGS optimizer and trajectory types. |
| `overlap`, `overlap_numeric` | Analytic/AD and numerical Slater overlap kernels. |
| `params` | MNDO/MNDO-d parameter parsing and derived multipoles. |
| `repulsion` | NDDO core-core energies and gradients. |
| `rotations` | Two-center local-frame rotation matrices. |
| `scf` | NDDO RHF/UHF driver and memory controls. |
| `system` | Molecules, atoms, XYZ parsing, and element symbols. |
| `xndo` | Unified native-method dispatch. |
| `zindo` | ZINDO/S SCF, CIS/UCIS properties, and state derivatives. |
| `python` | PyO3 extension registration when feature `python` is enabled. |

These modules intentionally expose specialized helpers that are easier to use
through the crate-root dispatchers. Run `cargo doc --all-features --no-deps`
for generated signatures and source-linked details of every low-level item.
