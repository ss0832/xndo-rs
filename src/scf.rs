// SPDX-License-Identifier: GPL-3.0-or-later

//! NDDO SCF driver: restricted (RHF) closed shell and unrestricted (UHF) open shell.
//!
//! Because NDDO assumes an orthonormal AO basis, the working equations are the plain
//! eigenproblem `F C = C epsilon` (no `S`); overlap enters only the resonance term of `H_core`.
//! The initial density is a **superposition of atomic densities** (`sad_density`): the
//! exact free-atom density in a minimal valence basis, far better than the bare-core guess
//! and charge convergence is accelerated with the A-DIIS/CDIIS hybrid on the `[F,P]` commutator.
//!
//! PROVENANCE: derived from MOPAC (Molecular Orbital PACkage) v23.2.5,
//! Copyright 2021 Virginia Polytechnic Institute and State University,
//! licensed under the Apache License, Version 2.0.
//! UPSTREAM: src/SCF and src/models/parameters_C.F90 (shell occupancies).
//! MODIFIED for xndo-rs v0.3.0 on 2026-09-14:
//! A-DIIS/CDIIS acceleration and the superposition-of-atomic-densities guess
//! are this project's; the shell occupancies, the helium exception and the
//! heat-of-formation convention follow upstream.
//! Retained notices: NOTICE; per-file record: THIRD_PARTY_NOTICES.md.

use crate::basis::Basis;
use crate::error::{Result, XndoError};
use crate::fock::{build_fock, build_fock_spin};
use crate::hamiltonian::{build_core_limited, CoreHamiltonian};
use crate::linalg::{symmetric_eigen, Matrix};
use crate::math::Vec3;
use crate::orbitals::OrbitalEnergies;
use crate::params::NddoParameters;
use crate::repulsion::core_core_energy;
pub use crate::scf_accel::ScfAccelerator;
use crate::scf_accel::{lowest_solution, open_shell_starts, uniform_valence_density, Accelerator};
use crate::system::Molecule;

/// Choice of SCF reference (restricted vs unrestricted), independent of the
/// spin multiplicity. `Auto` picks RHF for a closed shell and UHF for an open
/// shell; `Rhf`/`Uhf` force the reference even for a singlet (a UHF singlet can
/// break spin symmetry, e.g. stretched bonds or diradicals).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Reference {
    #[default]
    Auto,
    Rhf,
    Uhf,
}

#[derive(Clone, Debug)]
pub struct NddoOptions {
    pub charge: f64,
    pub multiplicity: usize,
    pub max_scf: usize,
    pub e_tol: f64,
    pub p_tol: f64,
    /// Legacy flag: `false` forces [`ScfAccelerator::None`] regardless of `accelerator`.
    pub use_diis: bool,
    pub accelerator: ScfAccelerator,
    /// Commutator-error norm below which the ADIIS/CDIIS hybrid switches to CDIIS.
    pub adiis_switch: f64,
    /// Restricted/unrestricted reference. `Auto` (default) uses RHF for closed
    /// shells and UHF for open shells; `Uhf` forces the unrestricted path even
    /// for a singlet.
    pub reference: Reference,
    /// Level shift (eV) applied to the virtual space during diagonalization
    /// (SaundersHillier). It changes the SCF path/basin without changing the
    /// converged fixed point, preventing the variational collapse to
    /// unphysical charge-transfer states seen for open-shell transition metals.
    pub level_shift_ev: f64,
    /// Density damping factor in `[0,1)`: the density fed to the next Fock is
    /// `(1 - damping) * P_new + damping * P_old`. Slows early-iteration charge
    /// sloshing, which together with the level shift keeps the SCF in the physical basin.
    pub damping: f64,
    /// Transition-metal d-block penalty (eV) added to the metal d-orbital Fock
    /// diagonal at iteration 0 and linearly annealed to zero over
    /// `d_penalty_iters`. This raises the metal d orbitals early so electrons
    /// stay on the ligands (the physical charge distribution), then relaxes to
    /// the true Fock, steering open-shell transition-metal SCF away from the
    /// unphysical metal-anion collapse that a plain iteration / level shift
    /// cannot escape (the collapse fills the *low* metal-d orbitals, which a
    /// virtual-space shift leaves untouched).
    pub d_penalty_start: f64,
    /// Number of iterations over which `d_penalty_start` anneals to zero.
    pub d_penalty_iters: usize,
    /// Optional CPHF/Hessian interaction cutoff radius `r_off` (Bohr). When set,
    /// the coupled-perturbed response damps atom-pair contributions smoothly to
    /// zero between `r_off / 2` and `r_off` (a C2 quintic switch), so distant
    /// pairs whose response is negligible are skipped. `None` (the default)
    /// computes the exact response and is **bit-identical** to the uncut result.
    pub hessian_cutoff: Option<f64>,
    /// Soft memory budget (MiB) for Hessian workspaces. It bounds concurrent
    /// finite-difference jobs and CPHF response batches, and rejects an analytic
    /// Hessian before allocation if its resident occupied-virtual blocks alone
    /// exceed the limit. Set to `0` to disable the guard.
    pub hessian_memory_mb: usize,
    /// Soft memory budget (MiB) for the resident `O(N2)` two-electron pair
    /// cache built by [`crate::hamiltonian::build_core_limited`]. Exceeding it
    /// raises [`XndoError::ResourceLimit`] before the allocation is attempted,
    /// instead of letting the allocator abort the process. `0` disables the
    /// check. Overridden per-process by `XNDO_MAX_PAIR_CACHE_MB`.
    pub integral_memory_mb: usize,
    /// Soft memory budget (MiB) for the SCF convergence-accelerator history.
    /// Each retained iterate costs `O(nao2)`; the effective DIIS depth is
    /// reduced (never below 2) so the history fits. `0` disables the cap and
    /// always keeps the full depth. This is the difference between a
    /// 4800-basis-function job needing 4.4 GiB of history and 0.5 GiB.
    pub scf_memory_mb: usize,
}

impl Default for NddoOptions {
    fn default() -> Self {
        Self {
            charge: 0.0,
            multiplicity: 1,
            max_scf: 200,
            e_tol: 1.0e-8,
            p_tol: 1.0e-7,
            use_diis: true,
            accelerator: ScfAccelerator::AdiisCdiis,
            adiis_switch: 0.1,
            reference: Reference::Auto,
            level_shift_ev: 0.0,
            damping: 0.0,
            d_penalty_start: 0.0,
            d_penalty_iters: 0,
            hessian_cutoff: None,
            hessian_memory_mb: 1024,
            integral_memory_mb: 0,
            scf_memory_mb: 512,
        }
    }
}

/// Resolve the two-electron pair-cache ceiling: the explicit option when set,
/// otherwise the process default (`XNDO_MAX_PAIR_CACHE_MB`, else 4 GiB).
pub(crate) fn pair_cache_limit(options: &NddoOptions) -> usize {
    if options.integral_memory_mb > 0 {
        options.integral_memory_mb
    } else {
        crate::hamiltonian::default_pair_cache_limit_mb()
    }
}

/// Apply a SaundersHillier level shift to a Fock matrix before diagonalization:
/// `F' = F + shift(I  P)`, raising the virtual space by ~`shift`. `P` is the
/// (idempotent, orthonormal-basis) density projector for the spin channel
/// (occupation 1 per orbital); for RHF pass `12P_tot`.
fn level_shift(f: &Matrix, p_proj: &Matrix, shift: f64) -> Matrix {
    if shift == 0.0 {
        return f.clone();
    }
    let mut out = f.clone();
    let n = f.rows;
    for i in 0..n {
        for j in 0..n {
            let ident = if i == j { 1.0 } else { 0.0 };
            out[(i, j)] += shift * (ident - p_proj[(i, j)]);
        }
    }
    out
}

#[derive(Clone, Debug)]
pub struct NddoResult {
    pub density: Matrix,
    /// Spin density `P_  P_` (open-shell UHF only; `None` for RHF).
    pub spin_density: Option<Matrix>,
    pub mo_energies: Vec<f64>,
    pub mo_coeff: Matrix,
    /// Beta-channel eigenvalues (UHF only; `None` for RHF).
    pub mo_energies_beta: Option<Vec<f64>>,
    /// Beta-channel coefficients (UHF only; `None` for RHF).
    pub mo_coeff_beta: Option<Matrix>,
    /// Occupied alpha count; equals `n_alpha`, kept under the older name.
    pub n_occ: usize,
    pub n_alpha: usize,
    pub n_beta: usize,
    pub electronic_ev: f64,
    pub core_ev: f64,
    pub total_ev: f64,
    pub heat_of_formation_kcal: f64,
    pub charges: Vec<f64>,
    pub dipole_debye: Vec3,
    pub dipole_magnitude: f64,
    /// Highest occupied spin orbital over *both* channels; see
    /// [`crate::orbitals::OrbitalEnergies::homo_ev`].
    pub homo_ev: Option<f64>,
    /// Lowest unoccupied spin orbital over *both* channels.
    pub lumo_ev: Option<f64>,
    pub iterations: usize,
    pub converged: bool,
    /// True when the UHF (open-shell) path was used.
    pub unrestricted: bool,
}

impl NddoResult {
    /// The orbital energies and the frontier quantities read off them.
    pub fn orbitals(&self) -> OrbitalEnergies {
        OrbitalEnergies::new(
            self.mo_energies.clone(),
            self.n_alpha,
            self.mo_energies_beta.clone(),
            self.n_beta,
        )
    }

    /// HOMO-LUMO gap in eV; `None` unless both frontier orbitals exist.
    pub fn homo_lumo_gap_ev(&self) -> Option<f64> {
        Some(self.lumo_ev? - self.homo_ev?)
    }
}

pub struct NddoCalculator {
    pub params: NddoParameters,
    pub options: NddoOptions,
}

impl NddoCalculator {
    pub fn new(params: NddoParameters) -> Self {
        Self {
            params,
            options: NddoOptions::default(),
        }
    }
    pub fn with_options(params: NddoParameters, options: NddoOptions) -> Self {
        Self { params, options }
    }
    pub fn calculate(&self, molecule: &Molecule) -> Result<NddoResult> {
        run_nddo_with_parameters(molecule, &self.params, &self.options)
    }
}

struct ScfState {
    density: Matrix,
    spin_density: Option<Matrix>,
    mo_energies: Vec<f64>,
    mo_coeff: Matrix,
    mo_energies_beta: Option<Vec<f64>>,
    mo_coeff_beta: Option<Matrix>,
    n_occ: usize,
    n_beta: usize,
    electronic_ev: f64,
    converged: bool,
    iterations: usize,
    unrestricted: bool,
}

/// Enable faer's multi-threaded (rayon) parallelism once per process. faer 0.24
/// defaults to sequential, leaving the O(N3) matmul/eigendecomposition in the
/// SCF single-threaded; on a many-core machine this is the single biggest
/// large-system speedup.
fn enable_faer_parallelism() {
    use std::sync::Once;
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        faer::set_global_parallelism(faer::Par::rayon(0));
    });
}

pub fn run_nddo_with_parameters(
    molecule: &Molecule,
    params: &NddoParameters,
    options: &NddoOptions,
) -> Result<NddoResult> {
    enable_faer_parallelism();
    let charge =
        reconcile_scalar_metadata("charge", molecule.charge, options.charge, 0.0, 1.0e-12)?;
    let multiplicity = reconcile_usize_metadata(
        "multiplicity",
        molecule.multiplicity,
        options.multiplicity,
        1,
    )?;
    if multiplicity < 1 {
        return Err(XndoError::InvalidInput(
            "multiplicity must be >= 1".to_string(),
        ));
    }
    let timing = std::env::var("XNDO_TIMING").is_ok();
    let t_all = std::time::Instant::now();
    let basis = Basis::build(molecule, params)?;
    let core = build_core_limited(molecule, &basis, params, pair_cache_limit(options))?;
    if timing {
        eprintln!(
            "[timing] basis+core (integrals): {:.3}s  (nao={})",
            t_all.elapsed().as_secs_f64(),
            basis.nao
        );
    }

    let mut n_elec = 0.0;
    for atom in &molecule.atoms {
        n_elec += params.element(atom.z)?.core_charge;
    }
    n_elec -= charge;
    let n_elec_int = n_elec.round() as i64;
    if (n_elec - n_elec_int as f64).abs() > 1.0e-6 || n_elec_int < 0 {
        return Err(XndoError::InvalidInput(format!(
            "invalid electron count {n_elec}"
        )));
    }
    let n_unpaired = (multiplicity - 1) as i64;
    if (n_elec_int - n_unpaired) < 0 || (n_elec_int - n_unpaired) % 2 != 0 {
        return Err(XndoError::InvalidInput(format!(
            "electron count {n_elec_int} is incompatible with multiplicity {} (need same parity)",
            multiplicity
        )));
    }
    let n_alpha = ((n_elec_int + n_unpaired) / 2) as usize;
    let n_beta = ((n_elec_int - n_unpaired) / 2) as usize;

    // Restricted vs unrestricted: honour an explicit request, else pick by shell.
    let use_uhf = match options.reference {
        Reference::Auto => n_alpha != n_beta,
        Reference::Uhf => true,
        Reference::Rhf => {
            if n_alpha != n_beta {
                return Err(XndoError::InvalidInput(format!(
                    "RHF requested for an open-shell system (n_alpha={n_alpha} != n_beta={n_beta}); use UHF or Auto"
                )));
            }
            false
        }
    };
    // Open-shell transition-metal systems can collapse to an unphysical
    // charge-transfer minimum without a level shift. When the caller has not
    // set one, default to a moderate shift for the UHF path; it guides the SCF
    // to the aufbau state and does not change the converged energy of a
    // well-behaved case (SaundersHillier).
    let has_open_d = use_uhf
        && molecule
            .atoms
            .iter()
            .any(|a| params.element(a.z).map(|e| e.has_d()).unwrap_or(false));
    let run_scf = |shift: f64, damping: f64, d_pen: f64, d_pen_iters: usize| -> Result<ScfState> {
        let mut opts = options.clone();
        opts.level_shift_ev = shift;
        opts.damping = damping;
        opts.d_penalty_start = d_pen;
        opts.d_penalty_iters = d_pen_iters;
        if use_uhf {
            // One SCF per candidate start, keeping the lowest: an open shell's
            // unpaired electron is assigned by the *first* Fock matrix, which a
            // guess can order wrongly in a way converging harder never undoes.
            // See `scf_accel::open_shell_starts`.
            let sad = sad_density(molecule, &basis, params)?;
            let uniform =
                uniform_valence_density(basis.nao, &valence_shells(molecule, &basis, params)?);
            lowest_solution(
                open_shell_starts(&sad, &uniform, n_alpha, n_beta),
                |start| {
                    uhf_loop(
                        molecule, &basis, params, &core, n_alpha, n_beta, &opts, start,
                    )
                },
                |state| state.electronic_ev,
            )
        } else {
            rhf_loop(molecule, &basis, params, &core, n_alpha, &opts)
        }
    };
    // Open-shell transition-metal SCF collapses to an unphysical metal-anion
    // fixed point from any guess (a genuine lower/adjacent SCF solution; the
    // Hamiltonian itself is bit-exact vs MOPAC). A metal-d Fock penalty, held
    // high early and annealed to zero, keeps electrons on the ligands and
    // tracks the physical basin down to the true (penalty-free) fixed point.
    // Level shift + damping alone cannot escape it (the collapse fills the low
    // metal-d orbitals, which a virtual-space shift leaves untouched).
    let state = if options.level_shift_ev != 0.0 || options.damping != 0.0 {
        run_scf(options.level_shift_ev, options.damping, 0.0, 0)?
    } else if has_open_d {
        // Anneal a 60 eV metal-d penalty over 60 iterations, with mild damping.
        // Recovers the physical basin for several open-shell TM halides; for the
        // remainder it is at worst neutral (the penalty vanishes before
        // convergence, so any converged result is a true fixed point).
        run_scf(0.0, 0.3, 60.0, 60)?
    } else {
        run_scf(0.0, 0.0, 0.0, 0)?
    };
    if timing {
        eprintln!(
            "[timing] SCF ({} iters): {:.3}s",
            state.iterations,
            t_all.elapsed().as_secs_f64()
        );
    }

    let core_ev = core_core_energy(molecule, params)?;
    let electronic_ev = state.electronic_ev;
    let total_ev = electronic_ev + core_ev;

    let mut e_isol_sum = 0.0;
    let mut eheat_sum = 0.0;
    for atom in &molecule.atoms {
        let e = params.element(atom.z)?;
        e_isol_sum += e.e_isol;
        eheat_sum += e.eheat_ev;
    }
    let heat_of_formation_kcal = (total_ev - e_isol_sum + eheat_sum) * params.constants.ev_to_kcal;

    if !state.converged {
        return Err(XndoError::ScfNotConverged {
            iterations: state.iterations,
            error: f64::NAN,
        });
    }

    // Mulliken net charges from the total density.
    let mut charges = vec![0.0; molecule.atoms.len()];
    for (ia, atom) in molecule.atoms.iter().enumerate() {
        let off = basis.atom_offset[ia];
        let n = basis.atom_norb[ia];
        let mut pop = 0.0;
        for mu in 0..n {
            pop += state.density[(off + mu, off + mu)];
        }
        charges[ia] = params.element(atom.z)?.core_charge - pop;
    }

    // Dipole: point-charge term + sp hybrid polarization (both in eBohr).
    let mut dip = Vec3::zero();
    for (ia, atom) in molecule.atoms.iter().enumerate() {
        dip += atom.position * charges[ia];
        let elem = params.element(atom.z)?;
        if elem.has_p() {
            let off = basis.atom_offset[ia];
            let hyb = -2.0 * elem.dd;
            dip += Vec3::new(
                hyb * state.density[(off, off + 1)],
                hyb * state.density[(off, off + 2)],
                hyb * state.density[(off, off + 3)],
            );
            if elem.has_d() {
                // One-center p-d hybridization term used by MNDO/d. The AO
                // order is s, px, py, pz, d(x2-y2), d(xz), d(z2), d(yz), d(xy).
                let inv_sqrt3 = 1.0 / 3.0_f64.sqrt();
                let dx = state.density[(off + 3, off + 5)]
                    + state.density[(off + 1, off + 4)]
                    + state.density[(off + 2, off + 8)]
                    - inv_sqrt3 * state.density[(off + 1, off + 6)];
                let dy = state.density[(off + 3, off + 7)] - state.density[(off + 2, off + 4)]
                    + state.density[(off + 1, off + 8)]
                    - inv_sqrt3 * state.density[(off + 2, off + 6)];
                let dz = state.density[(off + 1, off + 5)]
                    + state.density[(off + 2, off + 7)]
                    + 2.0 * inv_sqrt3 * state.density[(off + 3, off + 6)];
                dip += Vec3::new(dx, dy, dz) * (-2.0 * elem.ddp[5]);
            }
        }
    }
    let dipole_debye = dip * params.constants.debye_per_e_bohr();
    let dipole_magnitude = dipole_debye.norm();

    // Frontier pair over both spin channels. Reading the alpha channel alone
    // (what this returned before v0.3.0) gives the lowest unoccupied *alpha*
    // orbital, which for an open shell is above the beta LUMO.
    let orbitals = OrbitalEnergies::new(
        state.mo_energies.clone(),
        state.n_occ,
        state.mo_energies_beta.clone(),
        state.n_beta,
    );
    let homo_ev = orbitals.homo_ev();
    let lumo_ev = orbitals.lumo_ev();

    Ok(NddoResult {
        density: state.density,
        spin_density: state.spin_density,
        mo_energies: state.mo_energies,
        mo_coeff: state.mo_coeff,
        mo_energies_beta: state.mo_energies_beta,
        mo_coeff_beta: state.mo_coeff_beta,
        n_occ: state.n_occ,
        n_alpha: state.n_occ,
        n_beta: state.n_beta,
        electronic_ev,
        core_ev,
        total_ev,
        heat_of_formation_kcal,
        charges,
        dipole_debye,
        dipole_magnitude,
        homo_ev,
        lumo_ev,
        iterations: state.iterations,
        converged: state.converged,
        unrestricted: state.unrestricted,
    })
}

/// **Superposition of Atomic Densities (SAD)** initial guess: a block-diagonal density built
/// from each atom's spherically-averaged neutral valence configuration `s^{min(2,Zv)} p^{...}`.
/// In a minimal valence NDDO basis a free atom's density *is* diagonal (spherical symmetry), so
/// this superposition is the exact isolated-atom density and a far better SCF start than the
/// bare core (zero-density) guess  fewer iterations, more robust for larger systems.
fn sad_density(molecule: &Molecule, basis: &Basis, params: &NddoParameters) -> Result<Matrix> {
    let nao = basis.nao;
    let mut p = Matrix::zeros(nao, nao);
    for (ia, atom) in molecule.atoms.iter().enumerate() {
        let elem = params.element(atom.z)?;
        let off = basis.atom_offset[ia];
        let n = basis.atom_norb[ia];
        if n == 0 {
            continue;
        }
        // MOPAC initial shell occupancies (ios/iop/iod), spherically averaged.
        p[(off, off)] = elem.occ_s;
        if n >= 4 {
            let per_p = elem.occ_p / 3.0;
            for k in 1..4 {
                p[(off + k, off + k)] = per_p;
            }
        }
        if n >= 9 {
            let per_d = elem.occ_d / 5.0;
            for k in 4..9 {
                p[(off + k, off + k)] = per_d;
            }
        }
    }
    Ok(p)
}

/// `(first AO index, AO count, core charge)` per atom, the shape the shared
/// guess builders in [`crate::scf_accel`] take.
fn valence_shells(
    molecule: &Molecule,
    basis: &Basis,
    params: &NddoParameters,
) -> Result<Vec<(usize, usize, f64)>> {
    molecule
        .atoms
        .iter()
        .enumerate()
        .map(|(ia, atom)| -> Result<(usize, usize, f64)> {
            Ok((
                basis.atom_offset[ia],
                basis.atom_norb[ia],
                params.element(atom.z)?.core_charge,
            ))
        })
        .collect()
}

/// Build a density `P = w _{k<n_occ} c_k c_kT` from MO coefficients (`w` = 2 for RHF, 1 for UHF).
fn density_from_coeff(c: &Matrix, n_occ: usize, weight: f64) -> Matrix {
    let nao = c.rows;
    if nao == 0 {
        return Matrix::zeros(0, 0);
    }
    if n_occ == 0 {
        return Matrix::zeros(nao, nao);
    }
    // P = weight  C_occ C_occT, directly into the row-major result.
    c.leading_columns_gram(n_occ, weight)
}

fn reconcile_scalar_metadata(
    name: &str,
    molecule_value: f64,
    option_value: f64,
    default: f64,
    tolerance: f64,
) -> Result<f64> {
    let molecule_set = (molecule_value - default).abs() > tolerance;
    let option_set = (option_value - default).abs() > tolerance;
    if molecule_set && option_set && (molecule_value - option_value).abs() > tolerance {
        return Err(XndoError::InvalidInput(format!(
            "conflicting {name}: molecule={molecule_value}, options={option_value}"
        )));
    }
    Ok(if option_set {
        option_value
    } else {
        molecule_value
    })
}

fn reconcile_usize_metadata(
    name: &str,
    molecule_value: usize,
    option_value: usize,
    default: usize,
) -> Result<usize> {
    let molecule_set = molecule_value != default;
    let option_set = option_value != default;
    if molecule_set && option_set && molecule_value != option_value {
        return Err(XndoError::InvalidInput(format!(
            "conflicting {name}: molecule={molecule_value}, options={option_value}"
        )));
    }
    Ok(if option_set {
        option_value
    } else {
        molecule_value
    })
}

fn rms_diff(a: &Matrix, b: &Matrix) -> f64 {
    a.rms_difference(b)
}

fn rhf_loop(
    molecule: &Molecule,
    basis: &Basis,
    params: &NddoParameters,
    core: &CoreHamiltonian,
    n_occ: usize,
    options: &NddoOptions,
) -> Result<ScfState> {
    let nao = basis.nao;
    let mut density = sad_density(molecule, basis, params)?; // SAD initial guess
    let mut e_old = 0.0;
    let mut mo_energies = vec![0.0; nao];
    let mut mo_coeff = Matrix::zeros(nao, nao);
    let mut converged = false;
    let mut iterations = 0;
    // MOPAC's MNDO He parameterization has formal p AOs but an s-only core
    // attraction.  On charged He systems DIIS can converge to a different
    // charge-transfer fixed point; the unaccelerated Roothaan path follows the
    // MOPAC basin reproducibly.
    let contains_helium = molecule.atoms.iter().any(|atom| atom.z == 2);
    let accel = if options.use_diis && !contains_helium {
        options.accelerator
    } else {
        ScfAccelerator::None
    };
    let timing = std::env::var("XNDO_TIMING").is_ok();
    let (mut t_fock, mut t_eigen, mut t_accel, mut t_density) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut accel = Accelerator::new(nao, accel, options.adiis_switch, options.scf_memory_mb);

    for iter in 0..options.max_scf {
        iterations = iter + 1;
        let t0 = std::time::Instant::now();
        let f = build_fock(molecule, basis, params, core, &density)?;
        if timing {
            t_fock += t0.elapsed().as_secs_f64();
        }
        let e_elec = 0.5 * (density.frobenius_dot(&core.h_core) + density.frobenius_dot(&f));

        let t_a = std::time::Instant::now();
        let f_use = accel.step(f, &density);
        let t1 = std::time::Instant::now();
        if timing {
            t_accel += t_a.elapsed().as_secs_f64();
        }
        let (eps, c) = symmetric_eigen(&f_use)?;
        if timing {
            t_eigen += t1.elapsed().as_secs_f64();
        }
        let t2 = std::time::Instant::now();
        let p_new = density_from_coeff(&c, n_occ, 2.0);
        let dp = rms_diff(&p_new, &density);
        if timing {
            t_density += t2.elapsed().as_secs_f64();
        }
        let de = (e_elec - e_old).abs();

        mo_energies = eps;
        mo_coeff = c;
        density = p_new;
        e_old = e_elec;
        // The commutator above belongs to the *pre-extrapolation* density.  In
        // a DIIS step it can be small even when the extrapolated Fock produces
        // a materially different density, which falsely accepted an excited
        // UHF fixed point.  Density self-consistency is the authoritative
        // condition here; the next iteration evaluates its commutator.
        if iter > 0 && de < options.e_tol && dp < options.p_tol {
            converged = true;
            break;
        }
    }
    if timing {
        eprintln!(
            "[timing]   RHF fock: {t_fock:.3}s  accel: {t_accel:.3}s  eigen: {t_eigen:.3}s  density: {t_density:.3}s  ({iterations} iters, diis depth {})",
            accel.depth()
        );
    }
    let f_final = build_fock(molecule, basis, params, core, &density)?;
    let electronic_ev =
        0.5 * (density.frobenius_dot(&core.h_core) + density.frobenius_dot(&f_final));

    Ok(ScfState {
        density,
        spin_density: None,
        mo_energies,
        mo_coeff,
        mo_energies_beta: None,
        mo_coeff_beta: None,
        n_occ,
        n_beta: n_occ,
        electronic_ev,
        converged,
        iterations,
        unrestricted: false,
    })
}

#[allow(clippy::too_many_arguments)]
fn uhf_loop(
    molecule: &Molecule,
    basis: &Basis,
    params: &NddoParameters,
    core: &CoreHamiltonian,
    n_alpha: usize,
    n_beta: usize,
    options: &NddoOptions,
    start: (Matrix, Matrix),
) -> Result<ScfState> {
    let nao = basis.nao;
    // Transition-metal d AOs to penalize during annealing (open-shell d collapse fix).
    let d_pen_aos: Vec<usize> = if options.d_penalty_start > 0.0 {
        basis
            .aos
            .iter()
            .enumerate()
            .filter(|(_, ao)| {
                ao.orb >= 4 && params.element(ao.z).map(|e| !e.main_group).unwrap_or(false)
            })
            .map(|(i, _)| i)
            .collect()
    } else {
        Vec::new()
    };
    // The starting densities come from the caller, which runs this loop once per
    // candidate start and keeps the lowest solution. The spin split itself puts
    // the excess spin where each AO has room for it rather than in proportion to
    // how occupied it already is; both are `scf_accel`'s business.
    let (mut pa, mut pb) = start;
    let mut e_old = 0.0;
    let mut eps_a = vec![0.0; nao];
    let mut c_a = Matrix::zeros(nao, nao);
    let mut eps_b = vec![0.0; nao];
    let mut c_b = Matrix::zeros(nao, nao);
    let mut converged = false;
    let mut iterations = 0;

    // Same helium caveat as the RHF path: MOPAC's MNDO He parameterization has
    // formal p AOs but an s-only core attraction, and an extrapolated path can
    // land on a different charge-transfer fixed point.
    let contains_helium = molecule.atoms.iter().any(|atom| atom.z == 2);
    let accel_mode = if options.use_diis && !contains_helium {
        options.accelerator
    } else {
        ScfAccelerator::None
    };
    // One accelerator for both channels, extrapolating the stacked pair with a
    // single set of coefficients. A previous comment here claimed A-DIIS could
    // not be reused for UHF "without a spin-resolved energy functional", and so
    // left UHF unaccelerated unless the caller asked for plain CDIIS by name.
    // That is not so: the Frobenius product of the stacked matrices is
    // `Tr[P_a F_a] + Tr[P_b F_b]`, which *is* the UHF pairing, so the A-DIIS
    // functional is correct on the stack with no change.
    let mut accel =
        Accelerator::new_uhf(nao, accel_mode, options.adiis_switch, options.scf_memory_mb);

    for iter in 0..options.max_scf {
        iterations = iter + 1;
        let mut p_tot = pa.clone();
        for (t, b) in p_tot.as_mut_slice().iter_mut().zip(pb.as_slice()) {
            *t += *b;
        }
        let fa = build_fock_spin(molecule, basis, params, core, &p_tot, &pa)?;
        let fb = build_fock_spin(molecule, basis, params, core, &p_tot, &pb)?;

        let e_elec = 0.5
            * (p_tot.frobenius_dot(&core.h_core) + pa.frobenius_dot(&fa) + pb.frobenius_dot(&fb));

        let (fa_use, fb_use) = accel.step_uhf(fa, fb, &pa, &pb);
        // Level shift (path only; the converged fixed point is unchanged).
        let shift = options.level_shift_ev;
        let mut fa_diag = level_shift(&fa_use, &pa, shift);
        let mut fb_diag = level_shift(&fb_use, &pb, shift);
        // Transition-metal d-block penalty, linearly annealed to zero (path
        // only; vanishes before convergence so the fixed point is the true one).
        if !d_pen_aos.is_empty() {
            let frac = 1.0 - (iter as f64) / (options.d_penalty_iters.max(1) as f64);
            let pen = options.d_penalty_start * frac.max(0.0);
            if pen > 0.0 {
                for &k in &d_pen_aos {
                    fa_diag[(k, k)] += pen;
                    fb_diag[(k, k)] += pen;
                }
            }
        }
        let (ea_eps, ca) = symmetric_eigen(&fa_diag)?;
        let (eb_eps, cb) = symmetric_eigen(&fb_diag)?;
        let mut pa_new = density_from_coeff(&ca, n_alpha, 1.0);
        let mut pb_new = density_from_coeff(&cb, n_beta, 1.0);

        let dp = rms_diff(&pa_new, &pa) + rms_diff(&pb_new, &pb);
        // Density damping (path only): blend with the previous density.
        let lambda = options.damping;
        if lambda > 0.0 {
            for (n, o) in pa_new.as_mut_slice().iter_mut().zip(pa.as_slice()) {
                *n = (1.0 - lambda) * *n + lambda * *o;
            }
            for (n, o) in pb_new.as_mut_slice().iter_mut().zip(pb.as_slice()) {
                *n = (1.0 - lambda) * *n + lambda * *o;
            }
        }
        let de = (e_elec - e_old).abs();
        eps_a = ea_eps;
        c_a = ca;
        eps_b = eb_eps;
        c_b = cb;
        pa = pa_new;
        pb = pb_new;
        e_old = e_elec;
        // See the RHF loop: `err_norm` is evaluated before DIIS extrapolation,
        // so it must not by itself accept the newly generated density. Never
        // accept convergence while the d-block penalty is still active.
        let penalty_active = !d_pen_aos.is_empty() && iter < options.d_penalty_iters;
        if iter > 0 && !penalty_active && de < options.e_tol && dp < options.p_tol {
            converged = true;
            break;
        }
    }

    let mut density = pa.clone();
    for (t, b) in density.as_mut_slice().iter_mut().zip(pb.as_slice()) {
        *t += *b;
    }
    let mut spin = pa.clone();
    for (s, b) in spin.as_mut_slice().iter_mut().zip(pb.as_slice()) {
        *s -= *b;
    }
    // Final energy.
    let fa = build_fock_spin(molecule, basis, params, core, &density, &pa)?;
    let fb = build_fock_spin(molecule, basis, params, core, &density, &pb)?;
    let electronic_ev =
        0.5 * (density.frobenius_dot(&core.h_core) + pa.frobenius_dot(&fa) + pb.frobenius_dot(&fb));

    Ok(ScfState {
        density,
        spin_density: Some(spin),
        mo_energies: eps_a,
        mo_coeff: c_a,
        // Kept rather than dropped: the beta eigenvalues are what make the
        // reported LUMO the lowest unoccupied orbital rather than the lowest
        // unoccupied *alpha* orbital, which for a doublet is a different one.
        mo_energies_beta: Some(eps_b),
        mo_coeff_beta: Some(c_b),
        n_occ: n_alpha,
        n_beta,
        electronic_ev,
        converged,
        iterations,
        unrestricted: true,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(xyz: &str, charge: f64) -> NddoResult {
        let mol = Molecule::from_xyz_str(xyz, charge).unwrap();
        let params = NddoParameters::mndo().unwrap();
        run_nddo_with_parameters(&mol, &params, &NddoOptions::default()).unwrap()
    }

    fn run_mult(xyz: &str, charge: f64, mult: usize) -> NddoResult {
        let mol = Molecule::from_xyz_str(xyz, charge).unwrap();
        let params = NddoParameters::mndo().unwrap();
        let opts = NddoOptions {
            charge,
            multiplicity: mult,
            ..NddoOptions::default()
        };
        run_nddo_with_parameters(&mol, &params, &opts).unwrap()
    }

    #[test]
    fn scf_memory_budget_does_not_change_the_converged_energy() {
        // Trimming the accelerator history is a path change only: the SCF fixed
        // point (and hence every reported quantity) must be identical.
        let xyz =
            "4\nformaldehyde\nC 0.0 0.0 0.0\nO 0.0 0.0 1.21\nH 0.94 0.0 -0.54\nH -0.94 0.0 -0.54\n";
        let molecule = Molecule::from_xyz_str(xyz, 0.0).unwrap();
        let params = NddoParameters::mndo().unwrap();
        let full = run_nddo_with_parameters(&molecule, &params, &NddoOptions::default()).unwrap();
        // 1 MiB forces the minimum depth and the A-DIIS  CDIIS degradation.
        let tight = run_nddo_with_parameters(
            &molecule,
            &params,
            &NddoOptions {
                scf_memory_mb: 1,
                ..NddoOptions::default()
            },
        )
        .unwrap();
        assert!(tight.converged);
        assert!(
            (full.total_ev - tight.total_ev).abs() < 1.0e-8,
            "full={} tight={}",
            full.total_ev,
            tight.total_ev
        );
    }

    #[test]
    fn pair_cache_budget_is_reported_not_aborted() {
        // 180 atoms  16 110 pairs  ~3 MiB of pair cache.
        let mut xyz = String::from("180\nwater lattice\n");
        for i in 0..60 {
            let (x, y, z) = (
                3.1 * (i % 5) as f64,
                3.1 * ((i / 5) % 4) as f64,
                3.1 * (i / 20) as f64,
            );
            xyz.push_str(&format!(
                "O {x} {y} {z}\nH {} {y} {z}\nH {} {} {z}\n",
                x + 0.96,
                x - 0.24,
                y + 0.93
            ));
        }
        let molecule = Molecule::from_xyz_str(&xyz, 0.0).unwrap();
        let params = NddoParameters::mndo().unwrap();
        let basis = Basis::build(&molecule, &params).unwrap();
        let required = crate::hamiltonian::pair_cache_bytes(&basis);
        assert!(required > 2 * 1024 * 1024, "cache estimate {required} B");

        // A budget below the requirement must fail cleanly, before allocating.
        let err = match crate::hamiltonian::build_core_limited(&molecule, &basis, &params, 1) {
            Err(err) => err,
            Ok(_) => panic!("a 1 MiB budget must reject a 3 MiB pair cache"),
        };
        let message = err.to_string();
        assert!(
            matches!(err, XndoError::ResourceLimit { .. }) && message.contains("pair cache"),
            "unexpected error: {message}"
        );
        // Zero disables the check; a large budget passes.
        assert!(crate::hamiltonian::build_core_limited(&molecule, &basis, &params, 0).is_ok());
        assert!(crate::hamiltonian::build_core_limited(&molecule, &basis, &params, 4096).is_ok());
    }

    #[test]
    fn molecule_charge_and_multiplicity_are_honored() {
        let xyz = "2\nh2+\nH 0.0 0.0 0.0\nH 0.74 0.0 0.0\n";
        let params = NddoParameters::mndo().unwrap();
        let molecule_metadata = Molecule::from_xyz_str(xyz, 1.0)
            .unwrap()
            .with_multiplicity(2);
        let from_molecule =
            run_nddo_with_parameters(&molecule_metadata, &params, &NddoOptions::default()).unwrap();

        let neutral_metadata = Molecule::from_xyz_str(xyz, 0.0).unwrap();
        let explicit = run_nddo_with_parameters(
            &neutral_metadata,
            &params,
            &NddoOptions {
                charge: 1.0,
                multiplicity: 2,
                ..NddoOptions::default()
            },
        )
        .unwrap();
        assert!(from_molecule.unrestricted);
        assert!((from_molecule.total_ev - explicit.total_ev).abs() < 1.0e-12);
        assert!((from_molecule.charges.iter().sum::<f64>() - 1.0).abs() < 1.0e-10);
    }

    #[test]
    fn conflicting_molecule_and_option_metadata_is_rejected() {
        let xyz = "2\nh2\nH 0.0 0.0 0.0\nH 0.74 0.0 0.0\n";
        let molecule = Molecule::from_xyz_str(xyz, 1.0).unwrap();
        let params = NddoParameters::mndo().unwrap();
        let err = run_nddo_with_parameters(
            &molecule,
            &params,
            &NddoOptions {
                charge: -1.0,
                multiplicity: 1,
                ..NddoOptions::default()
            },
        )
        .unwrap_err();
        assert!(err.to_string().contains("conflicting charge"));
    }

    #[test]
    fn water_heat_of_formation() {
        let xyz =
            "3\nwater\nO 0.0000 0.0000 0.0000\nH 0.9584 0.0000 0.0000\nH -0.2400 0.9278 0.0000\n";
        let r = run(xyz, 0.0);
        eprintln!(
            "H2O: dHf={:.3} kcal/mol  elec={:.4} eV core={:.4} eV  dipole={:.3} D  charges={:?}  iters={}",
            r.heat_of_formation_kcal, r.electronic_ev, r.core_ev, r.dipole_magnitude, r.charges, r.iterations
        );
        assert!(r.converged);
        // OpenMOPAC v23.2.5, MNDO PRECISE 1SCF, same geometry.
        assert!((r.heat_of_formation_kcal - (-60.5825049873711)).abs() < 0.005);
        assert!((r.dipole_magnitude - 1.79060758060072).abs() < 0.002);
        let qsum: f64 = r.charges.iter().sum();
        assert!(qsum.abs() < 1e-6);
        assert!(r.charges[0] < 0.0 && r.charges[1] > 0.0);
    }

    #[test]
    fn mndod_hcl_dipole_matches_openmopac() {
        let xyz = "2\nhydrogen chloride\nCl 0.0 0.0 0.0\nH 0.0 0.0 1.275\n";
        let molecule = Molecule::from_xyz_str(xyz, 0.0).unwrap();
        let params = NddoParameters::mndod().unwrap();
        let result = run_nddo_with_parameters(&molecule, &params, &NddoOptions::default()).unwrap();
        // OpenMOPAC v23.2.5, MNDOD PRECISE 1SCF, same geometry.
        assert!(
            (result.heat_of_formation_kcal - (-10.7256362989638)).abs() < 0.005,
            "heat of formation = {}",
            result.heat_of_formation_kcal
        );
        assert!(
            (result.dipole_magnitude - 1.012037782634).abs() < 0.002,
            "dipole = {}",
            result.dipole_magnitude
        );
    }

    #[test]
    fn water_no_diis_debug() {
        let xyz =
            "3\nwater\nO 0.0000 0.0000 0.0000\nH 0.9584 0.0000 0.0000\nH -0.2400 0.9278 0.0000\n";
        let mol = Molecule::from_xyz_str(xyz, 0.0).unwrap();
        let params = NddoParameters::mndo().unwrap();
        let opts = NddoOptions {
            use_diis: false,
            max_scf: 500,
            ..NddoOptions::default()
        };
        let r = run_nddo_with_parameters(&mol, &params, &opts).unwrap();
        eprintln!(
            "H2O(no-diis): dHf={:.3} conv={} iters={}",
            r.heat_of_formation_kcal, r.converged, r.iterations
        );
    }

    #[test]
    fn methane_heat_of_formation() {
        let xyz = "5\nmethane\nC 0.0000 0.0000 0.0000\nH 0.6276 0.6276 0.6276\nH -0.6276 -0.6276 0.6276\nH -0.6276 0.6276 -0.6276\nH 0.6276 -0.6276 -0.6276\n";
        let r = run(xyz, 0.0);
        eprintln!(
            "CH4: dHf={:.3} kcal/mol dipole={:.3} D iters={}",
            r.heat_of_formation_kcal, r.dipole_magnitude, r.iterations
        );
        assert!(r.converged);
    }

    #[test]
    fn accelerators_agree_on_energy() {
        // A-DIIS/CDIIS, plain CDIIS, and no acceleration must reach the same converged
        // energy (same SCF fixed point); the hybrid should not need more iterations.
        let xyz =
            "4\nformaldehyde\nC 0.0 0.0 0.0\nO 0.0 0.0 1.21\nH 0.94 0.0 -0.54\nH -0.94 0.0 -0.54\n";
        let mol = Molecule::from_xyz_str(xyz, 0.0).unwrap();
        let params = NddoParameters::mndo().unwrap();
        let run = |acc: ScfAccelerator| {
            let opts = NddoOptions {
                accelerator: acc,
                ..NddoOptions::default()
            };
            run_nddo_with_parameters(&mol, &params, &opts).unwrap()
        };
        let hybrid = run(ScfAccelerator::AdiisCdiis);
        let cdiis = run(ScfAccelerator::Cdiis);
        let none = run(ScfAccelerator::None);
        eprintln!(
            "iters: hybrid={} cdiis={} none={}  E(hybrid)={:.6}",
            hybrid.iterations, cdiis.iterations, none.iterations, hybrid.total_ev
        );
        assert!((hybrid.total_ev - cdiis.total_ev).abs() < 1e-6);
        assert!((hybrid.total_ev - none.total_ev).abs() < 1e-6);
        // Both accelerated paths should be far faster than plain iteration.
        assert!(hybrid.iterations < none.iterations);
    }

    #[test]
    fn methyl_radical_uhf() {
        // Planar CH3 radical (doublet): UHF must converge, be open-shell, and carry net spin.
        let xyz = "4\nmethyl\nC 0.0 0.0 0.0\nH 1.079 0.0 0.0\nH -0.5395 0.9344 0.0\nH -0.5395 -0.9344 0.0\n";
        let r = run_mult(xyz, 0.0, 2);
        eprintln!(
            "CH3.: dHf={:.3} kcal/mol unrestricted={} iters={}",
            r.heat_of_formation_kcal, r.unrestricted, r.iterations
        );
        assert!(r.converged);
        assert!(r.unrestricted);
        // Total spin population ( P_  P_) should be  1 unpaired electron.
        let spin = r.spin_density.as_ref().unwrap();
        let n_spin: f64 = (0..spin.rows).map(|i| spin[(i, i)]).sum();
        assert!((n_spin - 1.0).abs() < 1e-6, "net spin {n_spin}");
        // Frozen regression value at this geometry (r(C-H) = 1.079 A). The
        // MOPAC v23.2.5 UHF/DOUBLET oracle comparison for CH3 uses the
        // r = 1.078 A geometry and lives in `tests/molecules.rs`
        // (28.0055283157567 kcal/mol); see tools/oracle/MNDO_VALIDATION.md.
        assert!((r.heat_of_formation_kcal - 24.6110804799901).abs() < 0.005);
    }
}
