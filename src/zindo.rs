// SPDX-License-Identifier: GPL-3.0-or-later

//! Zerner INDO/S (ZINDO/S) spectroscopic Hamiltonian.
//!
//! The parameter table is extracted from OpenMOPAC v23.2.5's
//! `parameters_for_INDO_C.F90`.  OpenMOPAC's `INDO` keyword is INDO/S (also
//! known as ZINDO/S), not the original Pople ground-state INDO model.
//!
//! This implementation covers the ordinary 1-AO/4-AO s/p path: RHF/UHF ground
//! states, vertical spin-adapted singlet/triplet RHF-CIS, and spin-orbital UCIS.
//! In addition to energies and oscillator strengths, the spectrum API reports
//! transition/permanent dipoles, state-specific NDO charges, hole/particle
//! populations and centroids, and a simple charge-transfer distance.  The 9-AO
//! transition-metal/d-shell branch in OpenMOPAC uses additional irregular
//! Slater-Condon and double-zeta machinery and is deliberately rejected rather
//! than silently approximated. The s/p branch includes analytic RHF/UHF
//! ground-state gradients/Hessians and state-specific RHF-CIS/UCIS
//! gradients/Hessians.
//! Excited-state optimization and nonadiabatic derivative couplings are not
//! provided.

// AO, state, and Cartesian indices are intentionally explicit in the
// response-equation matrix contractions below.
#![allow(clippy::needless_range_loop)]

use std::collections::HashMap;
use std::sync::OnceLock;

use crate::data_tables::{self, CsvTable};
use crate::dual::{Dual, Scalar};
use crate::dual2::Dual2;
use crate::error::{Result, XndoError};
use crate::linalg::{symmetric_eigen, Matrix};
use crate::math::Vec3;
use crate::scf::Reference;
use crate::system::Molecule;
use crate::zdo_gradient::{
    atom_population, exchange_weight_rhf, exchange_weight_uhf, pair_energy, solve_rhf_responses,
    solve_rhf_second_response, solve_uhf_responses, solve_uhf_second_response, spin_densities,
    RhfResponse, UhfResponse, UhfSecondResponse,
};

/// Historical constants used by the OpenMOPAC INDO/S implementation.
/// These are kept method-local instead of using the NDDO CODATA constants so
/// the spectroscopic model follows the published/reference implementation.
pub const ZINDO_AU2EV: f64 = 27.2114;
pub const ZINDO_AU2ANG: f64 = 0.529_177;
pub const ZINDO_TOMK: f64 = 1.2;
pub const EV_NM: f64 = 1_239.841_984_332_002_6;
pub const EV_TO_WAVENUMBER_CM1: f64 = 8_065.544_005;
pub const DEBYE_PER_E_BOHR: f64 = 2.541_746_473;

/// OpenMOPAC stores FG14..FG21 multiplied by these integer convention factors
/// and divides them on entry to INDO/S (`switch.F90`).  FG22..FG24 are unscaled.
const FG_SCALE_14_24: [f64; 11] = [
    3.0, 25.0, 5.0, 15.0, 35.0, 245.0, 49.0, 441.0, 1.0, 1.0, 1.0,
];

/// Per-element Zerner INDO/S parameters.
#[derive(Clone, Debug)]
pub struct ZindoElement {
    pub z: u8,
    pub n_orb: usize,
    pub core_charge: f64,
    pub zeta_sp: f64,
    pub zeta_d: [f64; 2],
    pub zeta_weight: [f64; 2],
    pub beta: [f64; 3],
    /// Zerner FG parameters, index 1..=24 (index 0 unused).  FG14..FG21 are
    /// stored in their *runtime* OpenMOPAC convention after division by the
    /// integer factors documented above.
    pub fg: [f64; 25],
    pub n_s: u8,
    pub n_p: u8,
    pub n_d: u8,
    pub mass: f64,
}

impl ZindoElement {
    #[inline]
    fn beta_for_orb(&self, orb: usize) -> f64 {
        match orb {
            0 => self.beta[0],
            1..=3 => self.beta[1],
            _ => self.beta[2],
        }
    }

    #[inline]
    fn is_sp_supported(&self) -> bool {
        self.n_orb == 1 || self.n_orb == 4
    }
}

#[derive(Clone, Debug, Default)]
pub struct ZindoParameters {
    pub elements: HashMap<u8, ZindoElement>,
}

impl ZindoParameters {
    /// Return the process-wide parsed parameter table.
    pub fn cached() -> Result<&'static Self> {
        static PARAMS: OnceLock<ZindoParameters> = OnceLock::new();
        if let Some(params) = PARAMS.get() {
            return Ok(params);
        }
        let parsed = Self::standard()?;
        let _ = PARAMS.set(parsed);
        Ok(PARAMS.get().expect("ZINDO/S parameter cache initialized"))
    }

    pub fn standard() -> Result<Self> {
        let t = CsvTable::parse(data_tables::ZINDO_S_PARAM_CSV)
            .ok_or_else(|| XndoError::InvalidInput("empty ZINDO/S parameter table".into()))?;
        let col = |name: &str| -> Result<usize> {
            t.col(name).ok_or_else(|| {
                XndoError::MissingParameter(format!("ZINDO/S parameter column {name}"))
            })
        };
        let c_z = col("z")?;
        let c_norb = col("norb")?;
        let c_zcore = col("zcore")?;
        let c_zsp = col("zeta_sp")?;
        let c_zd1 = col("zeta_d1")?;
        let c_zd2 = col("zeta_d2")?;
        let c_zw1 = col("zeta_w1")?;
        let c_zw2 = col("zeta_w2")?;
        let c_bs = col("beta_s")?;
        let c_bp = col("beta_p")?;
        let c_bd = col("beta_d")?;
        let mut c_fg = [0usize; 25];
        for (i, slot) in c_fg.iter_mut().enumerate().skip(1) {
            *slot = col(&format!("fg{i}"))?;
        }
        let edata = data_tables::element_data();
        let mut elements = HashMap::new();
        for row in &t.rows {
            let z = t.f64_at(row, c_z) as u8;
            if z == 0 || z as usize >= edata.len() {
                continue;
            }
            let mut fg = [0.0; 25];
            for i in 1..=24 {
                fg[i] = t.f64_at(row, c_fg[i]);
            }
            // Mirror OpenMOPAC `switch`: the parameter source contains the
            // spectroscopic convention numerators for these integrals.
            for i in 14..=24 {
                fg[i] /= FG_SCALE_14_24[i - 14];
            }
            let ed = edata[z as usize];
            elements.insert(
                z,
                ZindoElement {
                    z,
                    n_orb: t.f64_at(row, c_norb) as usize,
                    core_charge: t.f64_at(row, c_zcore),
                    zeta_sp: t.f64_at(row, c_zsp),
                    zeta_d: [t.f64_at(row, c_zd1), t.f64_at(row, c_zd2)],
                    zeta_weight: [t.f64_at(row, c_zw1), t.f64_at(row, c_zw2)],
                    beta: [
                        t.f64_at(row, c_bs),
                        t.f64_at(row, c_bp),
                        t.f64_at(row, c_bd),
                    ],
                    fg,
                    n_s: ed.npq_s.max(1),
                    n_p: ed.npq_p.max(2),
                    n_d: ed.npq_d,
                    mass: ed.mass,
                },
            );
        }
        if elements.is_empty() {
            return Err(XndoError::InvalidInput(
                "no ZINDO/S elements parsed from embedded table".into(),
            ));
        }
        Ok(Self { elements })
    }

    pub fn element(&self, z: u8) -> Result<&ZindoElement> {
        self.elements.get(&z).ok_or(XndoError::MissingElement(z))
    }

    /// Elements in the embedded Zerner table that are currently enabled by the
    /// s/p implementation (H-like one-AO or s+p four-AO atoms).
    pub fn supported_sp_elements(&self) -> Vec<u8> {
        let mut z: Vec<u8> = self
            .elements
            .values()
            .filter(|e| e.is_sp_supported())
            .map(|e| e.z)
            .collect();
        z.sort_unstable();
        z
    }
}

#[derive(Clone, Debug)]
pub struct ZindoOptions {
    pub charge: f64,
    pub multiplicity: usize,
    pub reference: Reference,
    pub max_scf: usize,
    pub e_tol_ev: f64,
    pub p_tol: f64,
    /// Linear density damping: `P_next = (1 - damping) P_new + damping P_old`.
    pub damping: f64,
    pub n_states: usize,
    /// Optional number of occupied orbitals nearest the HOMO retained in CIS.
    pub active_occupied: Option<usize>,
    /// Optional number of virtual orbitals nearest the LUMO retained in CIS.
    pub active_virtual: Option<usize>,
}

impl Default for ZindoOptions {
    fn default() -> Self {
        Self {
            charge: 0.0,
            multiplicity: 1,
            reference: Reference::Auto,
            max_scf: 250,
            e_tol_ev: 1.0e-8,
            p_tol: 1.0e-7,
            damping: 0.20,
            n_states: 10,
            active_occupied: None,
            active_virtual: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct ZindoResult {
    pub density: Matrix,
    pub fock: Matrix,
    pub fock_beta: Option<Matrix>,
    pub mo_coeff: Matrix,
    pub mo_coeff_beta: Option<Matrix>,
    pub mo_energies_ev: Vec<f64>,
    pub mo_energies_beta_ev: Option<Vec<f64>>,
    pub n_occ: usize,
    pub n_alpha: usize,
    pub n_beta: usize,
    pub spin_density: Option<Matrix>,
    pub unrestricted: bool,
    pub electronic_ev: f64,
    pub core_ev: f64,
    pub total_ev: f64,
    pub charges: Vec<f64>,
    pub iterations: usize,
    pub converged: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CisSpin {
    Singlet,
    Triplet,
    /// Spin-orbital, spin-conserving CIS on an unrestricted reference.
    Unrestricted,
}

impl CisSpin {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Singlet => "singlet",
            Self::Triplet => "triplet",
            Self::Unrestricted => "unrestricted",
        }
    }

    pub fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "singlet" | "s" | "1" => Some(Self::Singlet),
            "triplet" | "t" | "3" => Some(Self::Triplet),
            "unrestricted" | "ucis" | "u" => Some(Self::Unrestricted),
            _ => None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct CiContribution {
    /// 1-based occupied MO index.
    pub occupied: usize,
    /// 1-based virtual MO index.
    pub virtual_orbital: usize,
    pub coefficient: f64,
}

#[derive(Clone, Debug)]
pub struct ExcitedState {
    /// Spin-adapted CIS sector used for this root.
    pub spin: CisSpin,
    /// Exact spin multiplicity of the spin-adapted CIS sector (1 or 3) and
    /// corresponding S(S+1) expectation value (0 or 2).
    pub spin_multiplicity: usize,
    pub s2_expectation: f64,
    /// Excitation energy relative to the selected ground-state reference.
    pub energy_ev: f64,
    /// Approximate absolute CIS-state electronic + core energy obtained as
    /// E_ground + omega.  This is useful for state-to-state bookkeeping; it
    /// does not imply an excited-state geometry relaxation.
    pub state_total_energy_ev: f64,
    pub state_total_energy_hartree: f64,
    pub energy_cm1: f64,
    pub wavelength_nm: Option<f64>,
    /// Electric-dipole oscillator strength. Triplet values are exactly zero in
    /// the nonrelativistic spin-free model.
    pub oscillator_strength: f64,
    /// Ground-to-excited transition dipole in e*bohr and Debye.
    pub transition_dipole_au: [f64; 3],
    pub transition_dipole_debye: [f64; 3],
    pub transition_dipole_magnitude_au: f64,
    pub dipole_strength_au2: f64,
    /// Excited-state permanent dipole from the CIS one-particle density.
    pub permanent_dipole_au: [f64; 3],
    pub permanent_dipole_debye: [f64; 3],
    pub permanent_dipole_magnitude_debye: f64,
    /// Excited-state minus ground-state permanent dipole.  Unlike an absolute
    /// dipole of a charged molecule, this difference is origin independent
    /// because both states carry the same total charge.
    pub difference_dipole_au: [f64; 3],
    pub difference_dipole_debye: [f64; 3],
    pub difference_dipole_magnitude_debye: f64,
    /// State-specific NDO net atomic charges.
    pub charges: Vec<f64>,
    /// Hole/electron populations associated with the CIS transition. Each set
    /// sums to approximately one electron for a normalized single-excitation
    /// vector in the orthogonal NDO AO basis.
    pub hole_population: Vec<f64>,
    pub electron_population: Vec<f64>,
    pub hole_centroid_angstrom: [f64; 3],
    pub electron_centroid_angstrom: [f64; 3],
    pub charge_transfer_distance_angstrom: f64,
    pub dominant: Vec<CiContribution>,
}

#[derive(Clone, Debug)]
pub struct ZindoSpectrum {
    pub ground: ZindoResult,
    pub spin: CisSpin,
    pub ground_dipole_au: [f64; 3],
    pub ground_dipole_debye: [f64; 3],
    pub states: Vec<ExcitedState>,
}

#[derive(Clone, Debug)]
pub struct ExcitedStateGradient {
    /// One-based root index in the selected spin-adapted CIS sector.
    pub root: usize,
    pub spin: CisSpin,
    pub excitation_energy_ev: f64,
    pub state_total_energy_ev: f64,
    /// Derivative of the vertical excitation energy, eV/Bohr.
    pub excitation_gradient: Vec<Vec3>,
    /// Derivative of `E_ground + excitation_energy`, eV/Bohr.
    pub state_gradient: Vec<Vec3>,
}

#[derive(Clone, Debug)]
pub struct ZindoCisGradientResult {
    pub ground: ZindoResult,
    pub ground_gradient: Vec<Vec3>,
    pub spin: CisSpin,
    pub states: Vec<ExcitedStateGradient>,
}

#[derive(Clone, Debug)]
pub struct ExcitedStateHessian {
    pub root: usize,
    pub spin: CisSpin,
    pub excitation_energy_ev: f64,
    pub state_total_energy_ev: f64,
    pub excitation_gradient: Vec<Vec3>,
    pub state_gradient: Vec<Vec3>,
    pub excitation_hessian: Matrix,
    pub state_hessian: Matrix,
}

#[derive(Clone, Debug)]
pub struct ZindoCisHessianResult {
    pub ground: ZindoResult,
    pub ground_gradient: Vec<Vec3>,
    pub ground_hessian: Matrix,
    pub spin: CisSpin,
    pub states: Vec<ExcitedStateHessian>,
}

#[derive(Clone, Copy, Debug)]
struct ZAo {
    atom: usize,
    orb: usize,
}

#[derive(Clone, Debug)]
struct ZBasis {
    aos: Vec<ZAo>,
    atom_offset: Vec<usize>,
    atom_norb: Vec<usize>,
}

impl ZBasis {
    fn build(mol: &Molecule, params: &ZindoParameters) -> Result<Self> {
        let mut aos = Vec::new();
        let mut atom_offset = Vec::with_capacity(mol.atoms.len());
        let mut atom_norb = Vec::with_capacity(mol.atoms.len());
        for (ia, atom) in mol.atoms.iter().enumerate() {
            let e = params.element(atom.z)?;
            if !e.is_sp_supported() {
                return Err(XndoError::InvalidInput(format!(
                    "ZINDO/S d-shell element Z={} requires OpenMOPAC's 9-AO transition-metal branch; xndo-rs v0.2.4 enables only the 1/4-AO s/p path",
                    atom.z
                )));
            }
            atom_offset.push(aos.len());
            atom_norb.push(e.n_orb);
            for orb in 0..e.n_orb {
                aos.push(ZAo { atom: ia, orb });
            }
        }
        Ok(Self {
            aos,
            atom_offset,
            atom_norb,
        })
    }

    fn nao(&self) -> usize {
        self.aos.len()
    }
}

/// OpenMOPAC INDO/S two-centre gamma function (eV), with distance in A.
#[inline]
pub fn gamma1(g1: f64, g2: f64, r_ang: f64) -> f64 {
    gamma1_g(g1, g2, r_ang)
}

#[inline]
fn gamma1_g<S: Scalar>(g1: f64, g2: f64, r_ang: S) -> S {
    if g1 <= 0.0 || g2 <= 0.0 {
        return S::cst(0.0);
    }
    let denom = r_ang / ZINDO_AU2ANG + ZINDO_TOMK * 2.0 * ZINDO_AU2EV / (g1 + g2);
    denom.recip() * (ZINDO_TOMK * ZINDO_AU2EV)
}

fn transpose4<S: Scalar>(a: [[S; 4]; 4]) -> [[S; 4]; 4] {
    let mut b = [[S::cst(0.0); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            b[j][i] = a[i][j];
        }
    }
    b
}

/// s/p Slater overlap block using the same orientation/sign convention as the
/// generic NDDO overlap kernel.  ZINDO/S uses one common exponent for s and p.
fn sp_overlap_block_g<S: Scalar>(ei: &ZindoElement, ej: &ZindoElement, d: [S; 3]) -> [[S; 4]; 4] {
    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if r.val() < 1.0e-14 {
        return [[S::cst(0.0); 4]; 4];
    }
    let inv_r = r.recip();
    let dir = [d[0] * inv_r, d[1] * inv_r, d[2] * inv_r];
    let (ea, eb, direction, swap) = if ei.n_s >= ej.n_s {
        (ei, ej, dir, false)
    } else {
        (ej, ei, [-dir[0], -dir[1], -dir[2]], true)
    };
    let (s111, s211, s121, s221, s222) = crate::overlap_numeric::slater_locals_numeric(
        ea.n_s, ea.n_p, ea.zeta_sp, ea.zeta_sp, eb.n_s, eb.n_p, eb.zeta_sp, eb.zeta_sp, r,
    );
    let di = crate::overlap::build_di_g::<S>([s111, s211, s121, s221, s222], direction);
    if swap {
        transpose4(di)
    } else {
        di
    }
}

fn sp_overlap_block(ei: &ZindoElement, ri: Vec3, ej: &ZindoElement, rj: Vec3) -> [[f64; 4]; 4] {
    let d = rj - ri;
    sp_overlap_block_g(ei, ej, [d.x, d.y, d.z])
}

fn build_overlap(mol: &Molecule, params: &ZindoParameters, basis: &ZBasis) -> Result<Matrix> {
    let n = basis.nao();
    let mut s = Matrix::identity(n);
    for ia in 0..mol.atoms.len() {
        let ei = params.element(mol.atoms[ia].z)?;
        let oi = basis.atom_offset[ia];
        let ni = basis.atom_norb[ia];
        for ja in (ia + 1)..mol.atoms.len() {
            let ej = params.element(mol.atoms[ja].z)?;
            let oj = basis.atom_offset[ja];
            let nj = basis.atom_norb[ja];
            let block = sp_overlap_block(ei, mol.atoms[ia].position, ej, mol.atoms[ja].position);
            for li in 0..ni {
                for lj in 0..nj {
                    s[(oi + li, oj + lj)] = block[li][lj];
                    s[(oj + lj, oi + li)] = block[li][lj];
                }
            }
        }
    }
    Ok(s)
}

/// Zerner U-form one-centre diagonal energy for the s/p path (eV).
fn uform_sp(e: &ZindoElement, orb: usize) -> f64 {
    let f0 = e.fg[11];
    let g1 = e.fg[14];
    let f2 = e.fg[15];
    let zc = e.core_charge.max(0.0);
    // The two spectroscopic atomic configurations and their mixing weights are
    // FG1/FG2 (s) or FG3/FG4 (p), with FG8/FG9 probabilities.  FG10 is the d
    // configuration weight and is irrelevant for the s/p branch.
    let weights = [e.fg[8], e.fg[9]];
    let wsum = (weights[0] + weights[1]).max(1.0e-30);
    let mut value = 0.0;
    if orb == 0 {
        let is = zc.min(2.0);
        let ip = zc - is;
        let ips = [e.fg[1], e.fg[2]];
        for k in 0..2 {
            if weights[k] == 0.0 {
                continue;
            }
            let u = ips[k] - (is - 1.0).max(0.0) * f0 - ip * (f0 - 0.5 * g1);
            value += weights[k] * u;
        }
    } else {
        let ip = (zc - 2.0).max(1.0);
        let is = zc - ip;
        let ipp = [e.fg[3], e.fg[4]];
        for k in 0..2 {
            if weights[k] == 0.0 {
                continue;
            }
            let u = ipp[k] - (ip - 1.0).max(0.0) * (f0 - 2.0 * f2) - is * (f0 - 0.5 * g1);
            value += weights[k] * u;
        }
    }
    value / wsum
}

/// Build the INDO/S Coulomb (J) and exchange (K) AO-pair kernels used by the
/// closed-shell SCF and by the CIS integral transformation.
fn build_jk(
    mol: &Molecule,
    params: &ZindoParameters,
    basis: &ZBasis,
) -> Result<(Matrix, Matrix, f64)> {
    let n = basis.nao();
    let mut j = Matrix::zeros(n, n);
    let mut k = Matrix::zeros(n, n);
    for mu in 0..n {
        let am = basis.aos[mu];
        let em = params.element(mol.atoms[am.atom].z)?;
        for nu in 0..n {
            let an = basis.aos[nu];
            let en = params.element(mol.atoms[an.atom].z)?;
            if am.atom != an.atom {
                let r_ang = (mol.atoms[am.atom].position - mol.atoms[an.atom].position).norm()
                    * ZINDO_AU2ANG;
                let g = gamma1(em.fg[11], en.fg[11], r_ang);
                j[(mu, nu)] = g;
                k[(mu, nu)] = g;
                continue;
            }
            let f0 = em.fg[11];
            let g1 = em.fg[14];
            let f2 = em.fg[15];
            match (am.orb, an.orb) {
                (a, b) if a == b => {
                    j[(mu, nu)] = if a == 0 { f0 } else { f0 + 4.0 * f2 };
                    k[(mu, nu)] = j[(mu, nu)];
                }
                (0, 1..=3) | (1..=3, 0) => {
                    j[(mu, nu)] = f0;
                    k[(mu, nu)] = g1;
                }
                (1..=3, 1..=3) => {
                    j[(mu, nu)] = f0 - 2.0 * f2;
                    k[(mu, nu)] = 3.0 * f2;
                }
                _ => {}
            }
        }
    }
    let mut ecore = 0.0;
    for ia in 0..mol.atoms.len() {
        let ei = params.element(mol.atoms[ia].z)?;
        for ja in (ia + 1)..mol.atoms.len() {
            let ej = params.element(mol.atoms[ja].z)?;
            let r_ang = (mol.atoms[ia].position - mol.atoms[ja].position).norm() * ZINDO_AU2ANG;
            ecore += ei.core_charge * ej.core_charge * gamma1(ei.fg[11], ej.fg[11], r_ang);
        }
    }
    Ok((j, k, ecore))
}

fn build_core(
    mol: &Molecule,
    params: &ZindoParameters,
    basis: &ZBasis,
    overlap: &Matrix,
) -> Result<Matrix> {
    let n = basis.nao();
    let mut h = Matrix::zeros(n, n);
    for mu in 0..n {
        let a = basis.aos[mu];
        let ea = params.element(mol.atoms[a.atom].z)?;
        let mut diag = uform_sp(ea, a.orb);
        for ib in 0..mol.atoms.len() {
            if ib == a.atom {
                continue;
            }
            let eb = params.element(mol.atoms[ib].z)?;
            let r_ang = (mol.atoms[a.atom].position - mol.atoms[ib].position).norm() * ZINDO_AU2ANG;
            diag -= eb.core_charge * gamma1(ea.fg[11], eb.fg[11], r_ang);
        }
        h[(mu, mu)] = diag;
        for nu in 0..mu {
            let b = basis.aos[nu];
            let eb = params.element(mol.atoms[b.atom].z)?;
            let hij = 0.5 * (ea.beta_for_orb(a.orb) + eb.beta_for_orb(b.orb)) * overlap[(mu, nu)];
            h[(mu, nu)] = hij;
            h[(nu, mu)] = hij;
        }
    }
    Ok(h)
}

fn build_fock_closed(h: &Matrix, p: &Matrix, j: &Matrix, k: &Matrix, basis: &ZBasis) -> Matrix {
    let n = basis.nao();
    let mut f = h.clone();
    for mu in 0..n {
        let am = basis.aos[mu];
        let mut d = h[(mu, mu)];
        for nu in 0..n {
            let an = basis.aos[nu];
            if mu == nu {
                d += 0.5 * p[(nu, nu)] * j[(mu, nu)];
            } else if am.atom == an.atom {
                d += p[(nu, nu)] * (j[(mu, nu)] - 0.5 * k[(mu, nu)]);
            } else {
                d += p[(nu, nu)] * j[(mu, nu)];
            }
        }
        f[(mu, mu)] = d;
        for nu in 0..mu {
            let an = basis.aos[nu];
            let v = if am.atom == an.atom {
                h[(mu, nu)] + (1.5 * k[(mu, nu)] - 0.5 * j[(mu, nu)]) * p[(mu, nu)]
            } else {
                h[(mu, nu)] - 0.5 * k[(mu, nu)] * p[(mu, nu)]
            };
            f[(mu, nu)] = v;
            f[(nu, mu)] = v;
        }
    }
    f
}

/// Spin-resolved INDO/S Fock matrices. The equations reduce exactly to the
/// closed-shell form when `P_alpha = P_beta = P / 2`.
fn build_fock_uhf(
    h: &Matrix,
    pa: &Matrix,
    pb: &Matrix,
    j: &Matrix,
    k: &Matrix,
    basis: &ZBasis,
) -> (Matrix, Matrix) {
    let n = basis.nao();
    let pt = add_matrix(pa, pb);
    let mut fa = h.clone();
    let mut fb = h.clone();
    for mu in 0..n {
        let am = basis.aos[mu];
        let mut da = h[(mu, mu)];
        let mut db = h[(mu, mu)];
        for nu in 0..n {
            let an = basis.aos[nu];
            if mu == nu || am.atom == an.atom {
                da += pt[(nu, nu)] * j[(mu, nu)] - pa[(nu, nu)] * k[(mu, nu)];
                db += pt[(nu, nu)] * j[(mu, nu)] - pb[(nu, nu)] * k[(mu, nu)];
            } else {
                da += pt[(nu, nu)] * j[(mu, nu)];
                db += pt[(nu, nu)] * j[(mu, nu)];
            }
        }
        fa[(mu, mu)] = da;
        fb[(mu, mu)] = db;
        for nu in 0..mu {
            let an = basis.aos[nu];
            let (va, vb) = if am.atom == an.atom {
                (
                    h[(mu, nu)] + 2.0 * k[(mu, nu)] * pt[(mu, nu)]
                        - (j[(mu, nu)] + k[(mu, nu)]) * pa[(mu, nu)],
                    h[(mu, nu)] + 2.0 * k[(mu, nu)] * pt[(mu, nu)]
                        - (j[(mu, nu)] + k[(mu, nu)]) * pb[(mu, nu)],
                )
            } else {
                (
                    h[(mu, nu)] - k[(mu, nu)] * pa[(mu, nu)],
                    h[(mu, nu)] - k[(mu, nu)] * pb[(mu, nu)],
                )
            };
            fa[(mu, nu)] = va;
            fa[(nu, mu)] = va;
            fb[(mu, nu)] = vb;
            fb[(nu, mu)] = vb;
        }
    }
    (fa, fb)
}

fn atomic_density(mol: &Molecule, params: &ZindoParameters, basis: &ZBasis) -> Result<Matrix> {
    let n = basis.nao();
    let mut p = Matrix::zeros(n, n);
    for ia in 0..mol.atoms.len() {
        let e = params.element(mol.atoms[ia].z)?;
        let o = basis.atom_offset[ia];
        if e.n_orb == 1 {
            p[(o, o)] = e.core_charge.min(2.0);
        } else {
            let ns = e.core_charge.min(2.0);
            p[(o, o)] = ns;
            let np = (e.core_charge - ns).max(0.0) / 3.0;
            for q in 1..=3 {
                p[(o + q, o + q)] = np;
            }
        }
    }
    Ok(p)
}

fn electron_count(mol: &Molecule, params: &ZindoParameters, charge: f64) -> Result<usize> {
    let q = charge.round();
    if (charge - q).abs() > 1.0e-8 {
        return Err(XndoError::InvalidInput(
            "ZINDO/S charge must be an integer number of electrons".into(),
        ));
    }
    let mut ne = -q as i64;
    for atom in &mol.atoms {
        ne += params.element(atom.z)?.core_charge.round() as i64;
    }
    if ne < 0 {
        return Err(XndoError::InvalidInput(
            "negative ZINDO/S electron count".into(),
        ));
    }
    Ok(ne as usize)
}

fn electron_partition(ne: usize, multiplicity: usize) -> Result<(usize, usize)> {
    if multiplicity == 0 {
        return Err(XndoError::InvalidInput(
            "ZINDO/S multiplicity must be at least one".into(),
        ));
    }
    let unpaired = multiplicity - 1;
    if unpaired > ne || !(ne - unpaired).is_multiple_of(2) {
        return Err(XndoError::InvalidInput(format!(
            "ZINDO/S electron count {ne} is incompatible with multiplicity {multiplicity}"
        )));
    }
    let n_beta = (ne - unpaired) / 2;
    Ok((n_beta + unpaired, n_beta))
}

fn damp_density(raw: &Matrix, old: &Matrix, damping: f64) -> Matrix {
    if damping <= 0.0 {
        return raw.clone();
    }
    let mut result = raw.clone();
    for index in 0..result.as_slice().len() {
        result.as_mut_slice()[index] =
            (1.0 - damping) * raw.as_slice()[index] + damping * old.as_slice()[index];
    }
    result
}

fn sub_matrix(a: &Matrix, b: &Matrix) -> Matrix {
    let mut result = Matrix::zeros(a.rows, a.cols);
    for index in 0..result.as_slice().len() {
        result.as_mut_slice()[index] = a.as_slice()[index] - b.as_slice()[index];
    }
    result
}

fn charges(
    mol: &Molecule,
    params: &ZindoParameters,
    basis: &ZBasis,
    p: &Matrix,
) -> Result<Vec<f64>> {
    let mut out = Vec::with_capacity(mol.atoms.len());
    for ia in 0..mol.atoms.len() {
        let e = params.element(mol.atoms[ia].z)?;
        let o = basis.atom_offset[ia];
        let mut pop = 0.0;
        for q in 0..basis.atom_norb[ia] {
            pop += p[(o + q, o + q)];
        }
        out.push(e.core_charge - pop);
    }
    Ok(out)
}

/// Zerner INDO/S RHF or UHF SCF for fixed molecular geometry.
pub fn run_zindo_s(
    mol: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
) -> Result<ZindoResult> {
    let basis = ZBasis::build(mol, params)?;
    let n = basis.nao();
    if n == 0 {
        return Err(XndoError::InvalidInput(
            "ZINDO/S molecule has no basis functions".into(),
        ));
    }
    let ne = electron_count(mol, params, options.charge)?;
    let multiplicity = if options.multiplicity == 1 && mol.multiplicity != 1 {
        mol.multiplicity
    } else {
        options.multiplicity
    };
    let (n_alpha, n_beta) = electron_partition(ne, multiplicity)?;
    if n_alpha > n || n_beta > n {
        return Err(XndoError::InvalidInput(format!(
            "ZINDO/S has {ne} electrons but only {n} spatial orbitals"
        )));
    }
    let reference = match options.reference {
        Reference::Auto => {
            if n_alpha == n_beta {
                Reference::Rhf
            } else {
                Reference::Uhf
            }
        }
        value => value,
    };
    if reference == Reference::Rhf && n_alpha != n_beta {
        return Err(XndoError::InvalidInput(format!(
            "ZINDO/S RHF is incompatible with {n_alpha} alpha and {n_beta} beta electrons; use UHF or Auto"
        )));
    }
    let overlap = build_overlap(mol, params, &basis)?;
    let h = build_core(mol, params, &basis, &overlap)?;
    let (j, k, core_ev) = build_jk(mol, params, &basis)?;
    let guess = atomic_density(mol, params, &basis)?;
    let damping = options.damping.clamp(0.0, 0.95);
    let mut last_e = f64::INFINITY;
    let mut last_err = f64::INFINITY;
    let (_, initial_coefficients) = symmetric_eigen(&guess)?;
    if reference == Reference::Rhf {
        let n_occ = n_alpha;
        let mut p = initial_coefficients.leading_columns_gram(n_occ, 2.0);
        for iter in 1..=options.max_scf {
            let f = build_fock_closed(&h, &p, &j, &k, &basis);
            let (_, c) = symmetric_eigen(&f)?;
            let p_next = damp_density(&c.leading_columns_gram(n_occ, 2.0), &p, damping);
            let f_next = build_fock_closed(&h, &p_next, &j, &k, &basis);
            let e = 0.5 * p_next.frobenius_dot(&add_matrix(&h, &f_next));
            let p_err = p_next.rms_difference(&p);
            let e_err = (e - last_e).abs();
            p = p_next;
            last_e = e;
            last_err = p_err;
            if p_err < options.p_tol && e_err < options.e_tol_ev {
                let f0 = build_fock_closed(&h, &p, &j, &k, &basis);
                let (eps, c) = symmetric_eigen(&f0)?;
                let p_final = c.leading_columns_gram(n_occ, 2.0);
                let f_final = build_fock_closed(&h, &p_final, &j, &k, &basis);
                let electronic_ev = 0.5 * p_final.frobenius_dot(&add_matrix(&h, &f_final));
                return Ok(ZindoResult {
                    density: p_final.clone(),
                    fock: f_final,
                    fock_beta: None,
                    mo_coeff: c,
                    mo_coeff_beta: None,
                    mo_energies_ev: eps,
                    mo_energies_beta_ev: None,
                    n_occ,
                    n_alpha,
                    n_beta,
                    spin_density: None,
                    unrestricted: false,
                    electronic_ev,
                    core_ev,
                    total_ev: electronic_ev + core_ev,
                    charges: charges(mol, params, &basis, &p_final)?,
                    iterations: iter,
                    converged: true,
                });
            }
        }
    } else {
        let mut pa = initial_coefficients.leading_columns_gram(n_alpha, 1.0);
        let mut pb = initial_coefficients.leading_columns_gram(n_beta, 1.0);
        for iter in 1..=options.max_scf {
            let (fa, fb) = build_fock_uhf(&h, &pa, &pb, &j, &k, &basis);
            let (_, ca) = symmetric_eigen(&fa)?;
            let (_, cb) = symmetric_eigen(&fb)?;
            let pa_next = damp_density(&ca.leading_columns_gram(n_alpha, 1.0), &pa, damping);
            let pb_next = damp_density(&cb.leading_columns_gram(n_beta, 1.0), &pb, damping);
            let (fa_next, fb_next) = build_fock_uhf(&h, &pa_next, &pb_next, &j, &k, &basis);
            let total_next = add_matrix(&pa_next, &pb_next);
            let e = 0.5
                * (total_next.frobenius_dot(&h)
                    + pa_next.frobenius_dot(&fa_next)
                    + pb_next.frobenius_dot(&fb_next));
            let p_err = pa_next.rms_difference(&pa).max(pb_next.rms_difference(&pb));
            let e_err = (e - last_e).abs();
            pa = pa_next;
            pb = pb_next;
            last_e = e;
            last_err = p_err;
            if p_err < options.p_tol && e_err < options.e_tol_ev {
                let (fa0, fb0) = build_fock_uhf(&h, &pa, &pb, &j, &k, &basis);
                let (eps_a, ca) = symmetric_eigen(&fa0)?;
                let (eps_b, cb) = symmetric_eigen(&fb0)?;
                let pa_final = ca.leading_columns_gram(n_alpha, 1.0);
                let pb_final = cb.leading_columns_gram(n_beta, 1.0);
                let density = add_matrix(&pa_final, &pb_final);
                let (fa_final, fb_final) = build_fock_uhf(&h, &pa_final, &pb_final, &j, &k, &basis);
                let electronic_ev = 0.5
                    * (density.frobenius_dot(&h)
                        + pa_final.frobenius_dot(&fa_final)
                        + pb_final.frobenius_dot(&fb_final));
                return Ok(ZindoResult {
                    density: density.clone(),
                    fock: fa_final,
                    fock_beta: Some(fb_final),
                    mo_coeff: ca,
                    mo_coeff_beta: Some(cb),
                    mo_energies_ev: eps_a,
                    mo_energies_beta_ev: Some(eps_b),
                    n_occ: n_alpha,
                    n_alpha,
                    n_beta,
                    spin_density: Some(sub_matrix(&pa_final, &pb_final)),
                    unrestricted: true,
                    electronic_ev,
                    core_ev,
                    total_ev: electronic_ev + core_ev,
                    charges: charges(mol, params, &basis, &density)?,
                    iterations: iter,
                    converged: true,
                });
            }
        }
    }
    Err(XndoError::ScfNotConverged {
        iterations: options.max_scf,
        error: last_err,
    })
}

/// Analytic Cartesian RHF/UHF ground-state gradient in eV/Bohr.
pub(crate) fn analytic_ground_gradient(
    mol: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
) -> Result<(ZindoResult, Vec<Vec3>)> {
    let result = run_zindo_s(mol, params, options)?;
    let basis = ZBasis::build(mol, params)?;
    let spin = result
        .spin_density
        .as_ref()
        .map(|density| spin_densities(&result.density, density));
    let mut gradient = vec![Vec3::zero(); mol.atoms.len()];
    for ia in 0..mol.atoms.len() {
        let ea = params.element(mol.atoms[ia].z)?;
        let oa = basis.atom_offset[ia];
        let na = basis.atom_norb[ia];
        let pop_a = atom_population(&result.density, oa, na);
        for ib in (ia + 1)..mol.atoms.len() {
            let eb = params.element(mol.atoms[ib].z)?;
            let ob = basis.atom_offset[ib];
            let nb = basis.atom_norb[ib];
            let pop_b = atom_population(&result.density, ob, nb);
            let d = mol.atoms[ib].position - mol.atoms[ia].position;
            if d.norm() < 1.0e-12 {
                return Err(XndoError::InvalidInput(
                    "coincident atoms in ZINDO/S".into(),
                ));
            }
            let dv = [Dual::var(d.x, 0), Dual::var(d.y, 1), Dual::var(d.z, 2)];
            let r_bohr = (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]).sqrt();
            let gamma = gamma1_g(ea.fg[11], eb.fg[11], r_bohr * ZINDO_AU2ANG);
            let overlap = sp_overlap_block_g(ea, eb, dv);
            let mut resonance = Dual::constant(0.0);
            for i in 0..na {
                for j in 0..nb {
                    let beta = 0.5 * (ea.beta_for_orb(i) + eb.beta_for_orb(j));
                    resonance =
                        resonance + overlap[i][j] * (result.density[(oa + i, ob + j)] * beta);
                }
            }
            let exchange = match &spin {
                Some((alpha, beta)) => exchange_weight_uhf(alpha, beta, oa, na, ob, nb),
                None => exchange_weight_rhf(&result.density, oa, na, ob, nb),
            };
            let core = gamma * (ea.core_charge * eb.core_charge);
            let pair_derivative = pair_energy(
                gamma,
                core,
                ea.core_charge,
                eb.core_charge,
                pop_a,
                pop_b,
                exchange,
                resonance,
            );
            let dg = Vec3::new(
                pair_derivative.d[0],
                pair_derivative.d[1],
                pair_derivative.d[2],
            );
            gradient[ia] -= dg;
            gradient[ib] += dg;
        }
    }
    Ok((result, gradient))
}

struct ZindoPairDerivatives {
    ia: usize,
    ib: usize,
    gamma: Dual2,
    overlap: [[Dual2; 4]; 4],
    core: Dual2,
    exchange: f64,
    population_a: f64,
    population_b: f64,
}

struct ZindoResponseData {
    basis: ZBasis,
    pairs: Vec<ZindoPairDerivatives>,
    responses: Vec<RhfResponse>,
}

struct ZindoUhfResponseData {
    basis: ZBasis,
    pairs: Vec<ZindoPairDerivatives>,
    responses: Vec<UhfResponse>,
    density_alpha: Matrix,
    density_beta: Matrix,
}

fn zindo_uhf_response_data(
    mol: &Molecule,
    params: &ZindoParameters,
    ground: &ZindoResult,
) -> Result<ZindoUhfResponseData> {
    if !ground.unrestricted {
        return Err(XndoError::InvalidInput(
            "UCIS response requires a ZINDO/S UHF reference".into(),
        ));
    }
    let basis = ZBasis::build(mol, params)?;
    let (j0, k0, _) = build_jk(mol, params, &basis)?;
    let spin = ground
        .spin_density
        .as_ref()
        .expect("UHF spin density must be present");
    let (density_alpha, density_beta) = spin_densities(&ground.density, spin);
    let ncoord = 3 * mol.atoms.len();
    let nao = basis.nao();
    let zero = Matrix::zeros(nao, nao);
    let mut skeleton_alpha = (0..ncoord)
        .map(|_| Matrix::zeros(nao, nao))
        .collect::<Vec<_>>();
    let mut skeleton_beta = (0..ncoord)
        .map(|_| Matrix::zeros(nao, nao))
        .collect::<Vec<_>>();
    let mut pairs = Vec::new();
    for ia in 0..mol.atoms.len() {
        let ea = params.element(mol.atoms[ia].z)?;
        let oa = basis.atom_offset[ia];
        let na = basis.atom_norb[ia];
        let population_a = atom_population(&ground.density, oa, na);
        for ib in (ia + 1)..mol.atoms.len() {
            let eb = params.element(mol.atoms[ib].z)?;
            let ob = basis.atom_offset[ib];
            let nb = basis.atom_norb[ib];
            let population_b = atom_population(&ground.density, ob, nb);
            let displacement = mol.atoms[ib].position - mol.atoms[ia].position;
            let coordinates = [
                Dual2::var(displacement.x, 0),
                Dual2::var(displacement.y, 1),
                Dual2::var(displacement.z, 2),
            ];
            let distance = (coordinates[0] * coordinates[0]
                + coordinates[1] * coordinates[1]
                + coordinates[2] * coordinates[2])
                .sqrt();
            let gamma = gamma1_g(ea.fg[11], eb.fg[11], distance * ZINDO_AU2ANG);
            let overlap = sp_overlap_block_g(ea, eb, coordinates);
            for axis in 0..3 {
                let dg = gamma.g[axis];
                for &(atom, sign) in &[(ia, -1.0), (ib, 1.0)] {
                    let coordinate = 3 * atom + axis;
                    for i in 0..na {
                        let value = sign * (-eb.core_charge + population_b) * dg;
                        skeleton_alpha[coordinate][(oa + i, oa + i)] += value;
                        skeleton_beta[coordinate][(oa + i, oa + i)] += value;
                    }
                    for q in 0..nb {
                        let value = sign * (-ea.core_charge + population_a) * dg;
                        skeleton_alpha[coordinate][(ob + q, ob + q)] += value;
                        skeleton_beta[coordinate][(ob + q, ob + q)] += value;
                    }
                    for i in 0..na {
                        for q in 0..nb {
                            let mu = oa + i;
                            let nu = ob + q;
                            let beta = 0.5 * (ea.beta_for_orb(i) + eb.beta_for_orb(q));
                            let resonance = beta * overlap[i][q].g[axis];
                            let alpha = resonance - density_alpha[(mu, nu)] * dg;
                            let beta_spin = resonance - density_beta[(mu, nu)] * dg;
                            skeleton_alpha[coordinate][(mu, nu)] += sign * alpha;
                            skeleton_alpha[coordinate][(nu, mu)] += sign * alpha;
                            skeleton_beta[coordinate][(mu, nu)] += sign * beta_spin;
                            skeleton_beta[coordinate][(nu, mu)] += sign * beta_spin;
                        }
                    }
                }
            }
            pairs.push(ZindoPairDerivatives {
                ia,
                ib,
                gamma,
                overlap,
                core: gamma * (ea.core_charge * eb.core_charge),
                exchange: exchange_weight_uhf(&density_alpha, &density_beta, oa, na, ob, nb),
                population_a,
                population_b,
            });
        }
    }
    let coefficients_beta = ground
        .mo_coeff_beta
        .as_ref()
        .expect("UHF beta coefficients");
    let energies_beta = ground
        .mo_energies_beta_ev
        .as_ref()
        .expect("UHF beta energies");
    let responses = solve_uhf_responses(
        &ground.mo_coeff,
        coefficients_beta,
        &ground.mo_energies_ev,
        energies_beta,
        ground.n_alpha,
        ground.n_beta,
        &skeleton_alpha,
        &skeleton_beta,
        |dpa, dpb| Ok(build_fock_uhf(&zero, dpa, dpb, &j0, &k0, &basis)),
    )?;
    Ok(ZindoUhfResponseData {
        basis,
        pairs,
        responses,
        density_alpha,
        density_beta,
    })
}

fn zindo_response_data(
    mol: &Molecule,
    params: &ZindoParameters,
    ground: &ZindoResult,
) -> Result<ZindoResponseData> {
    let basis = ZBasis::build(mol, params)?;
    let (j0, k0, _) = build_jk(mol, params, &basis)?;
    let ncoord = 3 * mol.atoms.len();
    let nao = basis.nao();
    let zero = Matrix::zeros(nao, nao);
    let mut skeleton = (0..ncoord)
        .map(|_| Matrix::zeros(nao, nao))
        .collect::<Vec<_>>();
    let mut pairs = Vec::new();
    for ia in 0..mol.atoms.len() {
        let ea = params.element(mol.atoms[ia].z)?;
        let oa = basis.atom_offset[ia];
        let na = basis.atom_norb[ia];
        let pop_a = atom_population(&ground.density, oa, na);
        for ib in (ia + 1)..mol.atoms.len() {
            let eb = params.element(mol.atoms[ib].z)?;
            let ob = basis.atom_offset[ib];
            let nb = basis.atom_norb[ib];
            let pop_b = atom_population(&ground.density, ob, nb);
            let d = mol.atoms[ib].position - mol.atoms[ia].position;
            let dv = [Dual2::var(d.x, 0), Dual2::var(d.y, 1), Dual2::var(d.z, 2)];
            let r_bohr = (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]).sqrt();
            let gamma = gamma1_g(ea.fg[11], eb.fg[11], r_bohr * ZINDO_AU2ANG);
            let overlap = sp_overlap_block_g(ea, eb, dv);
            for axis in 0..3 {
                let dg = gamma.g[axis];
                for &(atom, sign) in &[(ia, -1.0), (ib, 1.0)] {
                    let coord = 3 * atom + axis;
                    for i in 0..na {
                        skeleton[coord][(oa + i, oa + i)] += sign * (-eb.core_charge + pop_b) * dg;
                    }
                    for j in 0..nb {
                        skeleton[coord][(ob + j, ob + j)] += sign * (-ea.core_charge + pop_a) * dg;
                    }
                    for i in 0..na {
                        for j in 0..nb {
                            let mu = oa + i;
                            let nu = ob + j;
                            let beta = 0.5 * (ea.beta_for_orb(i) + eb.beta_for_orb(j));
                            let value =
                                beta * overlap[i][j].g[axis] - 0.5 * ground.density[(mu, nu)] * dg;
                            skeleton[coord][(mu, nu)] += sign * value;
                            skeleton[coord][(nu, mu)] += sign * value;
                        }
                    }
                }
            }
            pairs.push(ZindoPairDerivatives {
                ia,
                ib,
                gamma,
                overlap,
                core: gamma * (ea.core_charge * eb.core_charge),
                exchange: exchange_weight_rhf(&ground.density, oa, na, ob, nb),
                population_a: pop_a,
                population_b: pop_b,
            });
        }
    }
    let responses = solve_rhf_responses(
        &ground.mo_coeff,
        &ground.mo_energies_ev,
        ground.n_occ,
        &skeleton,
        |dp| Ok(build_fock_closed(&zero, dp, &j0, &k0, &basis)),
    )?;
    Ok(ZindoResponseData {
        basis,
        pairs,
        responses,
    })
}

fn pair_coordinate(pair: &ZindoPairDerivatives, coordinate: usize) -> Option<(usize, f64)> {
    let atom = coordinate / 3;
    let axis = coordinate % 3;
    if atom == pair.ia {
        Some((axis, -1.0))
    } else if atom == pair.ib {
        Some((axis, 1.0))
    } else {
        None
    }
}

fn zindo_jk_derivative(data: &ZindoResponseData, coordinate: usize) -> (Matrix, Matrix) {
    let mut jx = Matrix::zeros(data.basis.nao(), data.basis.nao());
    let mut kx = Matrix::zeros(data.basis.nao(), data.basis.nao());
    for pair in &data.pairs {
        let Some((axis, sign)) = pair_coordinate(pair, coordinate) else {
            continue;
        };
        let value = sign * pair.gamma.g[axis];
        let oa = data.basis.atom_offset[pair.ia];
        let ob = data.basis.atom_offset[pair.ib];
        for i in 0..data.basis.atom_norb[pair.ia] {
            for q in 0..data.basis.atom_norb[pair.ib] {
                let mu = oa + i;
                let nu = ob + q;
                jx[(mu, nu)] = value;
                jx[(nu, mu)] = value;
                kx[(mu, nu)] = value;
                kx[(nu, mu)] = value;
            }
        }
    }
    (jx, kx)
}

fn zindo_jk_second_derivative(
    data: &ZindoResponseData,
    coordinate_x: usize,
    coordinate_y: usize,
) -> (Matrix, Matrix) {
    let mut jxy = Matrix::zeros(data.basis.nao(), data.basis.nao());
    let mut kxy = Matrix::zeros(data.basis.nao(), data.basis.nao());
    for pair in &data.pairs {
        let (Some((axis_x, sign_x)), Some((axis_y, sign_y))) = (
            pair_coordinate(pair, coordinate_x),
            pair_coordinate(pair, coordinate_y),
        ) else {
            continue;
        };
        let value = sign_x * sign_y * pair.gamma.h[axis_x][axis_y];
        let oa = data.basis.atom_offset[pair.ia];
        let ob = data.basis.atom_offset[pair.ib];
        for i in 0..data.basis.atom_norb[pair.ia] {
            for q in 0..data.basis.atom_norb[pair.ib] {
                let mu = oa + i;
                let nu = ob + q;
                jxy[(mu, nu)] = value;
                jxy[(nu, mu)] = value;
                kxy[(mu, nu)] = value;
                kxy[(nu, mu)] = value;
            }
        }
    }
    (jxy, kxy)
}

fn zindo_known_fock_second(
    mol: &Molecule,
    params: &ZindoParameters,
    ground: &ZindoResult,
    data: &ZindoResponseData,
    coordinate_x: usize,
    coordinate_y: usize,
) -> Result<Matrix> {
    let basis = &data.basis;
    let dpx = &data.responses[coordinate_x].density;
    let dpy = &data.responses[coordinate_y].density;
    let mut fxy = Matrix::zeros(basis.nao(), basis.nao());
    for pair in &data.pairs {
        let ea = params.element(mol.atoms[pair.ia].z)?;
        let eb = params.element(mol.atoms[pair.ib].z)?;
        let oa = basis.atom_offset[pair.ia];
        let ob = basis.atom_offset[pair.ib];
        let na = basis.atom_norb[pair.ia];
        let nb = basis.atom_norb[pair.ib];
        let x = pair_coordinate(pair, coordinate_x);
        let y = pair_coordinate(pair, coordinate_y);
        let dgx = x
            .map(|(axis, sign)| sign * pair.gamma.g[axis])
            .unwrap_or(0.0);
        let dgy = y
            .map(|(axis, sign)| sign * pair.gamma.g[axis])
            .unwrap_or(0.0);
        let dgxy = match (x, y) {
            (Some((ax, sx)), Some((ay, sy))) => sx * sy * pair.gamma.h[ax][ay],
            _ => 0.0,
        };
        let dpop_ax = atom_population(dpx, oa, na);
        let dpop_bx = atom_population(dpx, ob, nb);
        let dpop_ay = atom_population(dpy, oa, na);
        let dpop_by = atom_population(dpy, ob, nb);
        for i in 0..na {
            fxy[(oa + i, oa + i)] +=
                (-eb.core_charge + pair.population_b) * dgxy + dpop_by * dgx + dpop_bx * dgy;
        }
        for q in 0..nb {
            fxy[(ob + q, ob + q)] +=
                (-ea.core_charge + pair.population_a) * dgxy + dpop_ay * dgx + dpop_ax * dgy;
        }
        for i in 0..na {
            for q in 0..nb {
                let mu = oa + i;
                let nu = ob + q;
                let beta = 0.5 * (ea.beta_for_orb(i) + eb.beta_for_orb(q));
                let dsxy = match (x, y) {
                    (Some((ax, sx)), Some((ay, sy))) => sx * sy * pair.overlap[i][q].h[ax][ay],
                    _ => 0.0,
                };
                let value = beta * dsxy
                    - 0.5 * ground.density[(mu, nu)] * dgxy
                    - 0.5 * (dpy[(mu, nu)] * dgx + dpx[(mu, nu)] * dgy);
                fxy[(mu, nu)] += value;
                fxy[(nu, mu)] += value;
            }
        }
    }
    Ok(fxy)
}

fn zindo_known_fock_second_uhf(
    mol: &Molecule,
    params: &ZindoParameters,
    data: &ZindoUhfResponseData,
    coordinate_x: usize,
    coordinate_y: usize,
) -> Result<(Matrix, Matrix)> {
    let basis = &data.basis;
    let response_x = &data.responses[coordinate_x];
    let response_y = &data.responses[coordinate_y];
    let dpt_x = add_matrix(&response_x.density_alpha, &response_x.density_beta);
    let dpt_y = add_matrix(&response_y.density_alpha, &response_y.density_beta);
    let mut alpha = Matrix::zeros(basis.nao(), basis.nao());
    let mut beta_spin = Matrix::zeros(basis.nao(), basis.nao());
    for pair in &data.pairs {
        let ea = params.element(mol.atoms[pair.ia].z)?;
        let eb = params.element(mol.atoms[pair.ib].z)?;
        let oa = basis.atom_offset[pair.ia];
        let ob = basis.atom_offset[pair.ib];
        let na = basis.atom_norb[pair.ia];
        let nb = basis.atom_norb[pair.ib];
        let x = pair_coordinate(pair, coordinate_x);
        let y = pair_coordinate(pair, coordinate_y);
        let dgx = x
            .map(|(axis, sign)| sign * pair.gamma.g[axis])
            .unwrap_or(0.0);
        let dgy = y
            .map(|(axis, sign)| sign * pair.gamma.g[axis])
            .unwrap_or(0.0);
        let dgxy = match (x, y) {
            (Some((axis_x, sign_x)), Some((axis_y, sign_y))) => {
                sign_x * sign_y * pair.gamma.h[axis_x][axis_y]
            }
            _ => 0.0,
        };
        let dpop_ax = atom_population(&dpt_x, oa, na);
        let dpop_bx = atom_population(&dpt_x, ob, nb);
        let dpop_ay = atom_population(&dpt_y, oa, na);
        let dpop_by = atom_population(&dpt_y, ob, nb);
        for i in 0..na {
            let value =
                (-eb.core_charge + pair.population_b) * dgxy + dpop_by * dgx + dpop_bx * dgy;
            alpha[(oa + i, oa + i)] += value;
            beta_spin[(oa + i, oa + i)] += value;
        }
        for q in 0..nb {
            let value =
                (-ea.core_charge + pair.population_a) * dgxy + dpop_ay * dgx + dpop_ax * dgy;
            alpha[(ob + q, ob + q)] += value;
            beta_spin[(ob + q, ob + q)] += value;
        }
        for i in 0..na {
            for q in 0..nb {
                let mu = oa + i;
                let nu = ob + q;
                let resonance = 0.5 * (ea.beta_for_orb(i) + eb.beta_for_orb(q));
                let dsxy = match (x, y) {
                    (Some((axis_x, sign_x)), Some((axis_y, sign_y))) => {
                        sign_x * sign_y * pair.overlap[i][q].h[axis_x][axis_y]
                    }
                    _ => 0.0,
                };
                let alpha_value = resonance * dsxy
                    - data.density_alpha[(mu, nu)] * dgxy
                    - response_y.density_alpha[(mu, nu)] * dgx
                    - response_x.density_alpha[(mu, nu)] * dgy;
                let beta_value = resonance * dsxy
                    - data.density_beta[(mu, nu)] * dgxy
                    - response_y.density_beta[(mu, nu)] * dgx
                    - response_x.density_beta[(mu, nu)] * dgy;
                alpha[(mu, nu)] += alpha_value;
                alpha[(nu, mu)] += alpha_value;
                beta_spin[(mu, nu)] += beta_value;
                beta_spin[(nu, mu)] += beta_value;
            }
        }
    }
    Ok((alpha, beta_spin))
}

fn analytic_ground_hessian_uhf(
    mol: &Molecule,
    params: &ZindoParameters,
    ground: ZindoResult,
) -> Result<(ZindoResult, Vec<Vec3>, Matrix)> {
    if mol.atoms.len() == 2 {
        return Err(XndoError::LinearAlgebra(
            "an open-shell diatomic state-specific UHF Hessian is not unique on a degenerate electronic surface; use a specified state-averaged or diabatic model"
                .into(),
        ));
    }
    let basis = ZBasis::build(mol, params)?;
    let (j0, k0, _) = build_jk(mol, params, &basis)?;
    let ncoord = 3 * mol.atoms.len();
    let nao = basis.nao();
    let zero = Matrix::zeros(nao, nao);
    let spin_density = ground
        .spin_density
        .as_ref()
        .expect("UHF ZINDO/S spin density");
    let (pa, pb) = spin_densities(&ground.density, spin_density);
    let mut skeleton_a = (0..ncoord)
        .map(|_| Matrix::zeros(nao, nao))
        .collect::<Vec<_>>();
    let mut skeleton_b = (0..ncoord)
        .map(|_| Matrix::zeros(nao, nao))
        .collect::<Vec<_>>();
    let mut pairs = Vec::new();
    for ia in 0..mol.atoms.len() {
        let ea = params.element(mol.atoms[ia].z)?;
        let oa = basis.atom_offset[ia];
        let na = basis.atom_norb[ia];
        let pop_a = atom_population(&ground.density, oa, na);
        for ib in (ia + 1)..mol.atoms.len() {
            let eb = params.element(mol.atoms[ib].z)?;
            let ob = basis.atom_offset[ib];
            let nb = basis.atom_norb[ib];
            let pop_b = atom_population(&ground.density, ob, nb);
            let displacement = mol.atoms[ib].position - mol.atoms[ia].position;
            let coordinates = [
                Dual2::var(displacement.x, 0),
                Dual2::var(displacement.y, 1),
                Dual2::var(displacement.z, 2),
            ];
            let distance = (coordinates[0] * coordinates[0]
                + coordinates[1] * coordinates[1]
                + coordinates[2] * coordinates[2])
                .sqrt();
            let gamma = gamma1_g(ea.fg[11], eb.fg[11], distance * ZINDO_AU2ANG);
            let overlap = sp_overlap_block_g(ea, eb, coordinates);
            for axis in 0..3 {
                let dg = gamma.g[axis];
                for &(atom, sign) in &[(ia, -1.0), (ib, 1.0)] {
                    let coordinate = 3 * atom + axis;
                    for i in 0..na {
                        let value = sign * (-eb.core_charge + pop_b) * dg;
                        skeleton_a[coordinate][(oa + i, oa + i)] += value;
                        skeleton_b[coordinate][(oa + i, oa + i)] += value;
                    }
                    for q in 0..nb {
                        let value = sign * (-ea.core_charge + pop_a) * dg;
                        skeleton_a[coordinate][(ob + q, ob + q)] += value;
                        skeleton_b[coordinate][(ob + q, ob + q)] += value;
                    }
                    for i in 0..na {
                        for q in 0..nb {
                            let mu = oa + i;
                            let nu = ob + q;
                            let beta = 0.5 * (ea.beta_for_orb(i) + eb.beta_for_orb(q));
                            let resonance = beta * overlap[i][q].g[axis];
                            let alpha = resonance - pa[(mu, nu)] * dg;
                            let beta_spin = resonance - pb[(mu, nu)] * dg;
                            skeleton_a[coordinate][(mu, nu)] += sign * alpha;
                            skeleton_a[coordinate][(nu, mu)] += sign * alpha;
                            skeleton_b[coordinate][(mu, nu)] += sign * beta_spin;
                            skeleton_b[coordinate][(nu, mu)] += sign * beta_spin;
                        }
                    }
                }
            }
            pairs.push(ZindoPairDerivatives {
                ia,
                ib,
                gamma,
                overlap,
                core: gamma * (ea.core_charge * eb.core_charge),
                exchange: exchange_weight_uhf(&pa, &pb, oa, na, ob, nb),
                population_a: pop_a,
                population_b: pop_b,
            });
        }
    }
    let cb = ground
        .mo_coeff_beta
        .as_ref()
        .expect("UHF beta coefficients");
    let eps_b = ground
        .mo_energies_beta_ev
        .as_ref()
        .expect("UHF beta orbital energies");
    let responses = solve_uhf_responses(
        &ground.mo_coeff,
        cb,
        &ground.mo_energies_ev,
        eps_b,
        ground.n_alpha,
        ground.n_beta,
        &skeleton_a,
        &skeleton_b,
        |dpa, dpb| Ok(build_fock_uhf(&zero, dpa, dpb, &j0, &k0, &basis)),
    )?;
    let mut gradient = vec![Vec3::zero(); mol.atoms.len()];
    let mut hessian = Matrix::zeros(ncoord, ncoord);
    for pair in &pairs {
        let ea = params.element(mol.atoms[pair.ia].z)?;
        let eb = params.element(mol.atoms[pair.ib].z)?;
        let oa = basis.atom_offset[pair.ia];
        let ob = basis.atom_offset[pair.ib];
        let na = basis.atom_norb[pair.ia];
        let nb = basis.atom_norb[pair.ib];
        let mut resonance = Dual2::constant(0.0);
        for i in 0..na {
            for q in 0..nb {
                let beta = 0.5 * (ea.beta_for_orb(i) + eb.beta_for_orb(q));
                resonance =
                    resonance + pair.overlap[i][q] * (ground.density[(oa + i, ob + q)] * beta);
            }
        }
        let fixed = pair_energy(
            pair.gamma,
            pair.core,
            ea.core_charge,
            eb.core_charge,
            pair.population_a,
            pair.population_b,
            pair.exchange,
            resonance,
        );
        let pair_gradient = Vec3::new(fixed.g[0], fixed.g[1], fixed.g[2]);
        gradient[pair.ia] -= pair_gradient;
        gradient[pair.ib] += pair_gradient;
        for row_axis in 0..3 {
            for column_axis in 0..3 {
                for &(row_atom, row_sign) in &[(pair.ia, -1.0), (pair.ib, 1.0)] {
                    for &(column_atom, column_sign) in &[(pair.ia, -1.0), (pair.ib, 1.0)] {
                        hessian[(3 * row_atom + row_axis, 3 * column_atom + column_axis)] +=
                            row_sign * column_sign * fixed.h[row_axis][column_axis];
                    }
                }
            }
        }
        for (coordinate, response) in responses.iter().enumerate() {
            let dpt = add_matrix(&response.density_alpha, &response.density_beta);
            let dpop_a = atom_population(&dpt, oa, na);
            let dpop_b = atom_population(&dpt, ob, nb);
            let mut dexchange = 0.0;
            for i in 0..na {
                for q in 0..nb {
                    let mu = oa + i;
                    let nu = ob + q;
                    dexchange += 2.0
                        * (pa[(mu, nu)] * response.density_alpha[(mu, nu)]
                            + pb[(mu, nu)] * response.density_beta[(mu, nu)]);
                }
            }
            let dweight = -eb.core_charge * dpop_a - ea.core_charge * dpop_b
                + dpop_a * pair.population_b
                + pair.population_a * dpop_b
                - dexchange;
            for axis in 0..3 {
                let mut value = pair.gamma.g[axis] * dweight;
                for i in 0..na {
                    for q in 0..nb {
                        let beta = 0.5 * (ea.beta_for_orb(i) + eb.beta_for_orb(q));
                        value += 2.0 * dpt[(oa + i, ob + q)] * beta * pair.overlap[i][q].g[axis];
                    }
                }
                hessian[(3 * pair.ia + axis, coordinate)] -= value;
                hessian[(3 * pair.ib + axis, coordinate)] += value;
            }
        }
    }
    for row in 0..ncoord {
        for column in 0..row {
            let value = 0.5 * (hessian[(row, column)] + hessian[(column, row)]);
            hessian[(row, column)] = value;
            hessian[(column, row)] = value;
        }
    }
    Ok((ground, gradient, hessian))
}

/// Fully analytic Cartesian ZINDO/S RHF/UHF ground-state Hessian in eV/Bohr^2.
pub fn analytic_ground_hessian(
    mol: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
) -> Result<(ZindoResult, Vec<Vec3>, Matrix)> {
    let ground = run_zindo_s(mol, params, options)?;
    if ground.unrestricted {
        return analytic_ground_hessian_uhf(mol, params, ground);
    }
    let response_data = zindo_response_data(mol, params, &ground)?;
    let basis = &response_data.basis;
    let ncoord = 3 * mol.atoms.len();
    let mut gradient = vec![Vec3::zero(); mol.atoms.len()];
    let mut hessian = Matrix::zeros(ncoord, ncoord);
    for pair_data in &response_data.pairs {
        let ea = params.element(mol.atoms[pair_data.ia].z)?;
        let eb = params.element(mol.atoms[pair_data.ib].z)?;
        let oa = basis.atom_offset[pair_data.ia];
        let ob = basis.atom_offset[pair_data.ib];
        let na = basis.atom_norb[pair_data.ia];
        let nb = basis.atom_norb[pair_data.ib];
        let mut resonance = Dual2::constant(0.0);
        for i in 0..na {
            for j in 0..nb {
                let beta = 0.5 * (ea.beta_for_orb(i) + eb.beta_for_orb(j));
                resonance =
                    resonance + pair_data.overlap[i][j] * (ground.density[(oa + i, ob + j)] * beta);
            }
        }
        let fixed = pair_energy(
            pair_data.gamma,
            pair_data.core,
            ea.core_charge,
            eb.core_charge,
            pair_data.population_a,
            pair_data.population_b,
            pair_data.exchange,
            resonance,
        );
        let gpair = Vec3::new(fixed.g[0], fixed.g[1], fixed.g[2]);
        gradient[pair_data.ia] -= gpair;
        gradient[pair_data.ib] += gpair;
        for q in 0..3 {
            for r in 0..3 {
                for &(row_atom, row_sign) in &[(pair_data.ia, -1.0), (pair_data.ib, 1.0)] {
                    for &(col_atom, col_sign) in &[(pair_data.ia, -1.0), (pair_data.ib, 1.0)] {
                        hessian[(3 * row_atom + q, 3 * col_atom + r)] +=
                            row_sign * col_sign * fixed.h[q][r];
                    }
                }
            }
        }
        for (coord, response) in response_data.responses.iter().enumerate() {
            let dpt = &response.density;
            let dpop_a = atom_population(dpt, oa, na);
            let dpop_b = atom_population(dpt, ob, nb);
            let mut dexchange = 0.0;
            for i in 0..na {
                for j in 0..nb {
                    let mu = oa + i;
                    let nu = ob + j;
                    dexchange += ground.density[(mu, nu)] * dpt[(mu, nu)];
                }
            }
            let dweight = -eb.core_charge * dpop_a - ea.core_charge * dpop_b
                + dpop_a * pair_data.population_b
                + pair_data.population_a * dpop_b
                - dexchange;
            for q in 0..3 {
                let mut response_term = pair_data.gamma.g[q] * dweight;
                for i in 0..na {
                    for j in 0..nb {
                        let beta = 0.5 * (ea.beta_for_orb(i) + eb.beta_for_orb(j));
                        response_term +=
                            2.0 * dpt[(oa + i, ob + j)] * beta * pair_data.overlap[i][j].g[q];
                    }
                }
                hessian[(3 * pair_data.ia + q, coord)] -= response_term;
                hessian[(3 * pair_data.ib + q, coord)] += response_term;
            }
        }
    }
    for i in 0..ncoord {
        for j in 0..i {
            let value = 0.5 * (hessian[(i, j)] + hessian[(j, i)]);
            hessian[(i, j)] = value;
            hessian[(j, i)] = value;
        }
    }
    Ok((ground, gradient, hessian))
}

fn add_matrix(a: &Matrix, b: &Matrix) -> Matrix {
    debug_assert_eq!(a.rows, b.rows);
    debug_assert_eq!(a.cols, b.cols);
    let mut out = Matrix::zeros(a.rows, a.cols);
    for i in 0..out.as_slice().len() {
        out.as_mut_slice()[i] = a.as_slice()[i] + b.as_slice()[i];
    }
    out
}

/// NDO AO integral transformed to four MOs.  Only `(|)` and `(|)`
/// survive; J and K contain those two classes respectively.
fn mo_eri(c: &Matrix, j: &Matrix, k: &Matrix, p: usize, q: usize, r: usize, s: usize) -> f64 {
    let n = c.rows;
    let mut v = 0.0;
    for mu in 0..n {
        let pq = c[(mu, p)] * c[(mu, q)];
        if pq == 0.0 {
            continue;
        }
        for nu in 0..n {
            v += pq * c[(nu, r)] * c[(nu, s)] * j[(mu, nu)];
        }
    }
    // The exchange-like surviving integrals have  != .  Grouping the four
    // symmetry-related ordered products in these two symmetric coefficients is
    // both cheaper and less error-prone than enumerating the AO quartet.
    for mu in 0..n {
        for nu in (mu + 1)..n {
            let a = c[(mu, p)] * c[(nu, q)] + c[(nu, p)] * c[(mu, q)];
            let b = c[(mu, r)] * c[(nu, s)] + c[(nu, r)] * c[(mu, s)];
            v += k[(mu, nu)] * a * b;
        }
    }
    v
}

fn mo_eri_mixed(coefficients: [&Matrix; 4], j: &Matrix, k: &Matrix, orbitals: [usize; 4]) -> f64 {
    let n = coefficients[0].rows;
    let [p, q, r, s] = orbitals;
    let mut result = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            result += coefficients[0][(mu, p)]
                * coefficients[1][(mu, q)]
                * coefficients[2][(nu, r)]
                * coefficients[3][(nu, s)]
                * j[(mu, nu)];
        }
    }
    for mu in 0..n {
        for nu in (mu + 1)..n {
            for ao in [
                [mu, nu, mu, nu],
                [mu, nu, nu, mu],
                [nu, mu, mu, nu],
                [nu, mu, nu, mu],
            ] {
                result += coefficients[0][(ao[0], p)]
                    * coefficients[1][(ao[1], q)]
                    * coefficients[2][(ao[2], r)]
                    * coefficients[3][(ao[3], s)]
                    * k[(mu, nu)];
            }
        }
    }
    result
}

#[allow(clippy::too_many_arguments)]
fn mo_eri_mixed_derivative(
    coefficients: [&Matrix; 4],
    derivatives: [&Matrix; 4],
    j: &Matrix,
    jx: &Matrix,
    k: &Matrix,
    kx: &Matrix,
    orbitals: [usize; 4],
) -> f64 {
    let n = coefficients[0].rows;
    let mut result = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            let ao = [mu, mu, nu, nu];
            let values: [f64; 4] =
                std::array::from_fn(|index| coefficients[index][(ao[index], orbitals[index])]);
            let xs: [f64; 4] =
                std::array::from_fn(|index| derivatives[index][(ao[index], orbitals[index])]);
            let zeros = [0.0; 4];
            let (value, x, _, _) = product_derivatives(&values, &xs, &zeros, &zeros);
            result += x * j[(mu, nu)] + value * jx[(mu, nu)];
        }
    }
    for mu in 0..n {
        for nu in (mu + 1)..n {
            for ao in [
                [mu, nu, mu, nu],
                [mu, nu, nu, mu],
                [nu, mu, mu, nu],
                [nu, mu, nu, mu],
            ] {
                let values: [f64; 4] =
                    std::array::from_fn(|index| coefficients[index][(ao[index], orbitals[index])]);
                let xs: [f64; 4] =
                    std::array::from_fn(|index| derivatives[index][(ao[index], orbitals[index])]);
                let zeros = [0.0; 4];
                let (value, x, _, _) = product_derivatives(&values, &xs, &zeros, &zeros);
                result += x * k[(mu, nu)] + value * kx[(mu, nu)];
            }
        }
    }
    result
}

#[allow(clippy::too_many_arguments)]
fn mo_eri_mixed_second_derivative(
    coefficients: [&Matrix; 4],
    derivatives_x: [&Matrix; 4],
    derivatives_y: [&Matrix; 4],
    derivatives_xy: [&Matrix; 4],
    j: &Matrix,
    jx: &Matrix,
    jy: &Matrix,
    jxy: &Matrix,
    k: &Matrix,
    kx: &Matrix,
    ky: &Matrix,
    kxy: &Matrix,
    orbitals: [usize; 4],
) -> f64 {
    let n = coefficients[0].rows;
    let contribution = |ao: [usize; 4], integral: f64, ix: f64, iy: f64, ixy: f64| {
        let values: [f64; 4] =
            std::array::from_fn(|index| coefficients[index][(ao[index], orbitals[index])]);
        let xs: [f64; 4] =
            std::array::from_fn(|index| derivatives_x[index][(ao[index], orbitals[index])]);
        let ys: [f64; 4] =
            std::array::from_fn(|index| derivatives_y[index][(ao[index], orbitals[index])]);
        let xys: [f64; 4] =
            std::array::from_fn(|index| derivatives_xy[index][(ao[index], orbitals[index])]);
        let (value, x, y, xy) = product_derivatives(&values, &xs, &ys, &xys);
        xy * integral + x * iy + y * ix + value * ixy
    };
    let mut result = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            result += contribution(
                [mu, mu, nu, nu],
                j[(mu, nu)],
                jx[(mu, nu)],
                jy[(mu, nu)],
                jxy[(mu, nu)],
            );
        }
    }
    for mu in 0..n {
        for nu in (mu + 1)..n {
            for ao in [
                [mu, nu, mu, nu],
                [mu, nu, nu, mu],
                [nu, mu, mu, nu],
                [nu, mu, nu, mu],
            ] {
                result += contribution(ao, k[(mu, nu)], kx[(mu, nu)], ky[(mu, nu)], kxy[(mu, nu)]);
            }
        }
    }
    result
}

#[allow(clippy::too_many_arguments)]
fn mo_eri_derivative(
    c: &Matrix,
    cx: &Matrix,
    j: &Matrix,
    jx: &Matrix,
    k: &Matrix,
    kx: &Matrix,
    p: usize,
    q: usize,
    r: usize,
    s: usize,
) -> f64 {
    let n = c.rows;
    let mut value = 0.0;
    for mu in 0..n {
        let pq = c[(mu, p)] * c[(mu, q)];
        let pqx = cx[(mu, p)] * c[(mu, q)] + c[(mu, p)] * cx[(mu, q)];
        for nu in 0..n {
            let rs = c[(nu, r)] * c[(nu, s)];
            let rsx = cx[(nu, r)] * c[(nu, s)] + c[(nu, r)] * cx[(nu, s)];
            value += pqx * rs * j[(mu, nu)] + pq * rsx * j[(mu, nu)] + pq * rs * jx[(mu, nu)];
        }
    }
    for mu in 0..n {
        for nu in (mu + 1)..n {
            let a = c[(mu, p)] * c[(nu, q)] + c[(nu, p)] * c[(mu, q)];
            let ax = cx[(mu, p)] * c[(nu, q)]
                + c[(mu, p)] * cx[(nu, q)]
                + cx[(nu, p)] * c[(mu, q)]
                + c[(nu, p)] * cx[(mu, q)];
            let b = c[(mu, r)] * c[(nu, s)] + c[(nu, r)] * c[(mu, s)];
            let bx = cx[(mu, r)] * c[(nu, s)]
                + c[(mu, r)] * cx[(nu, s)]
                + cx[(nu, r)] * c[(mu, s)]
                + c[(nu, r)] * cx[(mu, s)];
            value += kx[(mu, nu)] * a * b + k[(mu, nu)] * (ax * b + a * bx);
        }
    }
    value
}

fn product_derivatives(value: &[f64], dx: &[f64], dy: &[f64], dxy: &[f64]) -> (f64, f64, f64, f64) {
    let n = value.len();
    let product_except = |skip_a: usize, skip_b: usize| {
        let mut product = 1.0;
        for (idx, &item) in value.iter().enumerate() {
            if idx != skip_a && idx != skip_b {
                product *= item;
            }
        }
        product
    };
    let product = value.iter().product::<f64>();
    let mut x = 0.0;
    let mut y = 0.0;
    let mut xy = 0.0;
    for i in 0..n {
        let rest = product_except(i, i);
        x += dx[i] * rest;
        y += dy[i] * rest;
        xy += dxy[i] * rest;
        for j in 0..n {
            if i != j {
                xy += dx[i] * dy[j] * product_except(i, j);
            }
        }
    }
    (product, x, y, xy)
}

#[allow(clippy::too_many_arguments)]
fn mo_eri_second_derivative(
    c: &Matrix,
    cx: &Matrix,
    cy: &Matrix,
    cxy: &Matrix,
    j: &Matrix,
    jx: &Matrix,
    jy: &Matrix,
    jxy: &Matrix,
    k: &Matrix,
    kx: &Matrix,
    ky: &Matrix,
    kxy: &Matrix,
    p: usize,
    q: usize,
    r: usize,
    s: usize,
) -> f64 {
    let n = c.rows;
    let mut result = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            let indices = [(mu, p), (mu, q), (nu, r), (nu, s)];
            let values = indices.map(|idx| c[idx]);
            let xs = indices.map(|idx| cx[idx]);
            let ys = indices.map(|idx| cy[idx]);
            let xys = indices.map(|idx| cxy[idx]);
            let (v, vx, vy, vxy) = product_derivatives(&values, &xs, &ys, &xys);
            result += vxy * j[(mu, nu)] + vx * jy[(mu, nu)] + vy * jx[(mu, nu)] + v * jxy[(mu, nu)];
        }
    }
    for mu in 0..n {
        for nu in (mu + 1)..n {
            let pair = |left: usize, right: usize| {
                let indices = [(mu, left), (nu, right)];
                let values = indices.map(|idx| c[idx]);
                let xs = indices.map(|idx| cx[idx]);
                let ys = indices.map(|idx| cy[idx]);
                let xys = indices.map(|idx| cxy[idx]);
                product_derivatives(&values, &xs, &ys, &xys)
            };
            let a1 = pair(p, q);
            let a2 = {
                let indices = [(nu, p), (mu, q)];
                product_derivatives(
                    &indices.map(|idx| c[idx]),
                    &indices.map(|idx| cx[idx]),
                    &indices.map(|idx| cy[idx]),
                    &indices.map(|idx| cxy[idx]),
                )
            };
            let b1 = pair(r, s);
            let b2 = {
                let indices = [(nu, r), (mu, s)];
                product_derivatives(
                    &indices.map(|idx| c[idx]),
                    &indices.map(|idx| cx[idx]),
                    &indices.map(|idx| cy[idx]),
                    &indices.map(|idx| cxy[idx]),
                )
            };
            let a = (a1.0 + a2.0, a1.1 + a2.1, a1.2 + a2.2, a1.3 + a2.3);
            let b = (b1.0 + b2.0, b1.1 + b2.1, b1.2 + b2.2, b1.3 + b2.3);
            let ab = a.0 * b.0;
            let abx = a.1 * b.0 + a.0 * b.1;
            let aby = a.2 * b.0 + a.0 * b.2;
            let abxy = a.3 * b.0 + a.1 * b.2 + a.2 * b.1 + a.0 * b.3;
            result +=
                kxy[(mu, nu)] * ab + kx[(mu, nu)] * aby + ky[(mu, nu)] * abx + k[(mu, nu)] * abxy;
        }
    }
    result
}

fn ao_position_matrices(
    mol: &Molecule,
    params: &ZindoParameters,
    basis: &ZBasis,
) -> Result<[Matrix; 3]> {
    let n = basis.nao();
    let mut r = [
        Matrix::zeros(n, n),
        Matrix::zeros(n, n),
        Matrix::zeros(n, n),
    ];
    for mu in 0..n {
        let a = basis.aos[mu];
        let pos = mol.atoms[a.atom].position;
        r[0][(mu, mu)] = pos.x;
        r[1][(mu, mu)] = pos.y;
        r[2][(mu, mu)] = pos.z;
    }
    for ia in 0..mol.atoms.len() {
        let e = params.element(mol.atoms[ia].z)?;
        if e.n_orb < 4 || e.zeta_sp <= 0.0 {
            continue;
        }
        // OpenMOPAC `dipole`: <ns|r|np> = (2n+1)/(2 zeta sqrt(3)) in Bohr.
        let rsp = (2.0 * e.n_s as f64 + 1.0) / (2.0 * e.zeta_sp * 3.0_f64.sqrt());
        let o = basis.atom_offset[ia];
        for axis in 0..3 {
            let pidx = o + 1 + axis;
            r[axis][(o, pidx)] = rsp;
            r[axis][(pidx, o)] = rsp;
        }
    }
    Ok(r)
}

fn mo_position(c: &Matrix, r: &Matrix, i: usize, a: usize) -> f64 {
    let n = c.rows;
    let mut v = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            v += c[(mu, i)] * r[(mu, nu)] * c[(nu, a)];
        }
    }
    v
}

fn mo_to_ao_density(c: &Matrix, p_mo: &Matrix) -> Matrix {
    // AO density = C P_MO C^T. An explicit contraction keeps this helper
    // independent of storage/layout choices in the linalg wrapper.
    let nao = c.rows;
    let nmo = c.cols;
    let mut p_ao = Matrix::zeros(nao, nao);
    for mu in 0..nao {
        for nu in 0..nao {
            let mut v = 0.0;
            for i in 0..nmo {
                let cmi = c[(mu, i)];
                if cmi == 0.0 {
                    continue;
                }
                for j in 0..nmo {
                    v += cmi * p_mo[(i, j)] * c[(nu, j)];
                }
            }
            p_ao[(mu, nu)] = v;
        }
    }
    p_ao
}

fn permanent_dipole_au(
    mol: &Molecule,
    params: &ZindoParameters,
    basis: &ZBasis,
    density: &Matrix,
    r_ao: &[Matrix; 3],
) -> Result<[f64; 3]> {
    let mut mu = [0.0; 3];
    for (ia, atom) in mol.atoms.iter().enumerate() {
        let zc = params.element(atom.z)?.core_charge;
        let r = mol.atoms[ia].position;
        mu[0] += zc * r.x;
        mu[1] += zc * r.y;
        mu[2] += zc * r.z;
    }
    let n = basis.nao();
    for axis in 0..3 {
        let mut electronic = 0.0;
        for mu_idx in 0..n {
            for nu_idx in 0..n {
                electronic += density[(mu_idx, nu_idx)] * r_ao[axis][(nu_idx, mu_idx)];
            }
        }
        mu[axis] -= electronic;
    }
    Ok(mu)
}

fn scale3(v: [f64; 3], factor: f64) -> [f64; 3] {
    [v[0] * factor, v[1] * factor, v[2] * factor]
}

fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn atom_diagonal_population(basis: &ZBasis, p: &Matrix) -> Vec<f64> {
    let mut out = vec![0.0; basis.atom_offset.len()];
    for ia in 0..basis.atom_offset.len() {
        let o = basis.atom_offset[ia];
        for q in 0..basis.atom_norb[ia] {
            out[ia] += p[(o + q, o + q)];
        }
    }
    out
}

fn population_centroid_angstrom(mol: &Molecule, population: &[f64]) -> [f64; 3] {
    let total = population.iter().sum::<f64>();
    if total.abs() < 1.0e-14 {
        return [f64::NAN; 3];
    }
    let mut c = [0.0; 3];
    for (atom, &w) in mol.atoms.iter().zip(population.iter()) {
        c[0] += w * atom.position.x * ZINDO_AU2ANG;
        c[1] += w * atom.position.y * ZINDO_AU2ANG;
        c[2] += w * atom.position.z * ZINDO_AU2ANG;
    }
    [c[0] / total, c[1] / total, c[2] / total]
}

fn state_density_and_attachment(
    c: &Matrix,
    singles: &[(usize, usize)],
    x: &Matrix,
    state: usize,
    nocc: usize,
) -> (Matrix, Matrix, Matrix) {
    let nmo = c.cols;
    let mut hole_mo = Matrix::zeros(nmo, nmo);
    let mut particle_mo = Matrix::zeros(nmo, nmo);
    for pidx in 0..singles.len() {
        let (i, a) = singles[pidx];
        let cp = x[(pidx, state)];
        for qidx in 0..singles.len() {
            let (j, b) = singles[qidx];
            let w = cp * x[(qidx, state)];
            if a == b {
                hole_mo[(i, j)] += w;
            }
            if i == j {
                particle_mo[(a, b)] += w;
            }
        }
    }
    let mut state_mo = Matrix::zeros(nmo, nmo);
    for i in 0..nocc {
        state_mo[(i, i)] = 2.0;
    }
    for i in 0..nmo {
        for j in 0..nmo {
            state_mo[(i, j)] += particle_mo[(i, j)] - hole_mo[(i, j)];
        }
    }
    (
        mo_to_ao_density(c, &state_mo),
        mo_to_ao_density(c, &hole_mo),
        mo_to_ao_density(c, &particle_mo),
    )
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SpinChannel {
    Alpha,
    Beta,
}

#[derive(Clone, Copy, Debug)]
struct UcisSingle {
    spin: SpinChannel,
    occupied: usize,
    virtual_orbital: usize,
}

const STATE_SPECIFIC_ROOT_ISOLATION_EV: f64 = 1.0e-4;

fn ensure_state_specific_roots_are_isolated(
    energies: &[f64],
    requested_states: usize,
    sector: &str,
) -> Result<()> {
    let requested = energies
        .iter()
        .enumerate()
        .filter(|(_, energy)| energy.is_finite() && **energy > 0.0)
        .take(requested_states)
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    for &state in &requested {
        for (other, energy) in energies.iter().enumerate() {
            if state == other || !energy.is_finite() {
                continue;
            }
            let separation = (energies[state] - energy).abs();
            if separation < STATE_SPECIFIC_ROOT_ISOLATION_EV {
                return Err(XndoError::LinearAlgebra(format!(
                    "state-specific {sector} derivative is not unique: roots {} and {} are separated by {separation:.3e} eV; use a state-averaged or diabatic model",
                    state + 1,
                    other + 1,
                )));
            }
        }
    }
    Ok(())
}

fn ucis_singles(ground: &ZindoResult, options: &ZindoOptions) -> Vec<UcisSingle> {
    let nmo_alpha = ground.mo_coeff.cols;
    let nmo_beta = ground
        .mo_coeff_beta
        .as_ref()
        .map(|coefficients| coefficients.cols)
        .unwrap_or(0);
    let mut singles = Vec::new();
    for (spin, nocc, nmo) in [
        (SpinChannel::Alpha, ground.n_alpha, nmo_alpha),
        (SpinChannel::Beta, ground.n_beta, nmo_beta),
    ] {
        let occupied_start = options
            .active_occupied
            .map(|count| nocc.saturating_sub(count.min(nocc)))
            .unwrap_or(0);
        let virtual_end = options
            .active_virtual
            .map(|count| (nocc + count).min(nmo))
            .unwrap_or(nmo);
        for occupied in occupied_start..nocc {
            for virtual_orbital in nocc..virtual_end {
                singles.push(UcisSingle {
                    spin,
                    occupied,
                    virtual_orbital,
                });
            }
        }
    }
    singles
}

fn channel_coefficients(ground: &ZindoResult, spin: SpinChannel) -> &Matrix {
    match spin {
        SpinChannel::Alpha => &ground.mo_coeff,
        SpinChannel::Beta => ground
            .mo_coeff_beta
            .as_ref()
            .expect("UHF beta coefficients"),
    }
}

fn channel_energies(ground: &ZindoResult, spin: SpinChannel) -> &[f64] {
    match spin {
        SpinChannel::Alpha => &ground.mo_energies_ev,
        SpinChannel::Beta => ground
            .mo_energies_beta_ev
            .as_deref()
            .expect("UHF beta energies"),
    }
}

fn response_coefficients(response: &UhfResponse, spin: SpinChannel) -> &Matrix {
    match spin {
        SpinChannel::Alpha => &response.mo_coeff_alpha,
        SpinChannel::Beta => &response.mo_coeff_beta,
    }
}

fn response_energies(response: &UhfResponse, spin: SpinChannel) -> &[f64] {
    match spin {
        SpinChannel::Alpha => &response.mo_energies_alpha,
        SpinChannel::Beta => &response.mo_energies_beta,
    }
}

fn build_ucis_matrix(
    ground: &ZindoResult,
    singles: &[UcisSingle],
    j: &Matrix,
    k: &Matrix,
) -> Matrix {
    let mut matrix = Matrix::zeros(singles.len(), singles.len());
    for pidx in 0..singles.len() {
        let p = singles[pidx];
        let cp = channel_coefficients(ground, p.spin);
        for qidx in 0..=pidx {
            let q = singles[qidx];
            let cq = channel_coefficients(ground, q.spin);
            let direct = mo_eri_mixed(
                [cp, cp, cq, cq],
                j,
                k,
                [p.occupied, p.virtual_orbital, q.occupied, q.virtual_orbital],
            );
            let exchange = if p.spin == q.spin {
                mo_eri_mixed(
                    [cp, cp, cp, cp],
                    j,
                    k,
                    [p.occupied, q.occupied, p.virtual_orbital, q.virtual_orbital],
                )
            } else {
                0.0
            };
            let mut value = direct - exchange;
            if p.spin == q.spin
                && p.occupied == q.occupied
                && p.virtual_orbital == q.virtual_orbital
            {
                let energies = channel_energies(ground, p.spin);
                value += energies[p.virtual_orbital] - energies[p.occupied];
            }
            matrix[(pidx, qidx)] = value;
            matrix[(qidx, pidx)] = value;
        }
    }
    matrix
}

fn build_ucis_matrix_derivative(
    ground: &ZindoResult,
    response: &UhfResponse,
    singles: &[UcisSingle],
    j: &Matrix,
    jx: &Matrix,
    k: &Matrix,
    kx: &Matrix,
) -> Matrix {
    let mut matrix = Matrix::zeros(singles.len(), singles.len());
    for pidx in 0..singles.len() {
        let p = singles[pidx];
        let cp = channel_coefficients(ground, p.spin);
        let cpx = response_coefficients(response, p.spin);
        for qidx in 0..=pidx {
            let q = singles[qidx];
            let cq = channel_coefficients(ground, q.spin);
            let cqx = response_coefficients(response, q.spin);
            let direct = mo_eri_mixed_derivative(
                [cp, cp, cq, cq],
                [cpx, cpx, cqx, cqx],
                j,
                jx,
                k,
                kx,
                [p.occupied, p.virtual_orbital, q.occupied, q.virtual_orbital],
            );
            let exchange = if p.spin == q.spin {
                mo_eri_mixed_derivative(
                    [cp, cp, cp, cp],
                    [cpx, cpx, cpx, cpx],
                    j,
                    jx,
                    k,
                    kx,
                    [p.occupied, q.occupied, p.virtual_orbital, q.virtual_orbital],
                )
            } else {
                0.0
            };
            let mut value = direct - exchange;
            if p.spin == q.spin
                && p.occupied == q.occupied
                && p.virtual_orbital == q.virtual_orbital
            {
                let energies = response_energies(response, p.spin);
                value += energies[p.virtual_orbital] - energies[p.occupied];
            }
            matrix[(pidx, qidx)] = value;
            matrix[(qidx, pidx)] = value;
        }
    }
    matrix
}

fn second_response_coefficients(response: &UhfSecondResponse, spin: SpinChannel) -> &Matrix {
    match spin {
        SpinChannel::Alpha => &response.mo_coeff_alpha,
        SpinChannel::Beta => &response.mo_coeff_beta,
    }
}

fn second_response_energies(response: &UhfSecondResponse, spin: SpinChannel) -> &[f64] {
    match spin {
        SpinChannel::Alpha => &response.mo_energies_alpha,
        SpinChannel::Beta => &response.mo_energies_beta,
    }
}

#[allow(clippy::too_many_arguments)]
fn build_ucis_matrix_second_derivative(
    ground: &ZindoResult,
    response_x: &UhfResponse,
    response_y: &UhfResponse,
    response_xy: &UhfSecondResponse,
    singles: &[UcisSingle],
    j: &Matrix,
    jx: &Matrix,
    jy: &Matrix,
    jxy: &Matrix,
    k: &Matrix,
    kx: &Matrix,
    ky: &Matrix,
    kxy: &Matrix,
) -> Matrix {
    let mut matrix = Matrix::zeros(singles.len(), singles.len());
    for pidx in 0..singles.len() {
        let p = singles[pidx];
        let cp = channel_coefficients(ground, p.spin);
        let cpx = response_coefficients(response_x, p.spin);
        let cpy = response_coefficients(response_y, p.spin);
        let cpxy = second_response_coefficients(response_xy, p.spin);
        for qidx in 0..=pidx {
            let q = singles[qidx];
            let cq = channel_coefficients(ground, q.spin);
            let cqx = response_coefficients(response_x, q.spin);
            let cqy = response_coefficients(response_y, q.spin);
            let cqxy = second_response_coefficients(response_xy, q.spin);
            let direct = mo_eri_mixed_second_derivative(
                [cp, cp, cq, cq],
                [cpx, cpx, cqx, cqx],
                [cpy, cpy, cqy, cqy],
                [cpxy, cpxy, cqxy, cqxy],
                j,
                jx,
                jy,
                jxy,
                k,
                kx,
                ky,
                kxy,
                [p.occupied, p.virtual_orbital, q.occupied, q.virtual_orbital],
            );
            let exchange = if p.spin == q.spin {
                mo_eri_mixed_second_derivative(
                    [cp, cp, cp, cp],
                    [cpx, cpx, cpx, cpx],
                    [cpy, cpy, cpy, cpy],
                    [cpxy, cpxy, cpxy, cpxy],
                    j,
                    jx,
                    jy,
                    jxy,
                    k,
                    kx,
                    ky,
                    kxy,
                    [p.occupied, q.occupied, p.virtual_orbital, q.virtual_orbital],
                )
            } else {
                0.0
            };
            let mut value = direct - exchange;
            if p.spin == q.spin
                && p.occupied == q.occupied
                && p.virtual_orbital == q.virtual_orbital
            {
                let energies = second_response_energies(response_xy, p.spin);
                value += energies[p.virtual_orbital] - energies[p.occupied];
            }
            matrix[(pidx, qidx)] = value;
            matrix[(qidx, pidx)] = value;
        }
    }
    matrix
}

fn zindo_uhf_jk_derivative(data: &ZindoUhfResponseData, coordinate: usize) -> (Matrix, Matrix) {
    let mut jx = Matrix::zeros(data.basis.nao(), data.basis.nao());
    let mut kx = Matrix::zeros(data.basis.nao(), data.basis.nao());
    for pair in &data.pairs {
        let Some((axis, sign)) = pair_coordinate(pair, coordinate) else {
            continue;
        };
        let value = sign * pair.gamma.g[axis];
        let oa = data.basis.atom_offset[pair.ia];
        let ob = data.basis.atom_offset[pair.ib];
        for i in 0..data.basis.atom_norb[pair.ia] {
            for q in 0..data.basis.atom_norb[pair.ib] {
                let mu = oa + i;
                let nu = ob + q;
                jx[(mu, nu)] = value;
                jx[(nu, mu)] = value;
                kx[(mu, nu)] = value;
                kx[(nu, mu)] = value;
            }
        }
    }
    (jx, kx)
}

fn zindo_uhf_jk_second_derivative(
    data: &ZindoUhfResponseData,
    coordinate_x: usize,
    coordinate_y: usize,
) -> (Matrix, Matrix) {
    let mut jxy = Matrix::zeros(data.basis.nao(), data.basis.nao());
    let mut kxy = Matrix::zeros(data.basis.nao(), data.basis.nao());
    for pair in &data.pairs {
        let (Some((axis_x, sign_x)), Some((axis_y, sign_y))) = (
            pair_coordinate(pair, coordinate_x),
            pair_coordinate(pair, coordinate_y),
        ) else {
            continue;
        };
        let value = sign_x * sign_y * pair.gamma.h[axis_x][axis_y];
        let oa = data.basis.atom_offset[pair.ia];
        let ob = data.basis.atom_offset[pair.ib];
        for i in 0..data.basis.atom_norb[pair.ia] {
            for q in 0..data.basis.atom_norb[pair.ib] {
                let mu = oa + i;
                let nu = ob + q;
                jxy[(mu, nu)] = value;
                jxy[(nu, mu)] = value;
                kxy[(mu, nu)] = value;
                kxy[(nu, mu)] = value;
            }
        }
    }
    (jxy, kxy)
}

fn zindo_uhf_ground_gradient(
    mol: &Molecule,
    params: &ZindoParameters,
    ground: &ZindoResult,
    data: &ZindoUhfResponseData,
) -> Result<Vec<Vec3>> {
    let mut gradient = vec![Vec3::zero(); mol.atoms.len()];
    for pair in &data.pairs {
        let ea = params.element(mol.atoms[pair.ia].z)?;
        let eb = params.element(mol.atoms[pair.ib].z)?;
        let oa = data.basis.atom_offset[pair.ia];
        let ob = data.basis.atom_offset[pair.ib];
        let mut resonance = Dual2::constant(0.0);
        for i in 0..data.basis.atom_norb[pair.ia] {
            for q in 0..data.basis.atom_norb[pair.ib] {
                let beta = 0.5 * (ea.beta_for_orb(i) + eb.beta_for_orb(q));
                resonance =
                    resonance + pair.overlap[i][q] * (ground.density[(oa + i, ob + q)] * beta);
            }
        }
        let fixed = pair_energy(
            pair.gamma,
            pair.core,
            ea.core_charge,
            eb.core_charge,
            pair.population_a,
            pair.population_b,
            pair.exchange,
            resonance,
        );
        let value = Vec3::new(fixed.g[0], fixed.g[1], fixed.g[2]);
        gradient[pair.ia] -= value;
        gradient[pair.ib] += value;
    }
    Ok(gradient)
}

/// Spin-orbital CIS gradients on a ZINDO/S UHF reference.
pub fn zindo_s_ucis_gradients(
    mol: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
) -> Result<ZindoCisGradientResult> {
    let ground = run_zindo_s(mol, params, options)?;
    if !ground.unrestricted {
        return Err(XndoError::InvalidInput(
            "zindo_s_ucis_gradients requires reference=UHF".into(),
        ));
    }
    if mol.atoms.len() == 2 {
        return Err(XndoError::LinearAlgebra(
            "an open-shell diatomic state-specific UCIS gradient is not unique on a degenerate electronic surface; use a specified state-averaged or diabatic model"
                .into(),
        ));
    }
    let data = zindo_uhf_response_data(mol, params, &ground)?;
    let (j, k, _) = build_jk(mol, params, &data.basis)?;
    let singles = ucis_singles(&ground, options);
    let ground_gradient = zindo_uhf_ground_gradient(mol, params, &ground, &data)?;
    if singles.is_empty() {
        return Ok(ZindoCisGradientResult {
            ground,
            ground_gradient,
            spin: CisSpin::Unrestricted,
            states: Vec::new(),
        });
    }
    let matrix = build_ucis_matrix(&ground, &singles, &j, &k);
    let (omega, vectors) = symmetric_eigen(&matrix)?;
    ensure_state_specific_roots_are_isolated(&omega, options.n_states, "UCIS")?;
    let ncoord = 3 * mol.atoms.len();
    let mut derivatives = vec![vec![0.0; ncoord]; omega.len()];
    for coordinate in 0..ncoord {
        let (jx, kx) = zindo_uhf_jk_derivative(&data, coordinate);
        let first = build_ucis_matrix_derivative(
            &ground,
            &data.responses[coordinate],
            &singles,
            &j,
            &jx,
            &k,
            &kx,
        );
        for state in 0..omega.len() {
            for p in 0..singles.len() {
                for q in 0..singles.len() {
                    derivatives[state][coordinate] +=
                        vectors[(p, state)] * first[(p, q)] * vectors[(q, state)];
                }
            }
        }
    }
    let mut states = Vec::new();
    for state in 0..omega.len() {
        if !omega[state].is_finite() || omega[state] <= 0.0 {
            continue;
        }
        if states.len() >= options.n_states {
            break;
        }
        let excitation_gradient = (0..mol.atoms.len())
            .map(|atom| {
                Vec3::new(
                    derivatives[state][3 * atom],
                    derivatives[state][3 * atom + 1],
                    derivatives[state][3 * atom + 2],
                )
            })
            .collect::<Vec<_>>();
        let state_gradient = ground_gradient
            .iter()
            .zip(&excitation_gradient)
            .map(|(ground_value, excitation)| *ground_value + *excitation)
            .collect();
        states.push(ExcitedStateGradient {
            root: state + 1,
            spin: CisSpin::Unrestricted,
            excitation_energy_ev: omega[state],
            state_total_energy_ev: ground.total_ev + omega[state],
            excitation_gradient,
            state_gradient,
        });
    }
    Ok(ZindoCisGradientResult {
        ground,
        ground_gradient,
        spin: CisSpin::Unrestricted,
        states,
    })
}

/// Fully analytic state Hessians for spin-orbital CIS on a ZINDO/S UHF
/// reference. Coupled first- and second-order alpha/beta orbital response is
/// evaluated without displaced SCF or CIS calculations.
pub fn zindo_s_ucis_hessians(
    mol: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
) -> Result<ZindoCisHessianResult> {
    let (ground, ground_gradient, ground_hessian) = analytic_ground_hessian(mol, params, options)?;
    if !ground.unrestricted {
        return Err(XndoError::InvalidInput(
            "zindo_s_ucis_hessians requires reference=UHF".into(),
        ));
    }
    let data = zindo_uhf_response_data(mol, params, &ground)?;
    let basis = &data.basis;
    let (j, k, _) = build_jk(mol, params, basis)?;
    let singles = ucis_singles(&ground, options);
    if singles.is_empty() {
        return Ok(ZindoCisHessianResult {
            ground,
            ground_gradient,
            ground_hessian,
            spin: CisSpin::Unrestricted,
            states: Vec::new(),
        });
    }
    let matrix = build_ucis_matrix(&ground, &singles, &j, &k);
    let (omega, vectors) = symmetric_eigen(&matrix)?;
    ensure_state_specific_roots_are_isolated(&omega, options.n_states, "UCIS")?;
    let ncoord = 3 * mol.atoms.len();
    let jk_first = (0..ncoord)
        .map(|coordinate| zindo_uhf_jk_derivative(&data, coordinate))
        .collect::<Vec<_>>();
    let a_first = (0..ncoord)
        .map(|coordinate| {
            let (jx, kx) = &jk_first[coordinate];
            build_ucis_matrix_derivative(
                &ground,
                &data.responses[coordinate],
                &singles,
                &j,
                jx,
                &k,
                kx,
            )
        })
        .collect::<Vec<_>>();
    let coefficients_beta = ground
        .mo_coeff_beta
        .as_ref()
        .expect("UHF beta coefficients");
    let energies_beta = ground
        .mo_energies_beta_ev
        .as_ref()
        .expect("UHF beta energies");
    let zero = Matrix::zeros(basis.nao(), basis.nao());
    let mut excitation_hessians = (0..omega.len())
        .map(|_| Matrix::zeros(ncoord, ncoord))
        .collect::<Vec<_>>();
    for coordinate_x in 0..ncoord {
        for coordinate_y in 0..=coordinate_x {
            let (known_alpha, known_beta) =
                zindo_known_fock_second_uhf(mol, params, &data, coordinate_x, coordinate_y)?;
            let second = solve_uhf_second_response(
                &ground.mo_coeff,
                coefficients_beta,
                &ground.mo_energies_ev,
                energies_beta,
                ground.n_alpha,
                ground.n_beta,
                &data.responses[coordinate_x],
                &data.responses[coordinate_y],
                &known_alpha,
                &known_beta,
                |dpa, dpb| Ok(build_fock_uhf(&zero, dpa, dpb, &j, &k, basis)),
            )?;
            let (jxy, kxy) = zindo_uhf_jk_second_derivative(&data, coordinate_x, coordinate_y);
            let (jx, kx) = &jk_first[coordinate_x];
            let (jy, ky) = &jk_first[coordinate_y];
            let a_second = build_ucis_matrix_second_derivative(
                &ground,
                &data.responses[coordinate_x],
                &data.responses[coordinate_y],
                &second,
                &singles,
                &j,
                jx,
                jy,
                &jxy,
                &k,
                kx,
                ky,
                &kxy,
            );
            for state in 0..omega.len() {
                let mut value = 0.0;
                for p in 0..singles.len() {
                    for q in 0..singles.len() {
                        value += vectors[(p, state)] * a_second[(p, q)] * vectors[(q, state)];
                    }
                }
                for other in 0..omega.len() {
                    if other == state {
                        continue;
                    }
                    let gap = omega[state] - omega[other];
                    if gap.abs() < 1.0e-9 {
                        continue;
                    }
                    let mut coupling_x = 0.0;
                    let mut coupling_y = 0.0;
                    for p in 0..singles.len() {
                        for q in 0..singles.len() {
                            coupling_x += vectors[(p, state)]
                                * a_first[coordinate_x][(p, q)]
                                * vectors[(q, other)];
                            coupling_y += vectors[(p, state)]
                                * a_first[coordinate_y][(p, q)]
                                * vectors[(q, other)];
                        }
                    }
                    value += 2.0 * coupling_x * coupling_y / gap;
                }
                excitation_hessians[state][(coordinate_x, coordinate_y)] = value;
                excitation_hessians[state][(coordinate_y, coordinate_x)] = value;
            }
        }
    }
    let mut states = Vec::new();
    for state in 0..omega.len() {
        if !omega[state].is_finite() || omega[state] <= 0.0 {
            continue;
        }
        if states.len() >= options.n_states {
            break;
        }
        let excitation_gradient = (0..mol.atoms.len())
            .map(|atom| {
                let component = |coordinate: usize| {
                    let mut value = 0.0;
                    for p in 0..singles.len() {
                        for q in 0..singles.len() {
                            value += vectors[(p, state)]
                                * a_first[coordinate][(p, q)]
                                * vectors[(q, state)];
                        }
                    }
                    value
                };
                Vec3::new(
                    component(3 * atom),
                    component(3 * atom + 1),
                    component(3 * atom + 2),
                )
            })
            .collect::<Vec<_>>();
        let state_gradient = ground_gradient
            .iter()
            .zip(&excitation_gradient)
            .map(|(ground_value, excitation)| *ground_value + *excitation)
            .collect();
        let mut state_hessian = ground_hessian.clone();
        for index in 0..state_hessian.as_slice().len() {
            state_hessian.as_mut_slice()[index] += excitation_hessians[state].as_slice()[index];
        }
        states.push(ExcitedStateHessian {
            root: state + 1,
            spin: CisSpin::Unrestricted,
            excitation_energy_ev: omega[state],
            state_total_energy_ev: ground.total_ev + omega[state],
            excitation_gradient,
            state_gradient,
            excitation_hessian: excitation_hessians[state].clone(),
            state_hessian,
        });
    }
    Ok(ZindoCisHessianResult {
        ground,
        ground_gradient,
        ground_hessian,
        spin: CisSpin::Unrestricted,
        states,
    })
}

fn ucis_state_density_and_attachment(
    ground: &ZindoResult,
    singles: &[UcisSingle],
    vectors: &Matrix,
    state: usize,
) -> (Matrix, Matrix, Matrix) {
    let nmo_alpha = ground.mo_coeff.cols;
    let nmo_beta = ground
        .mo_coeff_beta
        .as_ref()
        .expect("UHF beta coefficients")
        .cols;
    let mut hole_alpha = Matrix::zeros(nmo_alpha, nmo_alpha);
    let mut particle_alpha = Matrix::zeros(nmo_alpha, nmo_alpha);
    let mut hole_beta = Matrix::zeros(nmo_beta, nmo_beta);
    let mut particle_beta = Matrix::zeros(nmo_beta, nmo_beta);
    for pidx in 0..singles.len() {
        let p = singles[pidx];
        for qidx in 0..singles.len() {
            let q = singles[qidx];
            if p.spin != q.spin {
                continue;
            }
            let weight = vectors[(pidx, state)] * vectors[(qidx, state)];
            let (hole, particle) = match p.spin {
                SpinChannel::Alpha => (&mut hole_alpha, &mut particle_alpha),
                SpinChannel::Beta => (&mut hole_beta, &mut particle_beta),
            };
            if p.virtual_orbital == q.virtual_orbital {
                hole[(p.occupied, q.occupied)] += weight;
            }
            if p.occupied == q.occupied {
                particle[(p.virtual_orbital, q.virtual_orbital)] += weight;
            }
        }
    }
    let cb = ground
        .mo_coeff_beta
        .as_ref()
        .expect("UHF beta coefficients");
    let hole = add_matrix(
        &mo_to_ao_density(&ground.mo_coeff, &hole_alpha),
        &mo_to_ao_density(cb, &hole_beta),
    );
    let particle = add_matrix(
        &mo_to_ao_density(&ground.mo_coeff, &particle_alpha),
        &mo_to_ao_density(cb, &particle_beta),
    );
    let mut state_density = ground.density.clone();
    for index in 0..state_density.as_slice().len() {
        state_density.as_mut_slice()[index] += particle.as_slice()[index] - hole.as_slice()[index];
    }
    (state_density, hole, particle)
}

/// Spin-conserving spin-orbital CIS spectrum on a ZINDO/S UHF reference.
pub fn zindo_s_ucis(
    mol: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
) -> Result<ZindoSpectrum> {
    let ground = run_zindo_s(mol, params, options)?;
    if !ground.unrestricted {
        return Err(XndoError::InvalidInput(
            "zindo_s_ucis requires reference=UHF".into(),
        ));
    }
    let basis = ZBasis::build(mol, params)?;
    let (j, k, _) = build_jk(mol, params, &basis)?;
    let r_ao = ao_position_matrices(mol, params, &basis)?;
    let ground_dipole_au = permanent_dipole_au(mol, params, &basis, &ground.density, &r_ao)?;
    let ground_dipole_debye = scale3(ground_dipole_au, DEBYE_PER_E_BOHR);
    let singles = ucis_singles(&ground, options);
    if singles.is_empty() {
        return Ok(ZindoSpectrum {
            ground,
            spin: CisSpin::Unrestricted,
            ground_dipole_au,
            ground_dipole_debye,
            states: Vec::new(),
        });
    }
    let matrix = build_ucis_matrix(&ground, &singles, &j, &k);
    let (omega, vectors) = symmetric_eigen(&matrix)?;
    let mut transition_integrals = vec![[0.0; 3]; singles.len()];
    for (index, single) in singles.iter().enumerate() {
        let coefficients = channel_coefficients(&ground, single.spin);
        for axis in 0..3 {
            transition_integrals[index][axis] = mo_position(
                coefficients,
                &r_ao[axis],
                single.occupied,
                single.virtual_orbital,
            );
        }
    }
    let mut states = Vec::new();
    for state in 0..omega.len() {
        let energy = omega[state];
        if !energy.is_finite() || energy <= 0.0 {
            continue;
        }
        if states.len() >= options.n_states {
            break;
        }
        let mut transition = [0.0; 3];
        let mut coefficients = Vec::with_capacity(singles.len());
        for index in 0..singles.len() {
            let coefficient = vectors[(index, state)];
            for axis in 0..3 {
                transition[axis] += coefficient * transition_integrals[index][axis];
            }
            coefficients.push((coefficient.abs(), index, coefficient));
        }
        coefficients.sort_by(|left, right| {
            right
                .0
                .partial_cmp(&left.0)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let dominant = coefficients
            .into_iter()
            .take(5)
            .map(|(_, index, coefficient)| CiContribution {
                occupied: singles[index].occupied + 1,
                virtual_orbital: singles[index].virtual_orbital + 1,
                coefficient,
            })
            .collect();
        let dipole_strength = transition.iter().map(|value| value * value).sum::<f64>();
        let oscillator_strength = ((2.0 / 3.0) * (energy / ZINDO_AU2EV) * dipole_strength).max(0.0);
        let (state_density, hole_density, particle_density) =
            ucis_state_density_and_attachment(&ground, &singles, &vectors, state);
        let state_charges = charges(mol, params, &basis, &state_density)?;
        let hole_population = atom_diagonal_population(&basis, &hole_density);
        let electron_population = atom_diagonal_population(&basis, &particle_density);
        let hole_centroid_angstrom = population_centroid_angstrom(mol, &hole_population);
        let electron_centroid_angstrom = population_centroid_angstrom(mol, &electron_population);
        let delta = [
            electron_centroid_angstrom[0] - hole_centroid_angstrom[0],
            electron_centroid_angstrom[1] - hole_centroid_angstrom[1],
            electron_centroid_angstrom[2] - hole_centroid_angstrom[2],
        ];
        let state_dipole_au = permanent_dipole_au(mol, params, &basis, &state_density, &r_ao)?;
        let state_dipole_debye = scale3(state_dipole_au, DEBYE_PER_E_BOHR);
        let difference_dipole_au = [
            state_dipole_au[0] - ground_dipole_au[0],
            state_dipole_au[1] - ground_dipole_au[1],
            state_dipole_au[2] - ground_dipole_au[2],
        ];
        let difference_dipole_debye = scale3(difference_dipole_au, DEBYE_PER_E_BOHR);
        states.push(ExcitedState {
            spin: CisSpin::Unrestricted,
            spin_multiplicity: 0,
            s2_expectation: f64::NAN,
            energy_ev: energy,
            state_total_energy_ev: ground.total_ev + energy,
            state_total_energy_hartree: (ground.total_ev + energy) / ZINDO_AU2EV,
            energy_cm1: energy * EV_TO_WAVENUMBER_CM1,
            wavelength_nm: Some(EV_NM / energy),
            oscillator_strength,
            transition_dipole_au: transition,
            transition_dipole_debye: scale3(transition, DEBYE_PER_E_BOHR),
            transition_dipole_magnitude_au: norm3(transition),
            dipole_strength_au2: dipole_strength,
            permanent_dipole_au: state_dipole_au,
            permanent_dipole_debye: state_dipole_debye,
            permanent_dipole_magnitude_debye: norm3(state_dipole_debye),
            difference_dipole_au,
            difference_dipole_debye,
            difference_dipole_magnitude_debye: norm3(difference_dipole_debye),
            charges: state_charges,
            hole_population,
            electron_population,
            hole_centroid_angstrom,
            electron_centroid_angstrom,
            charge_transfer_distance_angstrom: norm3(delta),
            dominant,
        });
    }
    Ok(ZindoSpectrum {
        ground,
        spin: CisSpin::Unrestricted,
        ground_dipole_au,
        ground_dipole_debye,
        states,
    })
}

/// Analytic gradients of spin-adapted ZINDO/S CIS roots.
///
/// The ground-state part is the variational RHF gradient. Excitation-energy
/// derivatives include AO integral derivatives and the full static CPHF
/// orbital response; no displaced-geometry SCF or CIS calculation is used.
pub fn zindo_s_cis_gradients(
    mol: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
    spin: CisSpin,
) -> Result<ZindoCisGradientResult> {
    let ground = run_zindo_s(mol, params, options)?;
    let response_data = zindo_response_data(mol, params, &ground)?;
    let basis = &response_data.basis;
    let (j, k, _) = build_jk(mol, params, basis)?;
    let nmo = ground.mo_energies_ev.len();
    let nocc = ground.n_occ;
    let occ_start = options
        .active_occupied
        .map(|n| nocc.saturating_sub(n.min(nocc)))
        .unwrap_or(0);
    let virt_end = options
        .active_virtual
        .map(|n| (nocc + n).min(nmo))
        .unwrap_or(nmo);
    let singles = (occ_start..nocc)
        .flat_map(|i| (nocc..virt_end).map(move |a| (i, a)))
        .collect::<Vec<_>>();
    if singles.is_empty() {
        return Ok(ZindoCisGradientResult {
            ground,
            ground_gradient: vec![Vec3::zero(); mol.atoms.len()],
            spin,
            states: Vec::new(),
        });
    }
    let ns = singles.len();
    let mut a_matrix = Matrix::zeros(ns, ns);
    for pidx in 0..ns {
        let (i, a) = singles[pidx];
        for qidx in 0..=pidx {
            let (jmo, b) = singles[qidx];
            let direct = mo_eri(&ground.mo_coeff, &j, &k, i, a, jmo, b);
            let exchange = mo_eri(&ground.mo_coeff, &j, &k, i, jmo, a, b);
            let mut value = match spin {
                CisSpin::Singlet => 2.0 * direct - exchange,
                CisSpin::Triplet => -exchange,
                CisSpin::Unrestricted => unreachable!("UCIS uses the spin-orbital builder"),
            };
            if i == jmo && a == b {
                value += ground.mo_energies_ev[a] - ground.mo_energies_ev[i];
            }
            a_matrix[(pidx, qidx)] = value;
            a_matrix[(qidx, pidx)] = value;
        }
    }
    let (omega, x) = symmetric_eigen(&a_matrix)?;
    ensure_state_specific_roots_are_isolated(&omega, options.n_states, spin.as_str())?;
    let ncoord = 3 * mol.atoms.len();
    let mut omega_derivatives = vec![vec![0.0; ncoord]; omega.len()];
    for (coord, response) in response_data.responses.iter().enumerate() {
        let mut jx = Matrix::zeros(basis.nao(), basis.nao());
        let mut kx = Matrix::zeros(basis.nao(), basis.nao());
        for pair in &response_data.pairs {
            let axis = coord % 3;
            let atom = coord / 3;
            let sign = if atom == pair.ia {
                -1.0
            } else if atom == pair.ib {
                1.0
            } else {
                0.0
            };
            if sign == 0.0 {
                continue;
            }
            let value = sign * pair.gamma.g[axis];
            let oa = basis.atom_offset[pair.ia];
            let ob = basis.atom_offset[pair.ib];
            for i in 0..basis.atom_norb[pair.ia] {
                for q in 0..basis.atom_norb[pair.ib] {
                    let mu = oa + i;
                    let nu = ob + q;
                    jx[(mu, nu)] = value;
                    jx[(nu, mu)] = value;
                    kx[(mu, nu)] = value;
                    kx[(nu, mu)] = value;
                }
            }
        }
        let mut ax = Matrix::zeros(ns, ns);
        for pidx in 0..ns {
            let (i, a) = singles[pidx];
            for qidx in 0..=pidx {
                let (jmo, b) = singles[qidx];
                let direct = mo_eri_derivative(
                    &ground.mo_coeff,
                    &response.mo_coeff,
                    &j,
                    &jx,
                    &k,
                    &kx,
                    i,
                    a,
                    jmo,
                    b,
                );
                let exchange = mo_eri_derivative(
                    &ground.mo_coeff,
                    &response.mo_coeff,
                    &j,
                    &jx,
                    &k,
                    &kx,
                    i,
                    jmo,
                    a,
                    b,
                );
                let mut value = match spin {
                    CisSpin::Singlet => 2.0 * direct - exchange,
                    CisSpin::Triplet => -exchange,
                    CisSpin::Unrestricted => unreachable!("UCIS uses the spin-orbital builder"),
                };
                if i == jmo && a == b {
                    value += response.mo_energies[a] - response.mo_energies[i];
                }
                ax[(pidx, qidx)] = value;
                ax[(qidx, pidx)] = value;
            }
        }
        for state in 0..omega.len() {
            let mut derivative = 0.0;
            for p in 0..ns {
                for q in 0..ns {
                    derivative += x[(p, state)] * ax[(p, q)] * x[(q, state)];
                }
            }
            omega_derivatives[state][coord] = derivative;
        }
    }

    let mut ground_gradient = vec![Vec3::zero(); mol.atoms.len()];
    for pair in &response_data.pairs {
        let ea = params.element(mol.atoms[pair.ia].z)?;
        let eb = params.element(mol.atoms[pair.ib].z)?;
        let oa = basis.atom_offset[pair.ia];
        let ob = basis.atom_offset[pair.ib];
        let mut resonance = Dual2::constant(0.0);
        for i in 0..basis.atom_norb[pair.ia] {
            for q in 0..basis.atom_norb[pair.ib] {
                let beta = 0.5 * (ea.beta_for_orb(i) + eb.beta_for_orb(q));
                resonance =
                    resonance + pair.overlap[i][q] * (ground.density[(oa + i, ob + q)] * beta);
            }
        }
        let fixed = pair_energy(
            pair.gamma,
            pair.core,
            ea.core_charge,
            eb.core_charge,
            pair.population_a,
            pair.population_b,
            pair.exchange,
            resonance,
        );
        let value = Vec3::new(fixed.g[0], fixed.g[1], fixed.g[2]);
        ground_gradient[pair.ia] -= value;
        ground_gradient[pair.ib] += value;
    }
    let mut states = Vec::new();
    for state in 0..omega.len() {
        if !omega[state].is_finite() || omega[state] <= 0.0 {
            continue;
        }
        if states.len() >= options.n_states {
            break;
        }
        let excitation_gradient = (0..mol.atoms.len())
            .map(|atom| {
                Vec3::new(
                    omega_derivatives[state][3 * atom],
                    omega_derivatives[state][3 * atom + 1],
                    omega_derivatives[state][3 * atom + 2],
                )
            })
            .collect::<Vec<_>>();
        let state_gradient = ground_gradient
            .iter()
            .zip(&excitation_gradient)
            .map(|(g, xg)| *g + *xg)
            .collect();
        states.push(ExcitedStateGradient {
            root: state + 1,
            spin,
            excitation_energy_ev: omega[state],
            state_total_energy_ev: ground.total_ev + omega[state],
            excitation_gradient,
            state_gradient,
        });
    }
    Ok(ZindoCisGradientResult {
        ground,
        ground_gradient,
        spin,
        states,
    })
}

/// Fully analytic Hessians of spin-adapted ZINDO/S CIS roots.
///
/// Mixed second-order CPHF response and the ordinary symmetric-eigenvalue
/// second-derivative expression are evaluated directly. Finite differences
/// are used only by the regression tests.
pub fn zindo_s_cis_hessians(
    mol: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
    spin: CisSpin,
) -> Result<ZindoCisHessianResult> {
    let (ground, ground_gradient, ground_hessian) = analytic_ground_hessian(mol, params, options)?;
    let data = zindo_response_data(mol, params, &ground)?;
    let basis = &data.basis;
    let (j, k, _) = build_jk(mol, params, basis)?;
    let nmo = ground.mo_energies_ev.len();
    let nocc = ground.n_occ;
    let occ_start = options
        .active_occupied
        .map(|n| nocc.saturating_sub(n.min(nocc)))
        .unwrap_or(0);
    let virt_end = options
        .active_virtual
        .map(|n| (nocc + n).min(nmo))
        .unwrap_or(nmo);
    let singles = (occ_start..nocc)
        .flat_map(|i| (nocc..virt_end).map(move |a| (i, a)))
        .collect::<Vec<_>>();
    if singles.is_empty() {
        return Ok(ZindoCisHessianResult {
            ground,
            ground_gradient,
            ground_hessian,
            spin,
            states: Vec::new(),
        });
    }
    let ns = singles.len();
    let ncoord = 3 * mol.atoms.len();
    let mut a_matrix = Matrix::zeros(ns, ns);
    for pidx in 0..ns {
        let (i, a) = singles[pidx];
        for qidx in 0..=pidx {
            let (jmo, b) = singles[qidx];
            let direct = mo_eri(&ground.mo_coeff, &j, &k, i, a, jmo, b);
            let exchange = mo_eri(&ground.mo_coeff, &j, &k, i, jmo, a, b);
            let mut value = match spin {
                CisSpin::Singlet => 2.0 * direct - exchange,
                CisSpin::Triplet => -exchange,
                CisSpin::Unrestricted => unreachable!("UCIS uses the spin-orbital builder"),
            };
            if i == jmo && a == b {
                value += ground.mo_energies_ev[a] - ground.mo_energies_ev[i];
            }
            a_matrix[(pidx, qidx)] = value;
            a_matrix[(qidx, pidx)] = value;
        }
    }
    let (omega, x) = symmetric_eigen(&a_matrix)?;
    ensure_state_specific_roots_are_isolated(&omega, options.n_states, spin.as_str())?;
    let jk_first = (0..ncoord)
        .map(|coordinate| zindo_jk_derivative(&data, coordinate))
        .collect::<Vec<_>>();
    let mut a_first = Vec::with_capacity(ncoord);
    for coordinate in 0..ncoord {
        let (jx, kx) = &jk_first[coordinate];
        let response = &data.responses[coordinate];
        let mut derivative = Matrix::zeros(ns, ns);
        for pidx in 0..ns {
            let (i, a) = singles[pidx];
            for qidx in 0..=pidx {
                let (jmo, b) = singles[qidx];
                let direct = mo_eri_derivative(
                    &ground.mo_coeff,
                    &response.mo_coeff,
                    &j,
                    jx,
                    &k,
                    kx,
                    i,
                    a,
                    jmo,
                    b,
                );
                let exchange = mo_eri_derivative(
                    &ground.mo_coeff,
                    &response.mo_coeff,
                    &j,
                    jx,
                    &k,
                    kx,
                    i,
                    jmo,
                    a,
                    b,
                );
                let mut value = match spin {
                    CisSpin::Singlet => 2.0 * direct - exchange,
                    CisSpin::Triplet => -exchange,
                    CisSpin::Unrestricted => unreachable!("UCIS uses the spin-orbital builder"),
                };
                if i == jmo && a == b {
                    value += response.mo_energies[a] - response.mo_energies[i];
                }
                derivative[(pidx, qidx)] = value;
                derivative[(qidx, pidx)] = value;
            }
        }
        a_first.push(derivative);
    }
    let zero = Matrix::zeros(basis.nao(), basis.nao());
    let mut excitation_hessians = (0..omega.len())
        .map(|_| Matrix::zeros(ncoord, ncoord))
        .collect::<Vec<_>>();
    for coordinate_x in 0..ncoord {
        for coordinate_y in 0..=coordinate_x {
            let known =
                zindo_known_fock_second(mol, params, &ground, &data, coordinate_x, coordinate_y)?;
            let second = solve_rhf_second_response(
                &ground.mo_coeff,
                &ground.mo_energies_ev,
                ground.n_occ,
                &data.responses[coordinate_x],
                &data.responses[coordinate_y],
                &known,
                |density| Ok(build_fock_closed(&zero, density, &j, &k, basis)),
            )?;
            let (jxy, kxy) = zindo_jk_second_derivative(&data, coordinate_x, coordinate_y);
            let (jx, kx) = &jk_first[coordinate_x];
            let (jy, ky) = &jk_first[coordinate_y];
            let mut a_second = Matrix::zeros(ns, ns);
            for pidx in 0..ns {
                let (i, a) = singles[pidx];
                for qidx in 0..=pidx {
                    let (jmo, b) = singles[qidx];
                    let direct = mo_eri_second_derivative(
                        &ground.mo_coeff,
                        &data.responses[coordinate_x].mo_coeff,
                        &data.responses[coordinate_y].mo_coeff,
                        &second.mo_coeff,
                        &j,
                        jx,
                        jy,
                        &jxy,
                        &k,
                        kx,
                        ky,
                        &kxy,
                        i,
                        a,
                        jmo,
                        b,
                    );
                    let exchange = mo_eri_second_derivative(
                        &ground.mo_coeff,
                        &data.responses[coordinate_x].mo_coeff,
                        &data.responses[coordinate_y].mo_coeff,
                        &second.mo_coeff,
                        &j,
                        jx,
                        jy,
                        &jxy,
                        &k,
                        kx,
                        ky,
                        &kxy,
                        i,
                        jmo,
                        a,
                        b,
                    );
                    let mut value = match spin {
                        CisSpin::Singlet => 2.0 * direct - exchange,
                        CisSpin::Triplet => -exchange,
                        CisSpin::Unrestricted => unreachable!("UCIS uses the spin-orbital builder"),
                    };
                    if i == jmo && a == b {
                        value += second.mo_energies[a] - second.mo_energies[i];
                    }
                    a_second[(pidx, qidx)] = value;
                    a_second[(qidx, pidx)] = value;
                }
            }
            for state in 0..omega.len() {
                let mut value = 0.0;
                for p in 0..ns {
                    for q in 0..ns {
                        value += x[(p, state)] * a_second[(p, q)] * x[(q, state)];
                    }
                }
                for other in 0..omega.len() {
                    if other == state {
                        continue;
                    }
                    let gap = omega[state] - omega[other];
                    if gap.abs() < 1.0e-9 {
                        continue;
                    }
                    let mut x_coupling = 0.0;
                    let mut y_coupling = 0.0;
                    for p in 0..ns {
                        for q in 0..ns {
                            x_coupling +=
                                x[(p, state)] * a_first[coordinate_x][(p, q)] * x[(q, other)];
                            y_coupling +=
                                x[(p, state)] * a_first[coordinate_y][(p, q)] * x[(q, other)];
                        }
                    }
                    value += 2.0 * x_coupling * y_coupling / gap;
                }
                excitation_hessians[state][(coordinate_x, coordinate_y)] = value;
                excitation_hessians[state][(coordinate_y, coordinate_x)] = value;
            }
        }
    }
    let mut states = Vec::new();
    for state in 0..omega.len() {
        if !omega[state].is_finite() || omega[state] <= 0.0 {
            continue;
        }
        if states.len() >= options.n_states {
            break;
        }
        let excitation_gradient = (0..mol.atoms.len())
            .map(|atom| {
                let component = |coordinate: usize| {
                    let mut value = 0.0;
                    for p in 0..ns {
                        for q in 0..ns {
                            value += x[(p, state)] * a_first[coordinate][(p, q)] * x[(q, state)];
                        }
                    }
                    value
                };
                Vec3::new(
                    component(3 * atom),
                    component(3 * atom + 1),
                    component(3 * atom + 2),
                )
            })
            .collect::<Vec<_>>();
        let state_gradient = ground_gradient
            .iter()
            .zip(&excitation_gradient)
            .map(|(ground_value, excitation)| *ground_value + *excitation)
            .collect();
        let mut state_hessian = ground_hessian.clone();
        for index in 0..state_hessian.as_slice().len() {
            state_hessian.as_mut_slice()[index] += excitation_hessians[state].as_slice()[index];
        }
        states.push(ExcitedStateHessian {
            root: state + 1,
            spin,
            excitation_energy_ev: omega[state],
            state_total_energy_ev: ground.total_ev + omega[state],
            excitation_gradient,
            state_gradient,
            excitation_hessian: excitation_hessians[state].clone(),
            state_hessian,
        });
    }
    Ok(ZindoCisHessianResult {
        ground,
        ground_gradient,
        ground_hessian,
        spin,
        states,
    })
}

/// Backward-compatible singlet ZINDO/S CIS entry point.
pub fn zindo_s_cis(
    mol: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
) -> Result<ZindoSpectrum> {
    zindo_s_cis_spin(mol, params, options, CisSpin::Singlet)
}

/// Vertical electric-dipole UV-visible stick spectrum.
///
/// This is the explicitly named spectroscopic entry point for the singlet
/// ZINDO/S CIS sector. Each state includes excitation energy, wavelength,
/// oscillator and dipole strengths, transition/permanent/difference dipoles,
/// state charges, dominant configurations, and hole/electron descriptors.
/// Line broadening and solvent shifts are deliberately not imposed.
pub fn uv_vis_spectrum(
    mol: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
) -> Result<ZindoSpectrum> {
    zindo_s_cis_spin(mol, params, options, CisSpin::Singlet)
}

/// Vertical spin-adapted CIS spectrum on top of the converged ZINDO/S RHF state.
///
/// Singlet matrix: `Delta_e + 2(ia|jb) - (ij|ab)`.
/// Triplet matrix: `Delta_e - (ij|ab)`.  Electric-dipole transition moments
/// from the singlet RHF ground state are zeroed for triplets in this spin-free
/// implementation (no spin-orbit intensity borrowing).
pub fn zindo_s_cis_spin(
    mol: &Molecule,
    params: &ZindoParameters,
    options: &ZindoOptions,
    spin: CisSpin,
) -> Result<ZindoSpectrum> {
    let ground = run_zindo_s(mol, params, options)?;
    let basis = ZBasis::build(mol, params)?;
    let (j, k, _) = build_jk(mol, params, &basis)?;
    let r_ao = ao_position_matrices(mol, params, &basis)?;
    let ground_dipole_au = permanent_dipole_au(mol, params, &basis, &ground.density, &r_ao)?;
    let ground_dipole_debye = scale3(ground_dipole_au, DEBYE_PER_E_BOHR);
    let nmo = ground.mo_energies_ev.len();
    let nocc = ground.n_occ;
    if nocc == 0 || nocc >= nmo {
        return Ok(ZindoSpectrum {
            ground,
            spin,
            ground_dipole_au,
            ground_dipole_debye,
            states: Vec::new(),
        });
    }
    let occ_start = options
        .active_occupied
        .map(|n| nocc.saturating_sub(n.min(nocc)))
        .unwrap_or(0);
    let virt_end = options
        .active_virtual
        .map(|n| (nocc + n).min(nmo))
        .unwrap_or(nmo);
    let singles: Vec<(usize, usize)> = (occ_start..nocc)
        .flat_map(|i| (nocc..virt_end).map(move |a| (i, a)))
        .collect();
    if singles.is_empty() {
        return Ok(ZindoSpectrum {
            ground,
            spin,
            ground_dipole_au,
            ground_dipole_debye,
            states: Vec::new(),
        });
    }
    let ns = singles.len();
    let mut a_mat = Matrix::zeros(ns, ns);
    for pidx in 0..ns {
        let (i, a) = singles[pidx];
        for qidx in 0..=pidx {
            let (jmo, b) = singles[qidx];
            let direct = mo_eri(&ground.mo_coeff, &j, &k, i, a, jmo, b);
            let exchange = mo_eri(&ground.mo_coeff, &j, &k, i, jmo, a, b);
            let mut v = match spin {
                CisSpin::Singlet => 2.0 * direct - exchange,
                CisSpin::Triplet => -exchange,
                CisSpin::Unrestricted => unreachable!("UCIS uses the spin-orbital builder"),
            };
            if i == jmo && a == b {
                v += ground.mo_energies_ev[a] - ground.mo_energies_ev[i];
            }
            a_mat[(pidx, qidx)] = v;
            a_mat[(qidx, pidx)] = v;
        }
    }
    let (omega, x) = symmetric_eigen(&a_mat)?;
    let mut ria = vec![[0.0; 3]; ns];
    if spin == CisSpin::Singlet {
        for (idx, &(i, a)) in singles.iter().enumerate() {
            for axis in 0..3 {
                ria[idx][axis] = mo_position(&ground.mo_coeff, &r_ao[axis], i, a);
            }
        }
    }
    let mut states = Vec::new();
    for st in 0..omega.len().min(options.n_states) {
        let e = omega[st];
        if !e.is_finite() || e <= 0.0 {
            continue;
        }
        let mut mu = [0.0; 3];
        let mut coeffs = Vec::with_capacity(ns);
        for idx in 0..ns {
            let coeff = x[(idx, st)];
            if spin == CisSpin::Singlet {
                for axis in 0..3 {
                    mu[axis] += 2.0_f64.sqrt() * coeff * ria[idx][axis];
                }
            }
            coeffs.push((coeff.abs(), idx, coeff));
        }
        coeffs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
        let dominant = coeffs
            .into_iter()
            .take(5)
            .map(|(_, idx, coefficient)| {
                let (i, a) = singles[idx];
                CiContribution {
                    occupied: i + 1,
                    virtual_orbital: a + 1,
                    coefficient,
                }
            })
            .collect();
        let mu2 = mu.iter().map(|v| v * v).sum::<f64>();
        let fosc = if spin == CisSpin::Singlet {
            ((2.0 / 3.0) * (e / ZINDO_AU2EV) * mu2).max(0.0)
        } else {
            0.0
        };
        let (state_density, hole_density, particle_density) =
            state_density_and_attachment(&ground.mo_coeff, &singles, &x, st, nocc);
        let state_charges = charges(mol, params, &basis, &state_density)?;
        let hole_population = atom_diagonal_population(&basis, &hole_density);
        let electron_population = atom_diagonal_population(&basis, &particle_density);
        let hole_centroid_angstrom = population_centroid_angstrom(mol, &hole_population);
        let electron_centroid_angstrom = population_centroid_angstrom(mol, &electron_population);
        let delta = [
            electron_centroid_angstrom[0] - hole_centroid_angstrom[0],
            electron_centroid_angstrom[1] - hole_centroid_angstrom[1],
            electron_centroid_angstrom[2] - hole_centroid_angstrom[2],
        ];
        let state_mu_au = permanent_dipole_au(mol, params, &basis, &state_density, &r_ao)?;
        let state_mu_debye = scale3(state_mu_au, DEBYE_PER_E_BOHR);
        let difference_mu_au = [
            state_mu_au[0] - ground_dipole_au[0],
            state_mu_au[1] - ground_dipole_au[1],
            state_mu_au[2] - ground_dipole_au[2],
        ];
        let difference_mu_debye = scale3(difference_mu_au, DEBYE_PER_E_BOHR);
        let (spin_multiplicity, s2_expectation) = match spin {
            CisSpin::Singlet => (1, 0.0),
            CisSpin::Triplet => (3, 2.0),
            CisSpin::Unrestricted => (0, f64::NAN),
        };
        states.push(ExcitedState {
            spin,
            spin_multiplicity,
            s2_expectation,
            energy_ev: e,
            state_total_energy_ev: ground.total_ev + e,
            state_total_energy_hartree: (ground.total_ev + e) / ZINDO_AU2EV,
            energy_cm1: e * EV_TO_WAVENUMBER_CM1,
            wavelength_nm: Some(EV_NM / e),
            oscillator_strength: fosc,
            transition_dipole_au: mu,
            transition_dipole_debye: scale3(mu, DEBYE_PER_E_BOHR),
            transition_dipole_magnitude_au: mu2.sqrt(),
            dipole_strength_au2: mu2,
            permanent_dipole_au: state_mu_au,
            permanent_dipole_debye: state_mu_debye,
            permanent_dipole_magnitude_debye: norm3(state_mu_debye),
            difference_dipole_au: difference_mu_au,
            difference_dipole_debye: difference_mu_debye,
            difference_dipole_magnitude_debye: norm3(difference_mu_debye),
            charges: state_charges,
            hole_population,
            electron_population,
            hole_centroid_angstrom,
            electron_centroid_angstrom,
            charge_transfer_distance_angstrom: norm3(delta),
            dominant,
        });
    }
    Ok(ZindoSpectrum {
        ground,
        spin,
        ground_dipole_au,
        ground_dipole_debye,
        states,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Vec3;
    use crate::system::Atom;

    #[test]
    fn zerner_fg_runtime_scaling_matches_mopac_switch() {
        let p = ZindoParameters::standard().unwrap();
        let c = p.element(6).unwrap();
        assert!((c.fg[14] - 6.897842 / 3.0).abs() < 1.0e-12);
        assert!((c.fg[15] - 4.509913 / 25.0).abs() < 1.0e-12);
    }

    #[test]
    fn gamma_zero_distance_is_atomic_hardness_for_equal_atoms() {
        for g in [4.0, 8.0, 12.85] {
            assert!((gamma1(g, g, 0.0) - g).abs() < 1.0e-12);
        }
    }

    #[test]
    fn water_zindo_s_smoke() {
        let p = ZindoParameters::standard().unwrap();
        let a = 1.0 / ZINDO_AU2ANG;
        let mol = Molecule {
            atoms: vec![
                Atom {
                    z: 8,
                    position: Vec3::new(0.0, 0.0, 0.0),
                },
                Atom {
                    z: 1,
                    position: Vec3::new(0.7586 * a, 0.0, 0.5043 * a),
                },
                Atom {
                    z: 1,
                    position: Vec3::new(-0.7586 * a, 0.0, 0.5043 * a),
                },
            ],
            charge: 0.0,
            multiplicity: 1,
        };
        let o = ZindoOptions {
            n_states: 3,
            ..ZindoOptions::default()
        };
        let r = run_zindo_s(&mol, &p, &o).unwrap();
        assert!(r.total_ev.is_finite());
        assert_eq!(r.n_occ, 4);
        let s = zindo_s_cis(&mol, &p, &o).unwrap();
        assert!(s.states.iter().all(|x| x.energy_ev.is_finite()));
        assert!(s.states.iter().all(|x| x.energy_cm1.is_finite()));
        assert!(s
            .states
            .iter()
            .all(|x| x.charge_transfer_distance_angstrom.is_finite()));
        let t = zindo_s_cis_spin(&mol, &p, &o, CisSpin::Triplet).unwrap();
        assert!(t.states.iter().all(|x| x.oscillator_strength == 0.0));

        let uv = uv_vis_spectrum(&mol, &p, &o).unwrap();
        assert_eq!(uv.spin, CisSpin::Singlet);
        assert!(uv.states.iter().all(|x| x.wavelength_nm.is_some()));
        assert!(uv.states.iter().all(|x| x.oscillator_strength >= 0.0));
    }

    #[test]
    fn d_branch_is_explicitly_rejected() {
        let p = ZindoParameters::standard().unwrap();
        let mol = Molecule {
            atoms: vec![Atom {
                z: 26,
                position: Vec3::zero(),
            }],
            charge: 0.0,
            multiplicity: 1,
        };
        let e = run_zindo_s(&mol, &p, &ZindoOptions::default()).unwrap_err();
        assert!(e.to_string().contains("9-AO"));
    }
}
