// SPDX-License-Identifier: GPL-3.0-or-later

//! MolDS-compatible CNDO/2 and ground-state INDO single-point engines.
//!
//! This module is deliberately independent of Zerner INDO/S (`ZINDO/S`).
//! The numerical conventions and the H/Li/C/N/O/S atomic data are derived from
//! MolDS 0.3.1 (GPL-3.0-or-later); see `third_party/molds/NOTICE` and
//! `src/data/molds_cndo2_indo_parameters.csv`.
//!
//! Both methods use the ZDO secular equation `F C = C eps`; AO overlap enters
//! only the interatomic resonance integral. CNDO/2 is enabled for
//! H/Li/C/N/O/S and INDO for H/Li/C/N/O, exactly matching MolDS's enabled
//! atom lists. Sulfur is not accepted in INDO because its shared row lacks the
//! complete INDO one-centre coefficients.

// AO, atom, and Cartesian indices are intentionally explicit in the fixed-size
// matrix contractions below; iterator rewrites obscure the equations.
#![allow(clippy::needless_range_loop)]

use crate::constants::HARTREE_TO_EV;
use crate::data_tables::MOLDS_CNDO2_INDO_PARAM_CSV;
use crate::dual::{Dual, Scalar};
use crate::dual2::Dual2;
use crate::error::{Result, XndoError};
use crate::linalg::{symmetric_eigen, Matrix};
use crate::math::Vec3;
use crate::method::Method;
use crate::overlap_numeric::overlap_sto;
use crate::rotations::Rotation;
use crate::scf::Reference;
use crate::system::Molecule;
use crate::zdo_gradient::{
    atom_population, exchange_weight_rhf, exchange_weight_uhf, pair_energy, solve_rhf_responses,
    solve_uhf_responses, spin_densities,
};

const SPD_L: [u8; 9] = [0, 1, 1, 1, 2, 2, 2, 2, 2];
const SPD_M: [i8; 9] = [0, 0, 1, -1, 0, 1, -1, 2, -2];

#[derive(Clone, Debug)]
pub struct CndoIndoElement {
    pub z: u8,
    pub symbol: String,
    pub core_charge: f64,
    pub valence_electrons: usize,
    pub valence_shell: u8,
    pub n_orb: usize,
    pub has_d: bool,
    pub bonding_parameter_ev: f64,
    pub imu_s_ev: f64,
    pub imu_p_ev: f64,
    pub imu_d_ev: f64,
    pub zeta_s: f64,
    pub zeta_p: f64,
    pub zeta_d: f64,
    /// INDO Slater-Condon G1, converted from the MolDS native atomic-unit value.
    pub indo_g1_ev: f64,
    /// INDO Slater-Condon F2, converted from the MolDS native atomic-unit value.
    pub indo_f2_ev: f64,
    pub indo_f0_coeff_s: f64,
    pub indo_f0_coeff_p: f64,
    pub indo_g1_coeff_s: f64,
    pub indo_g1_coeff_p: f64,
    pub indo_f2_coeff_s: f64,
    pub indo_f2_coeff_p: f64,
}

impl CndoIndoElement {
    #[inline]
    fn principal_for_l(&self, l: u8) -> u8 {
        match l {
            0 | 1 => self.valence_shell,
            2 => self.valence_shell,
            _ => 0,
        }
    }

    #[inline]
    fn zeta_for_l(&self, l: u8) -> f64 {
        match l {
            0 => self.zeta_s,
            1 => self.zeta_p,
            2 => self.zeta_d,
            _ => 0.0,
        }
    }

    #[inline]
    fn imu_for_orb(&self, orb: usize) -> f64 {
        match orb {
            0 => self.imu_s_ev,
            1..=3 => self.imu_p_ev,
            _ => self.imu_d_ev,
        }
    }

    #[inline]
    fn indo_coeffs_for_orb(&self, orb: usize) -> (f64, f64, f64) {
        if orb == 0 {
            (
                self.indo_f0_coeff_s,
                self.indo_g1_coeff_s,
                self.indo_f2_coeff_s,
            )
        } else {
            (
                self.indo_f0_coeff_p,
                self.indo_g1_coeff_p,
                self.indo_f2_coeff_p,
            )
        }
    }
}

fn parse_num<T: std::str::FromStr>(text: &str, field: &str, z: u8) -> Result<T> {
    text.parse::<T>().map_err(|_| {
        XndoError::MissingParameter(format!(
            "invalid MolDS CNDO2/INDO field {field} for Z={z}: {text:?}"
        ))
    })
}

/// Read a MolDS CNDO2/INDO element row from the bundled machine-readable table.
pub fn element(z: u8) -> Result<CndoIndoElement> {
    for line in MOLDS_CNDO2_INDO_PARAM_CSV.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with("z,") {
            continue;
        }
        let f: Vec<&str> = line.split(',').collect();
        if f.len() < 25 {
            continue;
        }
        let row_z: u8 = match f[0].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if row_z != z {
            continue;
        }
        let shell = match f[4].trim().to_ascii_lowercase().as_str() {
            "k" => 1,
            "l" => 2,
            "m" => 3,
            other => {
                return Err(XndoError::MissingParameter(format!(
                    "invalid MolDS valence shell {other:?} for Z={z}"
                )))
            }
        };
        let has_p = parse_num::<u8>(f[6], "has_p", z)? != 0;
        let has_d = parse_num::<u8>(f[7], "has_d_cndo_indo", z)? != 0;
        let n_orb = if has_d {
            9
        } else if has_p {
            4
        } else {
            1
        };
        let z_eff_k = parse_num::<f64>(f[12], "z_eff_k", z)?;
        let z_eff_l = parse_num::<f64>(f[13], "z_eff_l", z)?;
        let z_eff_msp = parse_num::<f64>(f[14], "z_eff_msp", z)?;
        let z_eff_md = parse_num::<f64>(f[15], "z_eff_md", z)?;
        let zeta_sp = match shell {
            1 => z_eff_k,
            2 => z_eff_l / 2.0,
            3 => z_eff_msp / 3.0,
            _ => unreachable!(),
        };
        let zeta_d = if has_d { z_eff_md / shell as f64 } else { 0.0 };
        return Ok(CndoIndoElement {
            z,
            symbol: f[1].to_string(),
            core_charge: parse_num(f[2], "core_charge", z)?,
            valence_electrons: parse_num(f[3], "valence_electrons", z)?,
            valence_shell: shell,
            n_orb,
            has_d,
            bonding_parameter_ev: parse_num(f[8], "bonding_parameter_ev", z)?,
            imu_s_ev: parse_num(f[9], "imu_amu_s_ev", z)?,
            imu_p_ev: parse_num(f[10], "imu_amu_p_ev", z)?,
            imu_d_ev: parse_num(f[11], "imu_amu_d_ev", z)?,
            zeta_s: zeta_sp,
            zeta_p: if has_p { zeta_sp } else { 0.0 },
            zeta_d,
            indo_g1_ev: parse_num::<f64>(f[16], "indo_g1_native", z)? * HARTREE_TO_EV,
            indo_f2_ev: parse_num::<f64>(f[17], "indo_f2_native", z)? * HARTREE_TO_EV,
            indo_f0_coeff_s: parse_num(f[18], "indo_f0_coeff_s", z)?,
            indo_f0_coeff_p: parse_num(f[19], "indo_f0_coeff_p", z)?,
            indo_g1_coeff_s: parse_num(f[20], "indo_g1_coeff_s", z)?,
            indo_g1_coeff_p: parse_num(f[21], "indo_g1_coeff_p", z)?,
            indo_f2_coeff_s: parse_num(f[22], "indo_f2_coeff_s", z)?,
            indo_f2_coeff_p: parse_num(f[23], "indo_f2_coeff_p", z)?,
        });
    }
    Err(XndoError::MissingParameter(format!(
        "MolDS CNDO2/INDO element Z={z}; bundled donor subset is H, Li, C, N, O, S"
    )))
}

fn validate_element(method: Method, z: u8) -> Result<CndoIndoElement> {
    let e = element(z)?;
    match method {
        Method::Cndo2 => Ok(e),
        Method::Indo => {
            if matches!(z, 1 | 3 | 6 | 7 | 8) {
                Ok(e)
            } else {
                Err(XndoError::MissingParameter(format!(
                    "MolDS-compatible INDO in xndo-rs 0.2.4 is enabled for H/Li/C/N/O; Z={z} is not enabled because its complete INDO one-centre parameter row is not bundled"
                )))
            }
        }
        _ => Err(XndoError::InvalidInput(format!(
            "{} is not handled by the MolDS-compatible CNDO2/INDO engine",
            method
        ))),
    }
}

#[derive(Clone, Debug)]
pub struct CndoIndoOptions {
    pub charge: f64,
    pub multiplicity: usize,
    pub reference: Reference,
    pub max_scf: usize,
    pub e_tol_ev: f64,
    pub p_tol: f64,
    pub damping: f64,
}

impl Default for CndoIndoOptions {
    fn default() -> Self {
        Self {
            charge: 0.0,
            multiplicity: 1,
            reference: Reference::Auto,
            max_scf: 300,
            e_tol_ev: 1.0e-8,
            p_tol: 1.0e-7,
            damping: 0.20,
        }
    }
}

#[derive(Clone, Debug)]
pub struct CndoIndoResult {
    pub method: Method,
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

#[derive(Clone, Copy, Debug)]
struct Ao {
    atom: usize,
    orb: usize,
}

#[derive(Clone, Debug)]
struct Basis {
    aos: Vec<Ao>,
    offsets: Vec<usize>,
    norb: Vec<usize>,
    elements: Vec<CndoIndoElement>,
}

impl Basis {
    fn build(m: &Molecule, method: Method) -> Result<Self> {
        let mut aos = Vec::new();
        let mut offsets = Vec::with_capacity(m.atoms.len());
        let mut norb = Vec::with_capacity(m.atoms.len());
        let mut elements = Vec::with_capacity(m.atoms.len());
        for (ia, atom) in m.atoms.iter().enumerate() {
            let e = validate_element(method, atom.z)?;
            offsets.push(aos.len());
            norb.push(e.n_orb);
            for orb in 0..e.n_orb {
                aos.push(Ao { atom: ia, orb });
            }
            elements.push(e);
        }
        Ok(Self {
            aos,
            offsets,
            norb,
            elements,
        })
    }

    #[inline]
    fn nao(&self) -> usize {
        self.aos.len()
    }
}

#[inline]
fn orbital_l(orb: usize) -> u8 {
    SPD_L[orb]
}

fn overlap_block_g<S: Scalar>(
    ea: &CndoIndoElement,
    eb: &CndoIndoElement,
    dvec: [S; 3],
) -> [[S; 9]; 9] {
    let r = (dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2]).sqrt();
    let mut out = [[S::cst(0.0); 9]; 9];
    if r.val() < 1.0e-14 {
        for i in 0..ea.n_orb.min(eb.n_orb) {
            out[i][i] = S::cst(1.0);
        }
        return out;
    }
    let rot = Rotation::<S>::build(dvec, ea.has_d || eb.has_d);
    let mut local = [[S::cst(0.0); 9]; 9];
    for ia in 0..ea.n_orb {
        for ib in 0..eb.n_orb {
            if SPD_M[ia] != SPD_M[ib] {
                continue;
            }
            let la = orbital_l(ia);
            let lb = orbital_l(ib);
            let m = SPD_M[ia].unsigned_abs();
            local[ia][ib] = overlap_sto(
                ea.principal_for_l(la),
                ea.zeta_for_l(la),
                la,
                m,
                eb.principal_for_l(lb),
                eb.zeta_for_l(lb),
                lb,
                r,
            );
        }
    }
    for mu in 0..ea.n_orb {
        for nu in 0..eb.n_orb {
            let mut v = S::cst(0.0);
            for a in 0..ea.n_orb {
                let ca = rot.coeff(mu, a);
                for b in 0..eb.n_orb {
                    v = v + ca * rot.coeff(nu, b) * local[a][b];
                }
            }
            out[mu][nu] = v;
        }
    }
    out
}

fn overlap_block(ea: &CndoIndoElement, eb: &CndoIndoElement, dvec: [f64; 3]) -> [[f64; 9]; 9] {
    overlap_block_g(ea, eb, dvec)
}

fn build_overlap(m: &Molecule, b: &Basis) -> Result<Matrix> {
    let mut s = Matrix::identity(b.nao());
    for ia in 0..m.atoms.len() {
        for ja in ia + 1..m.atoms.len() {
            let d = m.atoms[ja].position - m.atoms[ia].position;
            if d.norm() < 1.0e-12 {
                return Err(XndoError::InvalidInput(
                    "coincident atoms in CNDO/2 or INDO".into(),
                ));
            }
            let block = overlap_block(&b.elements[ia], &b.elements[ja], [d.x, d.y, d.z]);
            for i in 0..b.norb[ia] {
                for j in 0..b.norb[ja] {
                    let mu = b.offsets[ia] + i;
                    let nu = b.offsets[ja] + j;
                    s[(mu, nu)] = block[i][j];
                    s[(nu, mu)] = block[i][j];
                }
            }
        }
    }
    Ok(s)
}

#[inline]
fn factorial(n: usize) -> f64 {
    (2..=n).fold(1.0, |a, x| a * x as f64)
}

#[inline]
fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut v = 1.0;
    for i in 0..k {
        v *= (n - i) as f64 / (i + 1) as f64;
    }
    v
}

/// Coefficient of eta^k in (xi+eta)^a (xi-eta)^b.
fn z_coefficient(a: usize, b: usize, k: usize) -> f64 {
    let mut v = 0.0;
    for i in 0..=a {
        if k < i {
            continue;
        }
        let j = k - i;
        if j > b {
            continue;
        }
        let sign = if j.is_multiple_of(2) { 1.0 } else { -1.0 };
        v += binomial(a, i) * binomial(b, j) * sign;
    }
    v
}

/// Mulliken auxiliary A_k(rho) = integral_1^inf x^k exp(-rho*x) dx.
fn aux_a<S: Scalar>(k: usize, rho: S) -> S {
    debug_assert!(rho.val() > 0.0);
    let e = (-rho).exp();
    let mut v = e / rho;
    for n in 1..=k {
        v = e / rho + v * (n as f64) / rho;
    }
    v
}

/// Mulliken auxiliary B_k(rho) = integral_-1^1 x^k exp(-rho*x) dx.
fn aux_b<S: Scalar>(k: usize, rho: S) -> S {
    if rho.val().abs() < 1.0 {
        // Entire power series.  This avoids the severe small-rho cancellation
        // of the upward recurrence and is inexpensive for the low orders here.
        let mut sum = S::cst(0.0);
        let mut coeff = S::cst(1.0); // (-rho)^m / m!
        for m in 0..120usize {
            if (k + m).is_multiple_of(2) {
                sum = sum + coeff * 2.0 / (k + m + 1) as f64;
            }
            coeff = coeff * (-rho) / (m + 1) as f64;
            if coeff.val().abs() < 1.0e-18 {
                break;
            }
        }
        return sum;
    }
    let ep = rho.exp();
    let em = (-rho).exp();
    let mut v = (ep - em) / rho;
    for n in 1..=k {
        let boundary = if n % 2 == 0 { ep - em } else { -ep - em };
        v = boundary / rho + v * (n as f64) / rho;
    }
    v
}

fn reduced_overlap<S: Scalar>(a: usize, b: usize, alpha: S, beta: S) -> S {
    let p = (alpha + beta) * 0.5;
    let q = (alpha - beta) * 0.5;
    let mut sum = S::cst(0.0);
    for k in 0..=a + b {
        sum = sum + aux_a(a + b - k, p) * aux_b(k, q) * z_coefficient(a, b, k);
    }
    sum * 0.5
}

fn gamma_zero_one_way(na: usize, za: f64, nb: usize, zb: f64) -> f64 {
    let mut value = factorial(2 * na - 1) / (2.0 * za).powi((2 * na) as i32);
    for l in 1..=2 * nb {
        let mut temp = l as f64 * (2.0 * zb).powi((2 * nb - l) as i32);
        temp *= factorial(2 * na + 2 * nb - l - 1);
        temp /= factorial(2 * nb - l) * (2 * nb) as f64;
        temp /= (2.0 * za + 2.0 * zb).powi((2 * (na + nb) - l) as i32);
        value -= temp;
    }
    value * (2.0 * za).powi((2 * na + 1) as i32) / factorial(2 * na)
}

fn gamma_r_one_way<S: Scalar>(na: usize, za: f64, nb: usize, zb: f64, r: S) -> S {
    if r.val() < 1.0e-10 {
        return S::cst(gamma_zero_one_way(na, za, nb, zb));
    }
    let half_r = r * 0.5;
    let mut value =
        half_r.powi((2 * na) as i32) * reduced_overlap(2 * na - 1, 0, r * (2.0 * za), S::cst(0.0));
    for l in 1..=2 * nb {
        let mut temp = l as f64 * (2.0 * zb).powi((2 * nb - l) as i32);
        temp /= factorial(2 * nb - l) * (2 * nb) as f64;
        let term = half_r.powi((2 * nb - l + 2 * na) as i32)
            * reduced_overlap(2 * na - 1, 2 * nb - l, r * (2.0 * za), r * (2.0 * zb))
            * temp;
        value = value - term;
    }
    value * (2.0 * za).powi((2 * na + 1) as i32) / factorial(2 * na)
}

/// CNDO/INDO gamma_AB in eV, using the MolDS/Pople Slater-charge expression.
fn gamma_ev<S: Scalar>(a: &CndoIndoElement, b: &CndoIndoElement, r_bohr: S) -> S {
    let na = a.valence_shell as usize;
    let nb = b.valence_shell as usize;
    let za = a.zeta_s;
    let zb = b.zeta_s;
    let gab = gamma_r_one_way(na, za, nb, zb, r_bohr);
    let gba = gamma_r_one_way(nb, zb, na, za, r_bohr);
    (gab + gba) * (0.5 * HARTREE_TO_EV)
}

fn gamma_matrix(m: &Molecule, b: &Basis) -> Result<Vec<Vec<f64>>> {
    let n = m.atoms.len();
    let mut g = vec![vec![0.0; n]; n];
    for ia in 0..n {
        g[ia][ia] = gamma_ev(&b.elements[ia], &b.elements[ia], 0.0);
        for ja in ia + 1..n {
            let r = (m.atoms[ia].position - m.atoms[ja].position).norm();
            if r < 1.0e-12 {
                return Err(XndoError::InvalidInput(
                    "coincident atoms in CNDO/2 or INDO".into(),
                ));
            }
            let v = gamma_ev(&b.elements[ia], &b.elements[ja], r);
            g[ia][ja] = v;
            g[ja][ia] = v;
        }
    }
    Ok(g)
}

#[inline]
fn bonding_adjust_k(a: &CndoIndoElement, b: &CndoIndoElement) -> f64 {
    if a.valence_shell >= 3 || b.valence_shell >= 3 {
        0.75
    } else {
        1.0
    }
}

fn build_guess_fock(m: &Molecule, b: &Basis, s: &Matrix) -> Matrix {
    let mut f = Matrix::zeros(b.nao(), b.nao());
    for mu in 0..b.nao() {
        let ao = b.aos[mu];
        f[(mu, mu)] = -b.elements[ao.atom].imu_for_orb(ao.orb);
        for nu in 0..mu {
            let bo = b.aos[nu];
            if ao.atom == bo.atom {
                continue;
            }
            let ea = &b.elements[ao.atom];
            let eb = &b.elements[bo.atom];
            let beta = 0.5
                * bonding_adjust_k(ea, eb)
                * (ea.bonding_parameter_ev + eb.bonding_parameter_ev);
            let v = beta * s[(mu, nu)];
            f[(mu, nu)] = v;
            f[(nu, mu)] = v;
        }
    }
    // m is intentionally kept in the signature to mirror build_hcore and make
    // the atom ordering explicit; no geometry term beyond S is required here.
    let _ = m;
    f
}

fn build_hcore(method: Method, m: &Molecule, b: &Basis, s: &Matrix, g: &[Vec<f64>]) -> Matrix {
    let mut h = Matrix::zeros(b.nao(), b.nao());
    for mu in 0..b.nao() {
        let ao = b.aos[mu];
        let ea = &b.elements[ao.atom];
        let gaa = g[ao.atom][ao.atom];
        let self_core = match method {
            Method::Cndo2 => -(ea.core_charge - 0.5) * gaa,
            Method::Indo => {
                let (cf0, cg1, cf2) = ea.indo_coeffs_for_orb(ao.orb);
                -(cf0 * gaa + cg1 * ea.indo_g1_ev + cf2 * ea.indo_f2_ev)
            }
            _ => unreachable!(),
        };
        let mut nuc_other = 0.0;
        for ib in 0..m.atoms.len() {
            if ib != ao.atom {
                nuc_other += b.elements[ib].core_charge * g[ao.atom][ib];
            }
        }
        h[(mu, mu)] = -ea.imu_for_orb(ao.orb) + self_core - nuc_other;

        for nu in 0..mu {
            let bo = b.aos[nu];
            if ao.atom == bo.atom {
                continue;
            }
            let eb = &b.elements[bo.atom];
            let beta = 0.5
                * bonding_adjust_k(ea, eb)
                * (ea.bonding_parameter_ev + eb.bonding_parameter_ev);
            let v = beta * s[(mu, nu)];
            h[(mu, nu)] = v;
            h[(nu, mu)] = v;
        }
    }
    h
}

fn atom_populations(b: &Basis, p: &Matrix) -> Vec<f64> {
    let mut q = vec![0.0; b.elements.len()];
    for ia in 0..b.elements.len() {
        for i in 0..b.norb[ia] {
            q[ia] += p[(b.offsets[ia] + i, b.offsets[ia] + i)];
        }
    }
    q
}

/// INDO one-centre Coulomb and exchange integrals in eV for s/p AOs.
fn indo_one_center_jk(e: &CndoIndoElement, mu: usize, nu: usize, gamma_aa: f64) -> (f64, f64) {
    if mu == 0 && nu == 0 {
        return (gamma_aa, gamma_aa);
    }
    if mu == 0 || nu == 0 {
        return (gamma_aa, e.indo_g1_ev / 3.0);
    }
    if mu == nu {
        let j = gamma_aa + 4.0 * e.indo_f2_ev / 25.0;
        return (j, j);
    }
    (
        gamma_aa - 2.0 * e.indo_f2_ev / 25.0,
        3.0 * e.indo_f2_ev / 25.0,
    )
}

fn build_fock_cndo_rhf(b: &Basis, h: &Matrix, g: &[Vec<f64>], p: &Matrix) -> Matrix {
    let pop = atom_populations(b, p);
    let mut f = h.clone();
    for mu in 0..b.nao() {
        let a = b.aos[mu].atom;
        let mut v = h[(mu, mu)] + pop[a] * g[a][a] - 0.5 * p[(mu, mu)] * g[a][a];
        for ib in 0..b.elements.len() {
            if ib != a {
                v += pop[ib] * g[a][ib];
            }
        }
        f[(mu, mu)] = v;
        for nu in 0..mu {
            let c = b.aos[nu].atom;
            let gab = g[a][c];
            let v = h[(mu, nu)] - 0.5 * p[(mu, nu)] * gab;
            f[(mu, nu)] = v;
            f[(nu, mu)] = v;
        }
    }
    f
}

fn build_fock_cndo_uhf(
    b: &Basis,
    h: &Matrix,
    g: &[Vec<f64>],
    pa: &Matrix,
    pb: &Matrix,
) -> (Matrix, Matrix) {
    let pt = add(pa, pb);
    let pop = atom_populations(b, &pt);
    let mut fa = h.clone();
    let mut fb = h.clone();
    for mu in 0..b.nao() {
        let a = b.aos[mu].atom;
        let mut common = h[(mu, mu)] + pop[a] * g[a][a];
        for ib in 0..b.elements.len() {
            if ib != a {
                common += pop[ib] * g[a][ib];
            }
        }
        fa[(mu, mu)] = common - pa[(mu, mu)] * g[a][a];
        fb[(mu, mu)] = common - pb[(mu, mu)] * g[a][a];
        for nu in 0..mu {
            let c = b.aos[nu].atom;
            let gab = g[a][c];
            let va = h[(mu, nu)] - pa[(mu, nu)] * gab;
            let vb = h[(mu, nu)] - pb[(mu, nu)] * gab;
            fa[(mu, nu)] = va;
            fa[(nu, mu)] = va;
            fb[(mu, nu)] = vb;
            fb[(nu, mu)] = vb;
        }
    }
    (fa, fb)
}

fn build_fock_indo_rhf(b: &Basis, h: &Matrix, g: &[Vec<f64>], p: &Matrix) -> Matrix {
    let pop = atom_populations(b, p);
    let mut f = h.clone();
    for mu in 0..b.nao() {
        let ao = b.aos[mu];
        let a = ao.atom;
        let ea = &b.elements[a];
        let mut v = h[(mu, mu)];
        for l in 0..b.norb[a] {
            let lam = b.offsets[a] + l;
            let (j, k) = indo_one_center_jk(ea, ao.orb, l, g[a][a]);
            v += p[(lam, lam)] * (j - 0.5 * k);
        }
        for ib in 0..b.elements.len() {
            if ib != a {
                v += pop[ib] * g[a][ib];
            }
        }
        f[(mu, mu)] = v;
        for nu in 0..mu {
            let bo = b.aos[nu];
            let c = bo.atom;
            let x = if a == c {
                let (j, k) = indo_one_center_jk(ea, ao.orb, bo.orb, g[a][a]);
                (1.5 * k - 0.5 * j) * p[(mu, nu)]
            } else {
                h[(mu, nu)] - 0.5 * p[(mu, nu)] * g[a][c]
            };
            f[(mu, nu)] = x;
            f[(nu, mu)] = x;
        }
    }
    f
}

fn build_fock_indo_uhf(
    b: &Basis,
    h: &Matrix,
    g: &[Vec<f64>],
    pa: &Matrix,
    pb: &Matrix,
) -> (Matrix, Matrix) {
    let pt = add(pa, pb);
    let pop = atom_populations(b, &pt);
    let mut fa = h.clone();
    let mut fb = h.clone();
    for mu in 0..b.nao() {
        let ao = b.aos[mu];
        let a = ao.atom;
        let ea = &b.elements[a];
        let mut va = h[(mu, mu)];
        let mut vb = h[(mu, mu)];
        for l in 0..b.norb[a] {
            let lam = b.offsets[a] + l;
            let (j, k) = indo_one_center_jk(ea, ao.orb, l, g[a][a]);
            va += pt[(lam, lam)] * j - pa[(lam, lam)] * k;
            vb += pt[(lam, lam)] * j - pb[(lam, lam)] * k;
        }
        for ib in 0..b.elements.len() {
            if ib != a {
                va += pop[ib] * g[a][ib];
                vb += pop[ib] * g[a][ib];
            }
        }
        fa[(mu, mu)] = va;
        fb[(mu, mu)] = vb;

        for nu in 0..mu {
            let bo = b.aos[nu];
            let c = bo.atom;
            let (xa, xb) = if a == c {
                let (j, k) = indo_one_center_jk(ea, ao.orb, bo.orb, g[a][a]);
                (
                    2.0 * k * pt[(mu, nu)] - (j + k) * pa[(mu, nu)],
                    2.0 * k * pt[(mu, nu)] - (j + k) * pb[(mu, nu)],
                )
            } else {
                (
                    h[(mu, nu)] - pa[(mu, nu)] * g[a][c],
                    h[(mu, nu)] - pb[(mu, nu)] * g[a][c],
                )
            };
            fa[(mu, nu)] = xa;
            fa[(nu, mu)] = xa;
            fb[(mu, nu)] = xb;
            fb[(nu, mu)] = xb;
        }
    }
    (fa, fb)
}

fn build_fock_rhf(method: Method, b: &Basis, h: &Matrix, g: &[Vec<f64>], p: &Matrix) -> Matrix {
    match method {
        Method::Cndo2 => build_fock_cndo_rhf(b, h, g, p),
        Method::Indo => build_fock_indo_rhf(b, h, g, p),
        _ => unreachable!(),
    }
}

fn build_fock_uhf(
    method: Method,
    b: &Basis,
    h: &Matrix,
    g: &[Vec<f64>],
    pa: &Matrix,
    pb: &Matrix,
) -> (Matrix, Matrix) {
    match method {
        Method::Cndo2 => build_fock_cndo_uhf(b, h, g, pa, pb),
        Method::Indo => build_fock_indo_uhf(b, h, g, pa, pb),
        _ => unreachable!(),
    }
}

fn electron_count(b: &Basis, charge: f64) -> Result<usize> {
    let q = charge.round();
    if (charge - q).abs() > 1.0e-8 {
        return Err(XndoError::InvalidInput(
            "CNDO/2 and INDO charge must be integral".into(),
        ));
    }
    let mut ne = -(q as i64);
    for e in &b.elements {
        ne += e.valence_electrons as i64;
    }
    if ne < 0 {
        return Err(XndoError::InvalidInput(
            "negative CNDO/2 or INDO electron count".into(),
        ));
    }
    Ok(ne as usize)
}

fn electron_partition(ne: usize, multiplicity: usize) -> Result<(usize, usize)> {
    if multiplicity == 0 {
        return Err(XndoError::InvalidInput(
            "CNDO/2 and INDO multiplicity must be >= 1".into(),
        ));
    }
    let unpaired = multiplicity - 1;
    if unpaired > ne || !(ne - unpaired).is_multiple_of(2) {
        return Err(XndoError::InvalidInput(format!(
            "electron count {ne} is incompatible with multiplicity {multiplicity}"
        )));
    }
    let n_beta = (ne - unpaired) / 2;
    Ok((n_beta + unpaired, n_beta))
}

fn core_energy(m: &Molecule, b: &Basis) -> Result<f64> {
    let mut e = 0.0;
    for ia in 0..m.atoms.len() {
        for ja in ia + 1..m.atoms.len() {
            let r = (m.atoms[ia].position - m.atoms[ja].position).norm();
            if r < 1.0e-12 {
                return Err(XndoError::InvalidInput(
                    "coincident atoms in CNDO/2 or INDO".into(),
                ));
            }
            e += b.elements[ia].core_charge * b.elements[ja].core_charge / r * HARTREE_TO_EV;
        }
    }
    Ok(e)
}

fn add(a: &Matrix, b: &Matrix) -> Matrix {
    let mut c = Matrix::zeros(a.rows, a.cols);
    for i in 0..c.as_slice().len() {
        c.as_mut_slice()[i] = a.as_slice()[i] + b.as_slice()[i];
    }
    c
}

fn sub(a: &Matrix, b: &Matrix) -> Matrix {
    let mut c = Matrix::zeros(a.rows, a.cols);
    for i in 0..c.as_slice().len() {
        c.as_mut_slice()[i] = a.as_slice()[i] - b.as_slice()[i];
    }
    c
}

fn damp_density(raw: &Matrix, old: &Matrix, damping: f64) -> Matrix {
    if damping <= 0.0 {
        return raw.clone();
    }
    let mut out = raw.clone();
    for i in 0..out.as_slice().len() {
        out.as_mut_slice()[i] = (1.0 - damping) * raw.as_slice()[i] + damping * old.as_slice()[i];
    }
    out
}

fn charges(b: &Basis, p: &Matrix) -> Vec<f64> {
    let pop = atom_populations(b, p);
    b.elements
        .iter()
        .zip(pop)
        .map(|(e, n)| e.core_charge - n)
        .collect()
}

fn initial_densities(
    guess: &Matrix,
    n_alpha: usize,
    n_beta: usize,
    unrestricted: bool,
) -> Result<(Matrix, Option<Matrix>)> {
    let (_, c0) = symmetric_eigen(guess)?;
    if unrestricted {
        Ok((
            c0.leading_columns_gram(n_alpha, 1.0),
            Some(c0.leading_columns_gram(n_beta, 1.0)),
        ))
    } else {
        Ok((c0.leading_columns_gram(n_alpha, 2.0), None))
    }
}

/// Execute MolDS-compatible CNDO/2 or ground-state INDO.
///
/// `Method::Cndo2` supports H/Li/C/N/O/S with the bundled MolDS donor table.
/// `Method::Indo` supports H/Li/C/N/O. Both RHF and a spin-resolved UHF extension
/// are available. The UHF equations reduce exactly to the donor RHF equations
/// when `P_alpha = P_beta = P/2`.
pub fn run_cndo_indo(
    m: &Molecule,
    method: Method,
    opt: &CndoIndoOptions,
) -> Result<CndoIndoResult> {
    if !matches!(method, Method::Cndo2 | Method::Indo) {
        return Err(XndoError::InvalidInput(format!(
            "{} is not CNDO/2 or ground-state INDO",
            method
        )));
    }
    if m.atoms.is_empty() {
        return Err(XndoError::InvalidInput(
            "CNDO/2 or INDO requires at least one atom".into(),
        ));
    }
    let b = Basis::build(m, method)?;
    let ne = electron_count(&b, opt.charge)?;
    let (n_alpha, n_beta) = electron_partition(ne, opt.multiplicity)?;
    if n_alpha > b.nao() || n_beta > b.nao() {
        return Err(XndoError::InvalidInput(format!(
            "{} electron count exceeds the {}-AO valence basis",
            method,
            b.nao()
        )));
    }
    let reference = match opt.reference {
        Reference::Auto => {
            if n_alpha == n_beta {
                Reference::Rhf
            } else {
                Reference::Uhf
            }
        }
        x => x,
    };
    if reference == Reference::Rhf && n_alpha != n_beta {
        return Err(XndoError::InvalidInput(format!(
            "{} RHF is incompatible with {n_alpha} alpha and {n_beta} beta electrons; use reference=UHF or Auto",
            method
        )));
    }

    let s = build_overlap(m, &b)?;
    let g = gamma_matrix(m, &b)?;
    let h = build_hcore(method, m, &b, &s, &g);
    let guess = build_guess_fock(m, &b, &s);
    let core = core_energy(m, &b)?;
    let damping = opt.damping.clamp(0.0, 0.95);
    let mut last_energy = f64::INFINITY;
    let mut last_error = f64::INFINITY;

    if reference == Reference::Rhf {
        let n_occ = n_alpha;
        let (mut p, _) = initial_densities(&guess, n_alpha, n_beta, false)?;
        for it in 1..=opt.max_scf {
            let f = build_fock_rhf(method, &b, &h, &g, &p);
            let (_, c) = symmetric_eigen(&f)?;
            let raw = c.leading_columns_gram(n_occ, 2.0);
            let p_next = damp_density(&raw, &p, damping);
            let f_next = build_fock_rhf(method, &b, &h, &g, &p_next);
            let electronic = 0.5 * p_next.frobenius_dot(&add(&h, &f_next));
            let p_error = p_next.rms_difference(&p);
            let e_error = (electronic - last_energy).abs();
            p = p_next;
            last_energy = electronic;
            last_error = p_error;
            if p_error < opt.p_tol && e_error < opt.e_tol_ev {
                let f_final = build_fock_rhf(method, &b, &h, &g, &p);
                let (eps, coeff) = symmetric_eigen(&f_final)?;
                let p_final = coeff.leading_columns_gram(n_occ, 2.0);
                let f_canonical = build_fock_rhf(method, &b, &h, &g, &p_final);
                let electronic_final = 0.5 * p_final.frobenius_dot(&add(&h, &f_canonical));
                return Ok(CndoIndoResult {
                    method,
                    density: p_final.clone(),
                    fock: f_canonical,
                    fock_beta: None,
                    mo_coeff: coeff,
                    mo_coeff_beta: None,
                    mo_energies_ev: eps,
                    mo_energies_beta_ev: None,
                    n_occ,
                    n_alpha,
                    n_beta,
                    spin_density: None,
                    unrestricted: false,
                    electronic_ev: electronic_final,
                    core_ev: core,
                    total_ev: electronic_final + core,
                    charges: charges(&b, &p_final),
                    iterations: it,
                    converged: true,
                });
            }
        }
    } else {
        let (mut pa, pb0) = initial_densities(&guess, n_alpha, n_beta, true)?;
        let mut pb = pb0.expect("UHF beta density must exist");
        for it in 1..=opt.max_scf {
            let (fa, fb) = build_fock_uhf(method, &b, &h, &g, &pa, &pb);
            let (_, ca) = symmetric_eigen(&fa)?;
            let (_, cb) = symmetric_eigen(&fb)?;
            let pa_raw = ca.leading_columns_gram(n_alpha, 1.0);
            let pb_raw = cb.leading_columns_gram(n_beta, 1.0);
            let pa_next = damp_density(&pa_raw, &pa, damping);
            let pb_next = damp_density(&pb_raw, &pb, damping);
            let (fa_next, fb_next) = build_fock_uhf(method, &b, &h, &g, &pa_next, &pb_next);
            let pt_next = add(&pa_next, &pb_next);
            let electronic = 0.5
                * (pt_next.frobenius_dot(&h)
                    + pa_next.frobenius_dot(&fa_next)
                    + pb_next.frobenius_dot(&fb_next));
            let p_error = pa_next.rms_difference(&pa).max(pb_next.rms_difference(&pb));
            let e_error = (electronic - last_energy).abs();
            pa = pa_next;
            pb = pb_next;
            last_energy = electronic;
            last_error = p_error;
            if p_error < opt.p_tol && e_error < opt.e_tol_ev {
                let (fa0, fb0) = build_fock_uhf(method, &b, &h, &g, &pa, &pb);
                let (eps_a, ca) = symmetric_eigen(&fa0)?;
                let (eps_b, cb) = symmetric_eigen(&fb0)?;
                let pa_final = ca.leading_columns_gram(n_alpha, 1.0);
                let pb_final = cb.leading_columns_gram(n_beta, 1.0);
                let pt_final = add(&pa_final, &pb_final);
                let (fa_final, fb_final) = build_fock_uhf(method, &b, &h, &g, &pa_final, &pb_final);
                let electronic_final = 0.5
                    * (pt_final.frobenius_dot(&h)
                        + pa_final.frobenius_dot(&fa_final)
                        + pb_final.frobenius_dot(&fb_final));
                return Ok(CndoIndoResult {
                    method,
                    density: pt_final.clone(),
                    fock: fa_final,
                    fock_beta: Some(fb_final),
                    mo_coeff: ca,
                    mo_coeff_beta: Some(cb),
                    mo_energies_ev: eps_a,
                    mo_energies_beta_ev: Some(eps_b),
                    n_occ: n_alpha,
                    n_alpha,
                    n_beta,
                    spin_density: Some(sub(&pa_final, &pb_final)),
                    unrestricted: true,
                    electronic_ev: electronic_final,
                    core_ev: core,
                    total_ev: electronic_final + core,
                    charges: charges(&b, &pt_final),
                    iterations: it,
                    converged: true,
                });
            }
        }
    }

    Err(XndoError::ScfNotConverged {
        iterations: opt.max_scf,
        error: last_error,
    })
}

/// Analytic Cartesian ground-state gradient in eV/Bohr.
pub(crate) fn analytic_ground_gradient(
    m: &Molecule,
    method: Method,
    opt: &CndoIndoOptions,
) -> Result<(CndoIndoResult, Vec<Vec3>)> {
    let result = run_cndo_indo(m, method, opt)?;
    let basis = Basis::build(m, method)?;
    let mut gradient = vec![Vec3::zero(); m.atoms.len()];
    let spin = result
        .spin_density
        .as_ref()
        .map(|s| spin_densities(&result.density, s));

    for ia in 0..m.atoms.len() {
        let ea = &basis.elements[ia];
        let oa = basis.offsets[ia];
        let na = basis.norb[ia];
        let pop_a = atom_population(&result.density, oa, na);
        for ib in (ia + 1)..m.atoms.len() {
            let eb = &basis.elements[ib];
            let ob = basis.offsets[ib];
            let nb = basis.norb[ib];
            let pop_b = atom_population(&result.density, ob, nb);
            let d = m.atoms[ib].position - m.atoms[ia].position;
            if d.norm() < 1.0e-12 {
                return Err(XndoError::InvalidInput(
                    "coincident atoms in CNDO/2 or INDO".into(),
                ));
            }
            let dv = [Dual::var(d.x, 0), Dual::var(d.y, 1), Dual::var(d.z, 2)];
            let r = (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]).sqrt();
            let gamma = gamma_ev(ea, eb, r);
            let overlap = overlap_block_g(ea, eb, dv);
            let beta = 0.5
                * bonding_adjust_k(ea, eb)
                * (ea.bonding_parameter_ev + eb.bonding_parameter_ev);
            let mut resonance = Dual::constant(0.0);
            for i in 0..na {
                for j in 0..nb {
                    resonance =
                        resonance + overlap[i][j] * (result.density[(oa + i, ob + j)] * beta);
                }
            }
            let exchange = match &spin {
                Some((alpha, beta_density)) => {
                    exchange_weight_uhf(alpha, beta_density, oa, na, ob, nb)
                }
                None => exchange_weight_rhf(&result.density, oa, na, ob, nb),
            };
            let core = Dual::constant(ea.core_charge * eb.core_charge * HARTREE_TO_EV) / r;
            let pair = pair_energy(
                gamma,
                core,
                ea.core_charge,
                eb.core_charge,
                pop_a,
                pop_b,
                exchange,
                resonance,
            );
            let dg = Vec3::new(pair.d[0], pair.d[1], pair.d[2]);
            gradient[ia] -= dg;
            gradient[ib] += dg;
        }
    }
    Ok((result, gradient))
}

/// Fully analytic Cartesian ground-state Hessian in eV/Bohr^2.
pub fn analytic_ground_hessian(
    m: &Molecule,
    method: Method,
    opt: &CndoIndoOptions,
) -> Result<(CndoIndoResult, Vec<Vec3>, Matrix)> {
    let result = run_cndo_indo(m, method, opt)?;
    let basis = Basis::build(m, method)?;
    let gamma0 = gamma_matrix(m, &basis)?;
    let ncoord = 3 * m.atoms.len();
    let nao = basis.nao();
    let zero = Matrix::zeros(nao, nao);
    let spin = result
        .spin_density
        .as_ref()
        .map(|s| spin_densities(&result.density, s));
    let mut skeleton_a = (0..ncoord)
        .map(|_| Matrix::zeros(nao, nao))
        .collect::<Vec<_>>();
    let mut skeleton_b = spin.as_ref().map(|_| {
        (0..ncoord)
            .map(|_| Matrix::zeros(nao, nao))
            .collect::<Vec<_>>()
    });

    struct PairData {
        ia: usize,
        ib: usize,
        gamma: Dual2,
        overlap: [[Dual2; 9]; 9],
        core: Dual2,
        beta: f64,
        exchange: f64,
        population_a: f64,
        population_b: f64,
    }
    let mut pairs = Vec::new();
    for ia in 0..m.atoms.len() {
        let ea = &basis.elements[ia];
        let oa = basis.offsets[ia];
        let na = basis.norb[ia];
        let pop_a = atom_population(&result.density, oa, na);
        for ib in (ia + 1)..m.atoms.len() {
            let eb = &basis.elements[ib];
            let ob = basis.offsets[ib];
            let nb = basis.norb[ib];
            let pop_b = atom_population(&result.density, ob, nb);
            let d = m.atoms[ib].position - m.atoms[ia].position;
            let dv = [Dual2::var(d.x, 0), Dual2::var(d.y, 1), Dual2::var(d.z, 2)];
            let r = (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]).sqrt();
            let gamma = gamma_ev(ea, eb, r);
            let overlap = overlap_block_g(ea, eb, dv);
            let beta = 0.5
                * bonding_adjust_k(ea, eb)
                * (ea.bonding_parameter_ev + eb.bonding_parameter_ev);
            let exchange = match &spin {
                Some((pa, pb)) => exchange_weight_uhf(pa, pb, oa, na, ob, nb),
                None => exchange_weight_rhf(&result.density, oa, na, ob, nb),
            };
            let core = Dual2::constant(ea.core_charge * eb.core_charge * HARTREE_TO_EV) / r;
            for axis in 0..3 {
                let dg = gamma.g[axis];
                for &(atom, sign) in &[(ia, -1.0), (ib, 1.0)] {
                    let coord = 3 * atom + axis;
                    for i in 0..na {
                        skeleton_a[coord][(oa + i, oa + i)] +=
                            sign * (-eb.core_charge + pop_b) * dg;
                        if let Some(sb) = skeleton_b.as_mut() {
                            sb[coord][(oa + i, oa + i)] += sign * (-eb.core_charge + pop_b) * dg;
                        }
                    }
                    for j in 0..nb {
                        skeleton_a[coord][(ob + j, ob + j)] +=
                            sign * (-ea.core_charge + pop_a) * dg;
                        if let Some(sb) = skeleton_b.as_mut() {
                            sb[coord][(ob + j, ob + j)] += sign * (-ea.core_charge + pop_a) * dg;
                        }
                    }
                    for i in 0..na {
                        for j in 0..nb {
                            let mu = oa + i;
                            let nu = ob + j;
                            let resonance = beta * overlap[i][j].g[axis];
                            let (xa, xb) = match &spin {
                                Some((pa, pb)) => {
                                    (resonance - pa[(mu, nu)] * dg, resonance - pb[(mu, nu)] * dg)
                                }
                                None => (resonance - 0.5 * result.density[(mu, nu)] * dg, 0.0),
                            };
                            skeleton_a[coord][(mu, nu)] += sign * xa;
                            skeleton_a[coord][(nu, mu)] += sign * xa;
                            if let Some(sb) = skeleton_b.as_mut() {
                                sb[coord][(mu, nu)] += sign * xb;
                                sb[coord][(nu, mu)] += sign * xb;
                            }
                        }
                    }
                }
            }
            pairs.push(PairData {
                ia,
                ib,
                gamma,
                overlap,
                core,
                beta,
                exchange,
                population_a: pop_a,
                population_b: pop_b,
            });
        }
    }

    let mut dp_total = Vec::with_capacity(ncoord);
    let mut dp_alpha = Vec::new();
    let mut dp_beta = Vec::new();
    if let Some((pa, pb)) = &spin {
        let cb = result
            .mo_coeff_beta
            .as_ref()
            .expect("UHF beta coefficients");
        let eps_b = result
            .mo_energies_beta_ev
            .as_ref()
            .expect("UHF beta energies");
        let sb = skeleton_b.as_ref().expect("UHF beta skeleton");
        let responses = solve_uhf_responses(
            &result.mo_coeff,
            cb,
            &result.mo_energies_ev,
            eps_b,
            result.n_alpha,
            result.n_beta,
            &skeleton_a,
            sb,
            |dpa, dpb| Ok(build_fock_uhf(method, &basis, &zero, &gamma0, dpa, dpb)),
        )?;
        for response in responses {
            let total = add(&response.density_alpha, &response.density_beta);
            dp_total.push(total);
            dp_alpha.push(response.density_alpha);
            dp_beta.push(response.density_beta);
        }
        let _ = (pa, pb);
    } else {
        for response in solve_rhf_responses(
            &result.mo_coeff,
            &result.mo_energies_ev,
            result.n_occ,
            &skeleton_a,
            |dp| Ok(build_fock_rhf(method, &basis, &zero, &gamma0, dp)),
        )? {
            dp_total.push(response.density);
        }
    }

    let mut gradient = vec![Vec3::zero(); m.atoms.len()];
    let mut hessian = Matrix::zeros(ncoord, ncoord);
    for pair in &pairs {
        let ea = &basis.elements[pair.ia];
        let eb = &basis.elements[pair.ib];
        let oa = basis.offsets[pair.ia];
        let ob = basis.offsets[pair.ib];
        let na = basis.norb[pair.ia];
        let nb = basis.norb[pair.ib];
        let mut resonance = Dual2::constant(0.0);
        for i in 0..na {
            for j in 0..nb {
                resonance =
                    resonance + pair.overlap[i][j] * (result.density[(oa + i, ob + j)] * pair.beta);
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
        let gpair = Vec3::new(fixed.g[0], fixed.g[1], fixed.g[2]);
        gradient[pair.ia] -= gpair;
        gradient[pair.ib] += gpair;
        for q in 0..3 {
            for r in 0..3 {
                for &(row_atom, row_sign) in &[(pair.ia, -1.0), (pair.ib, 1.0)] {
                    for &(col_atom, col_sign) in &[(pair.ia, -1.0), (pair.ib, 1.0)] {
                        hessian[(3 * row_atom + q, 3 * col_atom + r)] +=
                            row_sign * col_sign * fixed.h[q][r];
                    }
                }
            }
        }
        for coord in 0..ncoord {
            let dpt = &dp_total[coord];
            let dpop_a = atom_population(dpt, oa, na);
            let dpop_b = atom_population(dpt, ob, nb);
            let dexchange = if let Some((pa, pb)) = &spin {
                let mut value = 0.0;
                for i in 0..na {
                    for j in 0..nb {
                        let mu = oa + i;
                        let nu = ob + j;
                        value += 2.0
                            * (pa[(mu, nu)] * dp_alpha[coord][(mu, nu)]
                                + pb[(mu, nu)] * dp_beta[coord][(mu, nu)]);
                    }
                }
                value
            } else {
                let mut value = 0.0;
                for i in 0..na {
                    for j in 0..nb {
                        let mu = oa + i;
                        let nu = ob + j;
                        value += result.density[(mu, nu)] * dpt[(mu, nu)];
                    }
                }
                value
            };
            let dweight = -eb.core_charge * dpop_a - ea.core_charge * dpop_b
                + dpop_a * pair.population_b
                + pair.population_a * dpop_b
                - dexchange;
            for q in 0..3 {
                let mut response_term = pair.gamma.g[q] * dweight;
                for i in 0..na {
                    for j in 0..nb {
                        response_term +=
                            2.0 * dpt[(oa + i, ob + j)] * pair.beta * pair.overlap[i][j].g[q];
                    }
                }
                hessian[(3 * pair.ia + q, coord)] -= response_term;
                hessian[(3 * pair.ib + q, coord)] += response_term;
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
    Ok((result, gradient, hessian))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Vec3;
    use crate::system::Atom;

    fn mol(atoms: &[(u8, [f64; 3])]) -> Molecule {
        Molecule::new(
            atoms
                .iter()
                .map(|(z, r)| Atom {
                    z: *z,
                    position: Vec3::new(r[0], r[1], r[2]),
                })
                .collect(),
        )
    }

    #[test]
    fn donor_rows_and_orbital_counts_are_exactly_scoped() {
        assert_eq!(element(1).unwrap().n_orb, 1);
        assert_eq!(element(3).unwrap().n_orb, 4);
        assert_eq!(element(6).unwrap().n_orb, 4);
        assert_eq!(element(16).unwrap().n_orb, 9);
        assert!(validate_element(Method::Cndo2, 16).is_ok());
        assert!(validate_element(Method::Indo, 16).is_err());
        assert!(validate_element(Method::Cndo2, 3).is_ok());
        assert!(validate_element(Method::Indo, 3).is_ok());
    }

    #[test]
    fn gamma_is_symmetric_and_has_coulomb_tail() {
        let c = element(6).unwrap();
        let o = element(8).unwrap();
        for r in [0.5, 1.5, 4.0, 10.0] {
            let a = gamma_ev(&c, &o, r);
            let b = gamma_ev(&o, &c, r);
            assert!((a - b).abs() < 1.0e-11);
        }
        let far = gamma_ev(&c, &o, 20.0) / HARTREE_TO_EV;
        assert!((far - 1.0 / 20.0).abs() < 1.0e-5, "far gamma={far}");
    }

    #[test]
    fn indo_rhf_and_uhf_fock_reduce_to_same_equation_for_half_densities() {
        let m = mol(&[
            (8, [0.0, 0.0, 0.0]),
            (1, [0.0, 1.4, 0.0]),
            (1, [1.2, -0.4, 0.0]),
        ]);
        let b = Basis::build(&m, Method::Indo).unwrap();
        let s = build_overlap(&m, &b).unwrap();
        let g = gamma_matrix(&m, &b).unwrap();
        let h = build_hcore(Method::Indo, &m, &b, &s, &g);
        let guess = build_guess_fock(&m, &b, &s);
        let (_, c) = symmetric_eigen(&guess).unwrap();
        let p = c.leading_columns_gram(5, 2.0);
        let pa = c.leading_columns_gram(5, 1.0);
        let pb = pa.clone();
        let fr = build_fock_indo_rhf(&b, &h, &g, &p);
        let (fa, fb) = build_fock_indo_uhf(&b, &h, &g, &pa, &pb);
        assert!(fr.rms_difference(&fa) < 1.0e-11);
        assert!(fr.rms_difference(&fb) < 1.0e-11);
    }

    #[test]
    fn cndo2_rhf_and_uhf_fock_reduce_to_same_equation_for_half_densities() {
        let m = mol(&[
            (6, [0.0, 0.0, 0.0]),
            (1, [0.0, 0.0, 2.0]),
            (1, [1.7, 0.0, -0.6]),
            (1, [-0.8, 1.5, -0.6]),
            (1, [-0.8, -1.5, -0.6]),
        ]);
        let b = Basis::build(&m, Method::Cndo2).unwrap();
        let s = build_overlap(&m, &b).unwrap();
        let g = gamma_matrix(&m, &b).unwrap();
        let h = build_hcore(Method::Cndo2, &m, &b, &s, &g);
        let guess = build_guess_fock(&m, &b, &s);
        let (_, c) = symmetric_eigen(&guess).unwrap();
        let p = c.leading_columns_gram(4, 2.0);
        let pa = c.leading_columns_gram(4, 1.0);
        let pb = pa.clone();
        let fr = build_fock_cndo_rhf(&b, &h, &g, &p);
        let (fa, fb) = build_fock_cndo_uhf(&b, &h, &g, &pa, &pb);
        assert!(fr.rms_difference(&fa) < 1.0e-11);
        assert!(fr.rms_difference(&fb) < 1.0e-11);
    }

    #[test]
    fn cndo2_and_indo_have_independent_dispatchable_single_points() {
        let h2 = mol(&[(1, [0.0, 0.0, -0.7]), (1, [0.0, 0.0, 0.7])]);
        let cndo = run_cndo_indo(&h2, Method::Cndo2, &CndoIndoOptions::default()).unwrap();
        let indo = run_cndo_indo(&h2, Method::Indo, &CndoIndoOptions::default()).unwrap();
        assert!(cndo.total_ev.is_finite());
        assert!(indo.total_ev.is_finite());
        assert_eq!(cndo.method, Method::Cndo2);
        assert_eq!(indo.method, Method::Indo);
    }

    #[test]
    fn h2_matches_published_cndo2_indo_oracle() {
        // Independent published CINDO regression: H2 at R=0.74 Angstrom has
        // E_HF = -1.474625 Eh for both CNDO/2 and INDO.  Coordinates below
        // are in the internal Bohr unit.  The modest tolerance covers the
        // historical constants/rounding conventions while still catching
        // sign, gamma, resonance, core-repulsion, or dispatch errors.
        let r_bohr = 0.74 / crate::constants::BOHR_TO_ANGSTROM;
        let h2 = mol(&[
            (1, [0.0, 0.0, -0.5 * r_bohr]),
            (1, [0.0, 0.0, 0.5 * r_bohr]),
        ]);
        for method in [Method::Cndo2, Method::Indo] {
            let result = run_cndo_indo(&h2, method, &CndoIndoOptions::default()).unwrap();
            let energy_h = result.total_ev / HARTREE_TO_EV;
            assert!(
                (energy_h + 1.474625).abs() < 2.0e-4,
                "{method}: H2 oracle mismatch: {energy_h:.10} Eh"
            );
        }
    }

    #[test]
    fn open_shell_auto_selects_uhf() {
        let h = mol(&[(1, [0.0, 0.0, 0.0])]);
        let opt = CndoIndoOptions {
            multiplicity: 2,
            ..CndoIndoOptions::default()
        };
        let r = run_cndo_indo(&h, Method::Cndo2, &opt).unwrap();
        assert!(r.unrestricted);
        assert_eq!((r.n_alpha, r.n_beta), (1, 0));
        assert!(r.spin_density.is_some());
    }

    #[test]
    fn hydrogen_self_gamma_and_one_electron_energy_match_closed_form() {
        // For a normalized 1s Slater charge with zeta=1.2,
        // gamma_HH(0) = 5*zeta/8 = 0.75 Eh. This regression is independent
        // of the general Mulliken auxiliary-function implementation above.
        let hrow = element(1).unwrap();
        let gamma_h = gamma_ev(&hrow, &hrow, 0.0);
        let expected_gamma = 0.75 * HARTREE_TO_EV;
        assert!((gamma_h - expected_gamma).abs() < 1.0e-10);

        // A one-electron H atom has no pair interaction after same-spin
        // Coulomb/exchange cancellation. Both CNDO/2 and MolDS-labelled INDO
        // therefore reduce to h_1s = -I_s - gamma_HH/2.
        let h = mol(&[(1, [0.0, 0.0, 0.0])]);
        let opt = CndoIndoOptions {
            multiplicity: 2,
            reference: Reference::Uhf,
            ..CndoIndoOptions::default()
        };
        let expected = -7.176 - 0.5 * expected_gamma;
        for method in [Method::Cndo2, Method::Indo] {
            let r = run_cndo_indo(&h, method, &opt).unwrap();
            assert!(
                (r.total_ev - expected).abs() < 1.0e-9,
                "{method}: {} vs {expected}",
                r.total_ev
            );
        }
    }
}
