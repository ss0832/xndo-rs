// SPDX-License-Identifier: GPL-3.0-or-later

//! Dewar-Bingham-Lo MINDO/3 ground-state Hamiltonian.
//!
//! This is a native Rust RHF/UHF implementation of the 1975 MINDO/3
//! model.  It is deliberately separate from the NDDO/MNDO engine: MINDO/3
//! keeps the INDO-style one-centre electron repulsion structure, uses a
//! Mataga/Ohno-like atom-atom `gamma_AB` for all two-centre Coulomb terms, and
//! uses atom-pair resonance/core-core parameters.
//!
//! The numerical parameter set is audited against the public-domain MOPAC7
//! MINDO/3 tables and the original Bingham, Dewar & Lo formulation; see
//! `THIRD_PARTY_NOTICES.md` and `docs/parameter-provenance.md`.
//! Energies are evaluated internally in eV and distances in Angstrom.
//!
//! PROVENANCE: derived from MOPAC7 1.15,
//! public domain (no copyright asserted),
//! see third_party/mopac7/NOTICE.
//! UPSTREAM: fortran/block.f (MINDO/3 tables), analyt.f, delri.f, calpar.f, compfg.f.
//! MODIFIED for xndo-rs v0.3.0 on 2026-09-14:
//! native Rust RHF/UHF with analytic derivatives; evaluated with the
//! constants MOPAC7 itself uses.
//! Retained notices: NOTICE; per-file record: THIRD_PARTY_NOTICES.md.

// AO, atom, and Cartesian indices are intentionally explicit in the fixed-size
// matrix contractions below; iterator rewrites obscure the equations.
#![allow(clippy::needless_range_loop)]

use crate::constants::{ModelConstants, BOHR_TO_ANGSTROM};

use crate::dual::{Dual, Scalar};
use crate::dual2::Dual2;
use crate::error::{Result, XndoError};
use crate::linalg::{symmetric_eigen, Matrix};
use crate::math::Vec3;
use crate::orbitals::OrbitalEnergies;
use crate::scf::Reference;
use crate::scf_accel::{
    lowest_solution, open_shell_starts, sad_density_sp, uniform_valence_density, Accelerator,
    ScfAccelerator,
};
use crate::system::Molecule;
use crate::zdo_gradient::{
    atom_population, exchange_weight_rhf, exchange_weight_uhf, pair_energy, solve_rhf_responses,
    solve_uhf_responses, spin_densities,
};

/// MINDO/3 is evaluated with the constants MOPAC7 uses, which are the pre-2019
/// set ([`ModelConstants::HISTORICAL`]).
///
/// MOPAC7 1.15 writes them as literals rather than tabulating them:
/// `analyt.f:47` and `delri.f:22` have `A0 = 0.529167`; `calpar.f:129-165` and
/// `delri.f:23` use `27.21` for the Hartree; `compfg.f:147` converts the heat of
/// formation with `23.061`; and `analyt.f:147` uses `14.399` for the core-core
/// Coulomb term, which is [`COULOMB_EV_ANG`] below.
///
/// Same reasoning as MNDO/d (`tests/data/ORACLE_NOTES.md` item 13): the
/// parameters were fitted against these values, so they are part of the model
/// rather than a rounding choice.
const MOPAC7: ModelConstants = ModelConstants::HISTORICAL;

/// `e^2 * a0` in eV*Angstrom, MOPAC7's `14.399` (`analyt.f:147`).
///
/// The pre-2019 value; CODATA gives 14.399645478456. It is spelled out here
/// rather than taken from [`MOPAC7`] because MOPAC7 spells it out too.
const COULOMB_EV_ANG: f64 = 14.399;

#[derive(Clone, Copy, Debug)]
pub struct Mindo3Element {
    pub z: u8,
    pub n_s: u8,
    pub n_p: u8,
    pub n_orb: usize,
    pub core_charge: f64,
    pub uss: f64,
    pub upp: f64,
    pub vs: f64,
    pub vp: f64,
    pub zeta_s: f64,
    pub zeta_p: f64,
    pub gss: f64,
    pub gsp: f64,
    pub gpp: f64,
    pub gp2: f64,
    pub hsp: f64,
    /// The one-centre p-p' exchange integral `(pp'|pp')`, in eV.
    ///
    /// **Always exactly `(gpp - gp2) / 2`.** MOPAC7 does not store it: `fock1.f`
    /// writes `GPP(NI) - GP2(NI)` and `0.5*(GPP(NI) - GP2(NI))` into the Fock
    /// matrix directly, so the relation is structural rather than a coincidence
    /// of the published tables. It is a field here only because the rest of the
    /// one-centre set is, and `hp2_is_exactly_half_the_gpp_gp2_gap` holds it to
    /// the relation.
    ///
    /// Until v0.3.0 four of the nine values were stored rounded to two decimals
    /// -- N, Si, S and Cl, each 0.005 eV off -- which is a model change, not a
    /// formatting one. It cost 2.2e-2 eV per chlorine atom against MOPAC7:
    /// 8.6e-2 eV on CCl4, where every other molecule in the set agreed to about
    /// 8e-4. Chlorine showed it worst because `p^5` has the most p-p' pairs for
    /// the error to act on.
    pub hp2: f64,
    pub f03: f64,
    pub e_isol_ev: f64,
    pub e_heat_kcal: f64,
}

/// Canonical MINDO/3 element parameters used by the original main-group set.
// 3.14 is the published nitrogen hsp parameter, not an approximation to pi.
#[allow(clippy::approx_constant)]
pub fn element(z: u8) -> Result<Mindo3Element> {
    let e = match z {
        1 => Mindo3Element {
            z,
            n_s: 1,
            n_p: 1,
            n_orb: 1,
            core_charge: 1.0,
            uss: -12.505,
            upp: 0.0,
            vs: -13.605,
            vp: 0.0,
            zeta_s: 1.300000,
            zeta_p: 0.0,
            gss: 12.848,
            gsp: 0.0,
            gpp: 0.0,
            gp2: 0.0,
            hsp: 0.0,
            hp2: 0.0,
            f03: 12.848,
            e_isol_ev: -12.505,
            e_heat_kcal: 52.102,
        },
        5 => Mindo3Element {
            z,
            n_s: 2,
            n_p: 2,
            n_orb: 4,
            core_charge: 3.0,
            uss: -33.61,
            upp: -25.11,
            vs: -15.160,
            vp: -8.520,
            zeta_s: 1.211156,
            zeta_p: 0.972826,
            gss: 10.59,
            gsp: 9.56,
            gpp: 8.86,
            gp2: 7.86,
            hsp: 1.81,
            hp2: 0.50,
            f03: 8.958,
            e_isol_ev: -61.70,
            e_heat_kcal: 135.7,
        },
        6 => Mindo3Element {
            z,
            n_s: 2,
            n_p: 2,
            n_orb: 4,
            core_charge: 4.0,
            uss: -51.79,
            upp: -39.18,
            vs: -21.340,
            vp: -11.540,
            zeta_s: 1.739391,
            zeta_p: 1.709645,
            gss: 12.23,
            gsp: 11.47,
            gpp: 11.08,
            gp2: 9.84,
            hsp: 2.43,
            hp2: 0.62,
            f03: 10.833,
            e_isol_ev: -119.47,
            e_heat_kcal: 170.89,
        },
        7 => Mindo3Element {
            z,
            n_s: 2,
            n_p: 2,
            n_orb: 4,
            core_charge: 5.0,
            uss: -66.06,
            upp: -56.40,
            vs: -27.510,
            vp: -14.340,
            zeta_s: 2.704546,
            zeta_p: 1.870839,
            gss: 13.59,
            gsp: 12.66,
            gpp: 12.98,
            gp2: 11.59,
            hsp: 3.14,
            hp2: 0.695,
            f03: 12.377,
            e_isol_ev: -187.51,
            e_heat_kcal: 113.0,
        },
        8 => Mindo3Element {
            z,
            n_s: 2,
            n_p: 2,
            n_orb: 4,
            core_charge: 6.0,
            uss: -91.73,
            upp: -78.80,
            vs: -35.300,
            vp: -17.910,
            zeta_s: 3.640575,
            zeta_p: 2.168448,
            gss: 15.42,
            gsp: 14.48,
            gpp: 14.52,
            gp2: 12.98,
            hsp: 3.94,
            hp2: 0.77,
            f03: 13.985,
            e_isol_ev: -307.07,
            e_heat_kcal: 59.559,
        },
        9 => Mindo3Element {
            z,
            n_s: 2,
            n_p: 2,
            n_orb: 4,
            core_charge: 7.0,
            uss: -129.86,
            upp: -105.93,
            vs: -43.700,
            vp: -20.890,
            zeta_s: 3.111270,
            zeta_p: 1.419860,
            gss: 16.92,
            gsp: 17.25,
            gpp: 16.71,
            gp2: 14.91,
            hsp: 4.83,
            hp2: 0.90,
            f03: 16.250,
            e_isol_ev: -475.00,
            e_heat_kcal: 18.86,
        },
        14 => Mindo3Element {
            z,
            n_s: 3,
            n_p: 3,
            n_orb: 4,
            core_charge: 4.0,
            uss: -39.82,
            upp: -29.15,
            vs: -17.820,
            vp: -8.510,
            zeta_s: 1.629173,
            zeta_p: 1.381721,
            gss: 9.82,
            gsp: 8.36,
            gpp: 7.31,
            gp2: 6.54,
            hsp: 1.32,
            hp2: 0.385,
            f03: 7.57,
            e_isol_ev: -90.98,
            e_heat_kcal: 106.0,
        },
        15 => Mindo3Element {
            z,
            n_s: 3,
            n_p: 3,
            n_orb: 4,
            core_charge: 5.0,
            uss: -56.23,
            upp: -42.31,
            vs: -21.100,
            vp: -10.290,
            zeta_s: 1.926108,
            zeta_p: 1.590665,
            gss: 11.56,
            gsp: 10.08,
            gpp: 8.64,
            gp2: 7.68,
            hsp: 1.92,
            hp2: 0.48,
            f03: 9.00,
            e_isol_ev: -150.81,
            e_heat_kcal: 79.8,
        },
        16 => Mindo3Element {
            z,
            n_s: 3,
            n_p: 3,
            n_orb: 4,
            core_charge: 6.0,
            uss: -73.39,
            upp: -57.25,
            vs: -23.840,
            vp: -12.410,
            zeta_s: 1.719480,
            zeta_p: 1.403205,
            gss: 12.88,
            gsp: 11.26,
            gpp: 9.90,
            gp2: 8.83,
            hsp: 2.26,
            hp2: 0.535,
            f03: 10.20,
            e_isol_ev: -229.15,
            e_heat_kcal: 65.65,
        },
        17 => Mindo3Element {
            z,
            n_s: 3,
            n_p: 3,
            n_orb: 4,
            core_charge: 7.0,
            uss: -98.99,
            upp: -76.43,
            vs: -25.260,
            vp: -15.090,
            zeta_s: 3.430887,
            zeta_p: 1.627017,
            gss: 15.03,
            gsp: 13.16,
            gpp: 11.30,
            gp2: 9.97,
            hsp: 2.42,
            hp2: 0.665,
            f03: 11.73,
            e_isol_ev: -345.93,
            e_heat_kcal: 28.95,
        },
        _ => {
            return Err(XndoError::MissingParameter(format!(
                "MINDO/3 element Z={z}; canonical built-in set is H, B-F, Si-Cl"
            )))
        }
    };
    Ok(e)
}

/// MINDO/3 atom-pair resonance parameter beta_AB and core-core alpha_AB.
/// Values are symmetric and zero-valued pairs are intentionally rejected.
pub fn pair(za: u8, zb: u8) -> Result<(f64, f64)> {
    let (a, b) = if za >= zb { (za, zb) } else { (zb, za) };
    let v = match (a, b) {
        (1, 1) => (0.244770, 1.489450),
        (5, 1) => (0.185347, 2.090352),
        (5, 5) => (0.151324, 2.280544),
        (6, 1) => (0.315011, 1.475836),
        (6, 5) => (0.250031, 2.138291),
        (6, 6) => (0.419907, 1.371208),
        (7, 1) => (0.360776, 0.589380),
        (7, 5) => (0.310959, 1.909763),
        (7, 6) => (0.410886, 1.635259),
        (7, 7) => (0.377342, 2.029618),
        (8, 1) => (0.417759, 0.478901),
        (8, 5) => (0.349745, 2.484827),
        (8, 6) => (0.464514, 1.820975),
        (8, 7) => (0.458110, 1.873859),
        (8, 8) => (0.659407, 1.537190),
        (9, 1) => (0.195242, 3.771362),
        (9, 5) => (0.219591, 2.862183),
        (9, 6) => (0.247494, 2.725913),
        (9, 7) => (0.205347, 2.861667),
        (9, 8) => (0.334044, 2.266949),
        (9, 9) => (0.197464, 3.864997),
        (14, 1) => (0.289647, 0.940789),
        (14, 6) => (0.411377, 1.101382),
        (14, 14) => (0.291703, 0.918432),
        (15, 1) => (0.320118, 0.923170),
        (15, 6) => (0.457816, 1.029693),
        (15, 8) => (0.470000, 1.662500),
        (15, 9) => (0.300000, 1.750000),
        (15, 15) => (0.311790, 1.186652),
        (16, 1) => (0.220654, 1.700698),
        (16, 6) => (0.284620, 1.761370),
        (16, 7) => (0.313170, 1.878176),
        (16, 8) => (0.422890, 2.077240),
        (16, 16) => (0.202489, 1.751617),
        (17, 1) => (0.231653, 2.089404),
        (17, 6) => (0.315480, 1.676222),
        (17, 7) => (0.302298, 1.817064),
        (17, 15) => (0.277322, 1.543720),
        (17, 16) => (0.221764, 1.950318),
        (17, 17) => (0.258969, 1.792125),
        _ => {
            return Err(XndoError::MissingParameter(format!(
                "MINDO/3 atom-pair parameters Z={za},Z={zb}"
            )))
        }
    };
    Ok(v)
}

#[derive(Clone, Debug)]
pub struct Mindo3Options {
    pub charge: f64,
    pub multiplicity: usize,
    /// Restricted/unrestricted reference. `Auto` selects RHF for a closed shell
    /// and UHF otherwise. `Uhf` may also be forced for an even-electron singlet.
    pub reference: Reference,
    pub max_scf: usize,
    pub e_tol_ev: f64,
    pub p_tol: f64,
    pub damping: f64,
    /// SCF convergence accelerator. Defaults to A-DIIS then CDIIS, matching the
    /// NDDO driver; `damping` applies only when this is `ScfAccelerator::None`.
    pub accelerator: ScfAccelerator,
    /// Commutator norm below which A-DIIS hands over to CDIIS.
    pub adiis_switch: f64,
    /// Memory budget for the accelerator history, in MiB.
    pub scf_memory_mb: usize,
}
impl Default for Mindo3Options {
    fn default() -> Self {
        Self {
            charge: 0.0,
            multiplicity: 1,
            reference: Reference::Auto,
            max_scf: 250,
            e_tol_ev: 1.0e-8,
            p_tol: 1.0e-7,
            damping: 0.20,
            accelerator: ScfAccelerator::AdiisCdiis,
            adiis_switch: 0.1,
            scf_memory_mb: 512,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Mindo3Result {
    /// Total (alpha + beta) AO density.
    pub density: Matrix,
    /// Alpha Fock matrix for UHF, the ordinary Fock matrix for RHF.
    pub fock: Matrix,
    pub fock_beta: Option<Matrix>,
    /// Alpha MO coefficients for UHF, ordinary spatial MOs for RHF.
    pub mo_coeff: Matrix,
    pub mo_coeff_beta: Option<Matrix>,
    /// Alpha MO energies for UHF, ordinary spatial MO energies for RHF.
    pub mo_energies_ev: Vec<f64>,
    pub mo_energies_beta_ev: Option<Vec<f64>>,
    /// Legacy occupied count. For UHF this is the alpha occupied count.
    pub n_occ: usize,
    pub n_alpha: usize,
    pub n_beta: usize,
    pub spin_density: Option<Matrix>,
    pub unrestricted: bool,
    pub electronic_ev: f64,
    pub core_ev: f64,
    pub total_ev: f64,
    pub heat_of_formation_kcal: f64,
    pub charges: Vec<f64>,
    pub iterations: usize,
    pub converged: bool,
}

impl Mindo3Result {
    /// The orbital energies and the frontier quantities read off them.
    pub fn orbitals(&self) -> OrbitalEnergies {
        OrbitalEnergies::new(
            self.mo_energies_ev.clone(),
            self.n_alpha,
            self.mo_energies_beta_ev.clone(),
            self.n_beta,
        )
    }

    /// Highest occupied spin orbital over both channels, in eV.
    pub fn homo_ev(&self) -> Option<f64> {
        self.orbitals().homo_ev()
    }

    /// Lowest unoccupied spin orbital over both channels, in eV.
    pub fn lumo_ev(&self) -> Option<f64> {
        self.orbitals().lumo_ev()
    }

    /// HOMO-LUMO gap in eV; `None` unless both frontier orbitals exist.
    pub fn homo_lumo_gap_ev(&self) -> Option<f64> {
        self.orbitals().gap_ev()
    }
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
}
impl Basis {
    fn build(m: &Molecule) -> Result<Self> {
        let mut aos = Vec::new();
        let mut offsets = Vec::new();
        let mut norb = Vec::new();
        for (ia, a) in m.atoms.iter().enumerate() {
            let e = element(a.z)?;
            offsets.push(aos.len());
            norb.push(e.n_orb);
            for orb in 0..e.n_orb {
                aos.push(Ao { atom: ia, orb });
            }
        }
        Ok(Self { aos, offsets, norb })
    }
    fn nao(&self) -> usize {
        self.aos.len()
    }
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
fn overlap_block_g<S: Scalar>(ei: &Mindo3Element, ej: &Mindo3Element, d: [S; 3]) -> [[S; 4]; 4] {
    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if r.val() < 1e-14 {
        return [[S::cst(0.0); 4]; 4];
    }
    let inv_r = r.recip();
    let dir = [d[0] * inv_r, d[1] * inv_r, d[2] * inv_r];
    let (ea, eb, direction, swap) = if ei.n_s >= ej.n_s {
        (ei, ej, dir, false)
    } else {
        (ej, ei, [-dir[0], -dir[1], -dir[2]], true)
    };
    // Exponents in MOPAC7's Bohr, not CODATA's.
    //
    // The overlap depends on `zeta * R` and nothing else, so an exponent
    // tabulated against one Bohr radius and a distance measured in another is a
    // 1.93e-5 relative error in every exponential in the integral. MOPAC7 uses
    // a0 = 0.529167 A (`analyt.f:47`, `delri.f:22`), and these exponents were
    // fitted with it, so `zeta * k` with `k = a0_CODATA / a0_MOPAC7` is the
    // exact restatement of the same integral in the crate's internal Bohr --
    // not an approximation to it.
    //
    // This is the same defect ORACLE_NOTES item 13 fixed for MNDO/d, and
    // `params.rs` and `cndo_indo.rs` have applied it since; MINDO/3 was the
    // engine it was never carried to. Worth a few times 1e-4 eV per atom pair,
    // because the resonance energy it multiplies is several eV per bond.
    let k = MOPAC7.length_scale();
    let (s111, s211, s121, s221, s222) = crate::overlap_numeric::slater_locals_numeric(
        ea.n_s,
        ea.n_p.max(ea.n_s),
        ea.zeta_s * k,
        (ea.zeta_p * k).max(1e-12),
        eb.n_s,
        eb.n_p.max(eb.n_s),
        eb.zeta_s * k,
        (eb.zeta_p * k).max(1e-12),
        r,
    );
    let di = crate::overlap::build_di_g::<S>([s111, s211, s121, s221, s222], direction);
    if swap {
        transpose4(di)
    } else {
        di
    }
}

fn overlap_block(ei: &Mindo3Element, ri: Vec3, ej: &Mindo3Element, rj: Vec3) -> [[f64; 4]; 4] {
    let d = rj - ri;
    overlap_block_g(ei, ej, [d.x, d.y, d.z])
}
fn build_overlap(m: &Molecule, b: &Basis) -> Result<Matrix> {
    let mut s = Matrix::identity(b.nao());
    for ia in 0..m.atoms.len() {
        let ei = element(m.atoms[ia].z)?;
        for ja in ia + 1..m.atoms.len() {
            let ej = element(m.atoms[ja].z)?;
            let block = overlap_block(&ei, m.atoms[ia].position, &ej, m.atoms[ja].position);
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
fn gamma<S: Scalar>(ea: &Mindo3Element, eb: &Mindo3Element, r_ang: S) -> S {
    let ra = COULOMB_EV_ANG / ea.f03;
    let rb = COULOMB_EV_ANG / eb.f03;
    (r_ang * r_ang + 0.25 * (ra + rb) * (ra + rb))
        .sqrt()
        .recip()
        * COULOMB_EV_ANG
}

fn one_center_jk(e: &Mindo3Element) -> ([[f64; 4]; 4], [[f64; 4]; 4]) {
    let mut j = [[0.0; 4]; 4];
    let mut k = [[0.0; 4]; 4];
    j[0][0] = e.gss;
    if e.n_orb > 1 {
        for p in 1..4 {
            j[0][p] = e.gsp;
            j[p][0] = e.gsp;
            k[0][p] = e.hsp;
            k[p][0] = e.hsp;
            j[p][p] = e.gpp;
        }
        for p in 1..4 {
            for q in 1..4 {
                if p != q {
                    j[p][q] = e.gp2;
                    k[p][q] = e.hp2;
                }
            }
        }
    }
    (j, k)
}
fn build_hcore(m: &Molecule, b: &Basis, s: &Matrix) -> Result<Matrix> {
    let n = b.nao();
    let mut h = Matrix::zeros(n, n);
    let mut gatom = vec![vec![0.0; m.atoms.len()]; m.atoms.len()];
    for ia in 0..m.atoms.len() {
        let ei = element(m.atoms[ia].z)?;
        for ja in ia + 1..m.atoms.len() {
            let ej = element(m.atoms[ja].z)?;
            let r = (m.atoms[ia].position - m.atoms[ja].position).norm() * BOHR_TO_ANGSTROM;
            let g = gamma(&ei, &ej, r);
            gatom[ia][ja] = g;
            gatom[ja][ia] = g;
        }
    }
    for mu in 0..n {
        let a = b.aos[mu];
        let ea = element(m.atoms[a.atom].z)?;
        let u = if a.orb == 0 { ea.uss } else { ea.upp };
        let vnuc = (0..m.atoms.len())
            .filter(|&ib| ib != a.atom)
            .map(|ib| element(m.atoms[ib].z).map(|eb| gatom[a.atom][ib] * eb.core_charge))
            .collect::<Result<Vec<_>>>()?
            .iter()
            .sum::<f64>();
        h[(mu, mu)] = u - vnuc;
        for nu in 0..mu {
            let c = b.aos[nu];
            if a.atom == c.atom {
                continue;
            }
            let ec = element(m.atoms[c.atom].z)?;
            let ip_a = if a.orb == 0 { ea.vs } else { ea.vp };
            let ip_c = if c.orb == 0 { ec.vs } else { ec.vp };
            let (beta, _) = pair(ea.z, ec.z)?;
            let v = s[(mu, nu)] * (ip_a + ip_c) * beta;
            h[(mu, nu)] = v;
            h[(nu, mu)] = v;
        }
    }
    Ok(h)
}
fn build_jk_potentials(m: &Molecule, b: &Basis, p: &Matrix) -> Result<(Matrix, Matrix)> {
    let n = b.nao();
    let mut vj = Matrix::zeros(n, n);
    let mut vk = Matrix::zeros(n, n);
    // One-centre contributions. The construction is linear in P; this makes
    // the same J/K builder usable in both RHF and UHF.
    for ia in 0..m.atoms.len() {
        let e = element(m.atoms[ia].z)?;
        let (ji, ki) = one_center_jk(&e);
        let o = b.offsets[ia];
        let nb = b.norb[ia];
        for i in 0..nb {
            let mu = o + i;
            let mut jd = 0.0;
            let mut kd = 0.0;
            for j in 0..nb {
                let nu = o + j;
                jd += ji[i][j] * p[(nu, nu)];
                kd += ki[i][j] * p[(nu, nu)];
            }
            vj[(mu, mu)] = jd;
            vk[(mu, mu)] = kd;
        }
        for i in 0..nb {
            for j in 0..nb {
                let mu = o + i;
                let nu = o + j;
                vj[(mu, nu)] += 2.0 * ki[i][j] * p[(mu, nu)];
                vk[(mu, nu)] += (ji[i][j] + ki[i][j]) * p[(mu, nu)];
            }
        }
    }
    // Two-centre Mataga/Ohno gamma contributions.
    let mut gatom = vec![vec![0.0; m.atoms.len()]; m.atoms.len()];
    for ia in 0..m.atoms.len() {
        let ei = element(m.atoms[ia].z)?;
        for ja in ia + 1..m.atoms.len() {
            let ej = element(m.atoms[ja].z)?;
            let r = (m.atoms[ia].position - m.atoms[ja].position).norm() * BOHR_TO_ANGSTROM;
            let g = gamma(&ei, &ej, r);
            gatom[ia][ja] = g;
            gatom[ja][ia] = g;
        }
    }
    let mut pop = vec![0.0; m.atoms.len()];
    for ia in 0..m.atoms.len() {
        let o = b.offsets[ia];
        for i in 0..b.norb[ia] {
            pop[ia] += p[(o + i, o + i)];
        }
    }
    for ia in 0..m.atoms.len() {
        let shift = (0..m.atoms.len())
            .filter(|&j| j != ia)
            .map(|j| gatom[ia][j] * pop[j])
            .sum::<f64>();
        let o = b.offsets[ia];
        for i in 0..b.norb[ia] {
            vj[(o + i, o + i)] += shift;
        }
        for ja in 0..m.atoms.len() {
            if ia == ja {
                continue;
            }
            let g = gatom[ia][ja];
            let oi = b.offsets[ia];
            let oj = b.offsets[ja];
            for i in 0..b.norb[ia] {
                for j in 0..b.norb[ja] {
                    vk[(oi + i, oj + j)] += g * p[(oi + i, oj + j)];
                }
            }
        }
    }
    Ok((vj, vk))
}

fn build_fock_rhf(m: &Molecule, b: &Basis, h: &Matrix, p: &Matrix) -> Result<Matrix> {
    let (vj, vk) = build_jk_potentials(m, b, p)?;
    let mut f = h.clone();
    for idx in 0..f.as_slice().len() {
        f.as_mut_slice()[idx] = h.as_slice()[idx] + vj.as_slice()[idx] - 0.5 * vk.as_slice()[idx];
    }
    Ok(f)
}

fn build_fock_uhf(
    m: &Molecule,
    b: &Basis,
    h: &Matrix,
    pa: &Matrix,
    pb: &Matrix,
) -> Result<(Matrix, Matrix)> {
    let pt = add(pa, pb);
    let (jtot, _) = build_jk_potentials(m, b, &pt)?;
    let (_, ka) = build_jk_potentials(m, b, pa)?;
    let (_, kb) = build_jk_potentials(m, b, pb)?;
    let mut fa = h.clone();
    let mut fb = h.clone();
    for idx in 0..h.as_slice().len() {
        fa.as_mut_slice()[idx] = h.as_slice()[idx] + jtot.as_slice()[idx] - ka.as_slice()[idx];
        fb.as_mut_slice()[idx] = h.as_slice()[idx] + jtot.as_slice()[idx] - kb.as_slice()[idx];
    }
    Ok((fa, fb))
}

fn core_energy(m: &Molecule) -> Result<f64> {
    let mut e = 0.0;
    for ia in 0..m.atoms.len() {
        let a = element(m.atoms[ia].z)?;
        for ja in ia + 1..m.atoms.len() {
            let b = element(m.atoms[ja].z)?;
            let r = (m.atoms[ia].position - m.atoms[ja].position).norm() * BOHR_TO_ANGSTROM;
            if r < 1e-12 {
                return Err(XndoError::InvalidInput(
                    "coincident atoms in MINDO/3".into(),
                ));
            }
            let g = gamma(&a, &b, r);
            let (_, alpha) = pair(a.z, b.z)?;
            let hn = (a.z == 1 && (b.z == 7 || b.z == 8)) || (b.z == 1 && (a.z == 7 || a.z == 8));
            let scale = if hn {
                alpha * (-r).exp()
            } else {
                (-alpha * r).exp()
            };
            e += a.core_charge * b.core_charge * (g + scale * (COULOMB_EV_ANG / r - g));
        }
    }
    Ok(e)
}
fn electron_count(m: &Molecule, charge: f64) -> Result<usize> {
    let q = charge.round();
    if (charge - q).abs() > 1e-8 {
        return Err(XndoError::InvalidInput(
            "MINDO/3 charge must be integral".into(),
        ));
    }
    let mut ne = -(q as i64);
    for a in &m.atoms {
        ne += element(a.z)?.core_charge.round() as i64;
    }
    if ne < 0 {
        return Err(XndoError::InvalidInput(
            "negative MINDO/3 electron count".into(),
        ));
    }
    Ok(ne as usize)
}
/// Superposition of atomic densities.
///
/// This used to spread an atom's valence electrons uniformly over all of its
/// AOs, which puts s and p on an equal footing they do not have. The shared
/// `sad_density_sp` fills s first and spreads the remainder over p, which is
/// the usual SAD and is what every engine now starts from.
fn initial_density(m: &Molecule, b: &Basis) -> Result<Matrix> {
    Ok(sad_density_sp(b.nao(), &valence_shells(m, b)?))
}

/// `(first AO index, AO count, core charge)` per atom, the shape the shared
/// guess builders take.
fn valence_shells(m: &Molecule, b: &Basis) -> Result<Vec<(usize, usize, f64)>> {
    m.atoms
        .iter()
        .enumerate()
        .map(|(ia, a)| -> Result<(usize, usize, f64)> {
            let e = element(a.z)?;
            Ok((b.offsets[ia], e.n_orb, e.core_charge))
        })
        .collect()
}

fn electron_partition(ne: usize, multiplicity: usize) -> Result<(usize, usize)> {
    if multiplicity == 0 {
        return Err(XndoError::InvalidInput(
            "MINDO/3 multiplicity must be >= 1".into(),
        ));
    }
    let unpaired = multiplicity - 1;
    if unpaired > ne || !(ne - unpaired).is_multiple_of(2) {
        return Err(XndoError::InvalidInput(format!(
            "MINDO/3 electron count {ne} is incompatible with multiplicity {multiplicity}"
        )));
    }
    let nb = (ne - unpaired) / 2;
    Ok((nb + unpaired, nb))
}

fn damp_density(raw: &Matrix, old: &Matrix, damping: f64) -> Matrix {
    let mut next = raw.clone();
    if damping > 0.0 {
        for i in 0..next.as_slice().len() {
            next.as_mut_slice()[i] =
                (1.0 - damping) * raw.as_slice()[i] + damping * old.as_slice()[i];
        }
    }
    next
}

fn sub(a: &Matrix, b: &Matrix) -> Matrix {
    let mut c = Matrix::zeros(a.rows, a.cols);
    for i in 0..c.as_slice().len() {
        c.as_mut_slice()[i] = a.as_slice()[i] - b.as_slice()[i];
    }
    c
}
fn charges(m: &Molecule, b: &Basis, p: &Matrix) -> Result<Vec<f64>> {
    let mut q = Vec::new();
    for ia in 0..m.atoms.len() {
        let e = element(m.atoms[ia].z)?;
        let mut pop = 0.0;
        for i in 0..b.norb[ia] {
            pop += p[(b.offsets[ia] + i, b.offsets[ia] + i)];
        }
        q.push(e.core_charge - pop);
    }
    Ok(q)
}
fn reference_kcal(m: &Molecule) -> Result<f64> {
    let mut hf = 0.0;
    let mut eisol = 0.0;
    for a in &m.atoms {
        let e = element(a.z)?;
        hf += e.e_heat_kcal;
        eisol += e.e_isol_ev;
    }
    Ok(hf - eisol * MOPAC7.ev_to_kcal)
}
fn add(a: &Matrix, b: &Matrix) -> Matrix {
    let mut c = Matrix::zeros(a.rows, a.cols);
    for i in 0..c.as_slice().len() {
        c.as_mut_slice()[i] = a.as_slice()[i] + b.as_slice()[i];
    }
    c
}

/// RHF/UHF MINDO/3 single point.
///
/// The unrestricted branch uses the standard spin-separated Fock equations
/// `F^alpha = H + J[P] - K[P^alpha]` and
/// `F^beta = H + J[P] - K[P^beta]`, with `P=P^alpha+P^beta`.
pub fn run_mindo3(m: &Molecule, opt: &Mindo3Options) -> Result<Mindo3Result> {
    let b = Basis::build(m)?;
    let ne = electron_count(m, opt.charge)?;
    let (n_alpha, n_beta) = electron_partition(ne, opt.multiplicity)?;
    if n_alpha > b.nao() || n_beta > b.nao() {
        return Err(XndoError::InvalidInput(
            "MINDO/3 electron count exceeds basis capacity".into(),
        ));
    }
    let reference = match opt.reference {
        Reference::Auto => {
            if n_alpha == n_beta {
                Reference::Rhf
            } else {
                Reference::Uhf
            }
        }
        Reference::Rhf => Reference::Rhf,
        Reference::Uhf => Reference::Uhf,
    };
    if reference == Reference::Rhf && n_alpha != n_beta {
        return Err(XndoError::InvalidInput(format!(
            "MINDO/3 RHF is incompatible with {} alpha and {} beta electrons; use reference=UHF/Auto",
            n_alpha, n_beta
        )));
    }

    let s = build_overlap(m, &b)?;
    let h = build_hcore(m, &b, &s)?;
    let core = core_energy(m)?;
    let damping = opt.damping.clamp(0.0, 0.95);
    let mut last_e = f64::INFINITY;
    let mut last_err = f64::INFINITY;

    if reference == Reference::Rhf {
        let nocc = n_alpha;
        let mut p = initial_density(m, &b)?;
        let mut accel = Accelerator::new(
            b.nao(),
            opt.accelerator,
            opt.adiis_switch,
            opt.scf_memory_mb,
        );
        for it in 1..=opt.max_scf {
            let f = build_fock_rhf(m, &b, &h, &p)?;
            let f_use = accel.step(f.clone(), &p);
            let (_eps, c) = symmetric_eigen(&f_use)?;
            let raw = c.leading_columns_gram(nocc, 2.0);
            let pn = if accel.is_active() {
                raw
            } else {
                damp_density(&raw, &p, damping)
            };
            let eel = 0.5 * pn.frobenius_dot(&add(&h, &f));
            let perr = pn.rms_difference(&p);
            let eerr = (eel - last_e).abs();
            p = pn;
            last_e = eel;
            last_err = perr;
            if perr < opt.p_tol && eerr < opt.e_tol_ev {
                let ff = build_fock_rhf(m, &b, &h, &p)?;
                let (ee, cc) = symmetric_eigen(&ff)?;
                p = cc.leading_columns_gram(nocc, 2.0);
                let elec = 0.5 * p.frobenius_dot(&add(&h, &ff));
                let total = elec + core;
                return Ok(Mindo3Result {
                    density: p.clone(),
                    fock: ff,
                    fock_beta: None,
                    mo_coeff: cc,
                    mo_coeff_beta: None,
                    mo_energies_ev: ee,
                    mo_energies_beta_ev: None,
                    n_occ: nocc,
                    n_alpha,
                    n_beta,
                    spin_density: None,
                    unrestricted: false,
                    electronic_ev: elec,
                    core_ev: core,
                    total_ev: total,
                    heat_of_formation_kcal: total * MOPAC7.ev_to_kcal + reference_kcal(m)?,
                    charges: charges(m, &b, &p)?,
                    iterations: it,
                    converged: true,
                });
            }
        }
    } else {
        // One SCF per candidate start, keeping the lowest: which orbital the
        // unpaired electron lands in is settled by the first Fock matrix, and a
        // guess can get that ordering wrong in a way no amount of converging
        // undoes. See `scf_accel::open_shell_starts`.
        let solve = |start: (Matrix, Matrix)| -> Result<Mindo3Result> {
            let (mut pa, mut pb) = start;
            let mut last_e = f64::INFINITY;
            let mut last_err = f64::INFINITY;
            // One accelerator for both spin channels. They are *not* independent --
            // each Fock matrix is built from the total density, so it depends on the
            // other channel's -- and the textbook UHF-DIIS error vector is the
            // stacked pair extrapolated with a single set of coefficients.
            let mut accel = Accelerator::new_uhf(
                b.nao(),
                opt.accelerator,
                opt.adiis_switch,
                opt.scf_memory_mb,
            );
            for it in 1..=opt.max_scf {
                let (fa, fb) = build_fock_uhf(m, &b, &h, &pa, &pb)?;
                // The accelerator consumes the Fock matrix it extrapolates, but the
                // energy below has to be evaluated with the *unextrapolated* one --
                // the extrapolated matrix is a fit, not `F[P]`, and pairing it with
                // a density would report an energy the model never produced.
                let (fa_use, fb_use) = accel.step_uhf(fa.clone(), fb.clone(), &pa, &pb);
                let (_ea, ca) = symmetric_eigen(&fa_use)?;
                let (_eb, cb) = symmetric_eigen(&fb_use)?;
                let pa_raw = ca.leading_columns_gram(n_alpha, 1.0);
                let pb_raw = cb.leading_columns_gram(n_beta, 1.0);
                let (pa_next, pb_next) = if accel.is_active() {
                    (pa_raw, pb_raw)
                } else {
                    (
                        damp_density(&pa_raw, &pa, damping),
                        damp_density(&pb_raw, &pb, damping),
                    )
                };
                let pt_next = add(&pa_next, &pb_next);
                let eel = 0.5
                    * (pt_next.frobenius_dot(&h)
                        + pa_next.frobenius_dot(&fa)
                        + pb_next.frobenius_dot(&fb));
                let perr = pa_next.rms_difference(&pa).max(pb_next.rms_difference(&pb));
                let eerr = (eel - last_e).abs();
                pa = pa_next;
                pb = pb_next;
                last_e = eel;
                last_err = perr;
                if perr < opt.p_tol && eerr < opt.e_tol_ev {
                    let (fa2, fb2) = build_fock_uhf(m, &b, &h, &pa, &pb)?;
                    let (ea2, ca2) = symmetric_eigen(&fa2)?;
                    let (eb2, cb2) = symmetric_eigen(&fb2)?;
                    pa = ca2.leading_columns_gram(n_alpha, 1.0);
                    pb = cb2.leading_columns_gram(n_beta, 1.0);
                    let pt = add(&pa, &pb);
                    let (fa3, fb3) = build_fock_uhf(m, &b, &h, &pa, &pb)?;
                    let elec = 0.5
                        * (pt.frobenius_dot(&h) + pa.frobenius_dot(&fa3) + pb.frobenius_dot(&fb3));
                    let total = elec + core;
                    return Ok(Mindo3Result {
                        density: pt.clone(),
                        fock: fa3,
                        fock_beta: Some(fb3),
                        mo_coeff: ca2,
                        mo_coeff_beta: Some(cb2),
                        mo_energies_ev: ea2,
                        mo_energies_beta_ev: Some(eb2),
                        n_occ: n_alpha,
                        n_alpha,
                        n_beta,
                        spin_density: Some(sub(&pa, &pb)),
                        unrestricted: true,
                        electronic_ev: elec,
                        core_ev: core,
                        total_ev: total,
                        heat_of_formation_kcal: total * MOPAC7.ev_to_kcal + reference_kcal(m)?,
                        charges: charges(m, &b, &pt)?,
                        iterations: it,
                        converged: true,
                    });
                }
            }
            Err(XndoError::ScfNotConverged {
                iterations: opt.max_scf,
                error: last_err,
            })
        };
        let sad = initial_density(m, &b)?;
        let uniform = uniform_valence_density(b.nao(), &valence_shells(m, &b)?);
        return lowest_solution(
            open_shell_starts(&sad, &uniform, n_alpha, n_beta),
            solve,
            |r| r.total_ev,
        );
    }
    Err(XndoError::ScfNotConverged {
        iterations: opt.max_scf,
        error: last_err,
    })
}

/// Analytic Cartesian ground-state gradient in eV/Bohr.
pub(crate) fn analytic_ground_gradient(
    m: &Molecule,
    opt: &Mindo3Options,
) -> Result<(Mindo3Result, Vec<Vec3>)> {
    let result = run_mindo3(m, opt)?;
    let basis = Basis::build(m)?;
    let elements = m
        .atoms
        .iter()
        .map(|atom| element(atom.z))
        .collect::<Result<Vec<_>>>()?;
    let spin = result
        .spin_density
        .as_ref()
        .map(|s| spin_densities(&result.density, s));
    let mut gradient = vec![Vec3::zero(); m.atoms.len()];

    for ia in 0..m.atoms.len() {
        let ea = &elements[ia];
        let oa = basis.offsets[ia];
        let na = basis.norb[ia];
        let pop_a = atom_population(&result.density, oa, na);
        for ib in (ia + 1)..m.atoms.len() {
            let eb = &elements[ib];
            let ob = basis.offsets[ib];
            let nb = basis.norb[ib];
            let pop_b = atom_population(&result.density, ob, nb);
            let d = m.atoms[ib].position - m.atoms[ia].position;
            if d.norm() < 1.0e-12 {
                return Err(XndoError::InvalidInput(
                    "coincident atoms in MINDO/3".into(),
                ));
            }
            let dv = [Dual::var(d.x, 0), Dual::var(d.y, 1), Dual::var(d.z, 2)];
            let r_bohr = (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]).sqrt();
            let r_ang = r_bohr * BOHR_TO_ANGSTROM;
            let gamma_ab = gamma(ea, eb, r_ang);
            let overlap = overlap_block_g(ea, eb, dv);
            let (beta, alpha) = pair(ea.z, eb.z)?;
            let mut resonance = Dual::constant(0.0);
            for i in 0..na {
                let ip_a = if i == 0 { ea.vs } else { ea.vp };
                for j in 0..nb {
                    let ip_b = if j == 0 { eb.vs } else { eb.vp };
                    resonance = resonance
                        + overlap[i][j] * (result.density[(oa + i, ob + j)] * beta * (ip_a + ip_b));
                }
            }
            let exchange = match &spin {
                Some((alpha_density, beta_density)) => {
                    exchange_weight_uhf(alpha_density, beta_density, oa, na, ob, nb)
                }
                None => exchange_weight_rhf(&result.density, oa, na, ob, nb),
            };
            let hydrogen_hetero =
                (ea.z == 1 && (eb.z == 7 || eb.z == 8)) || (eb.z == 1 && (ea.z == 7 || ea.z == 8));
            let scale = if hydrogen_hetero {
                (-r_ang).exp() * alpha
            } else {
                (r_ang * (-alpha)).exp()
            };
            let core = (gamma_ab + scale * (Dual::constant(COULOMB_EV_ANG) / r_ang - gamma_ab))
                * (ea.core_charge * eb.core_charge);
            let pair_derivative = pair_energy(
                gamma_ab,
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

/// Fully analytic Cartesian RHF/UHF Hessian in eV/Bohr^2.
pub fn analytic_ground_hessian(
    m: &Molecule,
    opt: &Mindo3Options,
) -> Result<(Mindo3Result, Vec<Vec3>, Matrix)> {
    let result = run_mindo3(m, opt)?;
    let basis = Basis::build(m)?;
    let elements = m
        .atoms
        .iter()
        .map(|atom| element(atom.z))
        .collect::<Result<Vec<_>>>()?;
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
        overlap: [[Dual2; 4]; 4],
        core: Dual2,
        beta: f64,
        exchange: f64,
        population_a: f64,
        population_b: f64,
    }
    let mut pairs = Vec::new();
    for ia in 0..m.atoms.len() {
        let ea = &elements[ia];
        let oa = basis.offsets[ia];
        let na = basis.norb[ia];
        let pop_a = atom_population(&result.density, oa, na);
        for ib in (ia + 1)..m.atoms.len() {
            let eb = &elements[ib];
            let ob = basis.offsets[ib];
            let nb = basis.norb[ib];
            let pop_b = atom_population(&result.density, ob, nb);
            let d = m.atoms[ib].position - m.atoms[ia].position;
            let dv = [Dual2::var(d.x, 0), Dual2::var(d.y, 1), Dual2::var(d.z, 2)];
            let r_bohr = (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]).sqrt();
            let r_ang = r_bohr * BOHR_TO_ANGSTROM;
            let gamma = gamma(ea, eb, r_ang);
            let overlap = overlap_block_g(ea, eb, dv);
            let (beta, alpha) = pair(ea.z, eb.z)?;
            let exchange = match &spin {
                Some((pa, pb)) => exchange_weight_uhf(pa, pb, oa, na, ob, nb),
                None => exchange_weight_rhf(&result.density, oa, na, ob, nb),
            };
            let hydrogen_hetero =
                (ea.z == 1 && (eb.z == 7 || eb.z == 8)) || (eb.z == 1 && (ea.z == 7 || ea.z == 8));
            let scale = if hydrogen_hetero {
                (-r_ang).exp() * alpha
            } else {
                (r_ang * (-alpha)).exp()
            };
            let core = (gamma + scale * (Dual2::constant(COULOMB_EV_ANG) / r_ang - gamma))
                * (ea.core_charge * eb.core_charge);
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
                        let ip_a = if i == 0 { ea.vs } else { ea.vp };
                        for j in 0..nb {
                            let ip_b = if j == 0 { eb.vs } else { eb.vp };
                            let mu = oa + i;
                            let nu = ob + j;
                            let resonance = beta * (ip_a + ip_b) * overlap[i][j].g[axis];
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
    if spin.is_some() {
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
            |dpa, dpb| build_fock_uhf(m, &basis, &zero, dpa, dpb),
        )?;
        for response in responses {
            dp_total.push(add(&response.density_alpha, &response.density_beta));
            dp_alpha.push(response.density_alpha);
            dp_beta.push(response.density_beta);
        }
    } else {
        for response in solve_rhf_responses(
            &result.mo_coeff,
            &result.mo_energies_ev,
            result.n_occ,
            &skeleton_a,
            |dp| build_fock_rhf(m, &basis, &zero, dp),
        )? {
            dp_total.push(response.density);
        }
    }

    let mut gradient = vec![Vec3::zero(); m.atoms.len()];
    let mut hessian = Matrix::zeros(ncoord, ncoord);
    for pair_data in &pairs {
        let ea = &elements[pair_data.ia];
        let eb = &elements[pair_data.ib];
        let oa = basis.offsets[pair_data.ia];
        let ob = basis.offsets[pair_data.ib];
        let na = basis.norb[pair_data.ia];
        let nb = basis.norb[pair_data.ib];
        let mut resonance = Dual2::constant(0.0);
        for i in 0..na {
            let ip_a = if i == 0 { ea.vs } else { ea.vp };
            for j in 0..nb {
                let ip_b = if j == 0 { eb.vs } else { eb.vp };
                resonance = resonance
                    + pair_data.overlap[i][j]
                        * (result.density[(oa + i, ob + j)] * pair_data.beta * (ip_a + ip_b));
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
                + dpop_a * pair_data.population_b
                + pair_data.population_a * dpop_b
                - dexchange;
            for q in 0..3 {
                let mut response_term = pair_data.gamma.g[q] * dweight;
                for i in 0..na {
                    let ip_a = if i == 0 { ea.vs } else { ea.vp };
                    for j in 0..nb {
                        let ip_b = if j == 0 { eb.vs } else { eb.vp };
                        response_term += 2.0
                            * dpt[(oa + i, ob + j)]
                            * pair_data.beta
                            * (ip_a + ip_b)
                            * pair_data.overlap[i][j].g[q];
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
    Ok((result, gradient, hessian))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hp2_is_exactly_half_the_gpp_gp2_gap() {
        // MOPAC7 never stores this integral: `fock1.f` writes `GPP - GP2` and
        // `0.5*(GPP - GP2)` straight into the Fock matrix, so the relation is
        // part of the model and not a property of how the tables were printed.
        // Four values used to be stored rounded to two decimals, which moved
        // CCl4 by 8.6e-2 eV against MOPAC7.
        //
        // 1e-12 rather than exact equality. Not a hedge: none of these decimal
        // literals is representable in binary, so `8.86 - 7.86` is
        // 0.9999999999999991 and the difference from the literal `0.5` is real
        // arithmetic, about 4e-17. The bound sits five orders below the
        // smallest thing that could be a transcription error (a unit in the
        // third decimal, 1e-3) and four above the rounding, so it separates the
        // two without admitting anything in between.
        for z in [5, 6, 7, 8, 9, 14, 15, 16, 17] {
            let e = element(z).unwrap();
            let derived = 0.5 * (e.gpp - e.gp2);
            assert!(
                (e.hp2 - derived).abs() < 1.0e-12,
                "Z={z}: hp2 {} is not (gpp {} - gp2 {})/2 = {derived}",
                e.hp2,
                e.gpp,
                e.gp2,
            );
        }
    }

    #[test]
    fn canonical_elements_exist() {
        for z in [1, 5, 6, 7, 8, 9, 14, 15, 16, 17] {
            assert!(element(z).is_ok());
        }
    }
    #[test]
    fn canonical_pairs_are_symmetric() {
        for &(a, b) in &[(1, 6), (6, 8), (14, 6), (15, 17), (16, 17)] {
            assert_eq!(pair(a, b).ok(), pair(b, a).ok());
        }
    }
    #[test]
    fn mopac7_chlorine_pair_tail_is_not_shifted() {
        assert_eq!(pair(15, 17).unwrap(), (0.277322, 1.543720));
        assert_eq!(pair(16, 17).unwrap(), (0.221764, 1.950318));
        assert_eq!(pair(17, 17).unwrap(), (0.258969, 1.792125));
    }
    #[test]
    fn electron_partition_handles_doublet_and_triplet() {
        assert_eq!(electron_partition(7, 2).unwrap(), (4, 3));
        assert_eq!(electron_partition(8, 3).unwrap(), (5, 3));
        assert!(electron_partition(8, 2).is_err());
    }
}
