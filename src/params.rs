// SPDX-License-Identifier: GPL-3.0-or-later

//! MNDO element/pair parameters and the derived NDDO multipole quantities.
//!
//! Parses the embedded MNDO parameter tables (extracted from MOPAC v23.2.5,
//! see [`crate::data_tables`]) and, per element, precomputes the
//! dipole/quadrupole charge separations `dd`/`qq` and the KlopmanOhno
//! additive terms `rho0/rho1/rho2` used by the two-center two-electron
//! integrals. The closed forms and the `rho1`/`rho2` secant solves follow
//! MOPAC `calpar.F90`; the d-shell multipole extension (MNDO-d `ddpo`/`poij`)
//! is added by the d-orbital milestone.

use crate::constants::{EV_TO_KCAL, HARTREE_TO_EV};
use crate::data_tables::{self, CsvTable};
use crate::error::{Result, XndoError};
use crate::method::Method;
use std::collections::HashMap;

/// Per-element MNDO parameters plus derived NDDO quantities.
#[derive(Clone, Debug)]
pub struct NddoElement {
    pub z: u8,
    /// Principal quantum number of the valence s shell (`npq_s`); the sp
    /// overlap fast path uses this as "the" quantum number, which is exact for
    /// every element with `npq_s == npq_p` (all sp-gate elements).
    pub n: u8,
    pub n_s: u8,
    pub n_p: u8,
    pub n_d: u8,
    pub u_ss: f64,
    pub u_pp: f64,
    pub u_dd: f64,
    pub zeta_s: f64,
    pub zeta_p: f64,
    pub zeta_d: f64,
    pub beta_s: f64,
    pub beta_p: f64,
    pub beta_d: f64,
    pub g_ss: f64,
    pub g_sp: f64,
    pub g_pp: f64,
    pub g_p2: f64,
    pub h_sp: f64,
    /// Internal one-center exponents for transition metals (`zsn/zpn/zdn`).
    pub zsn: f64,
    pub zpn: f64,
    pub zdn: f64,
    /// Explicit SlaterCondon overrides for d elements.
    pub f0sd: f64,
    pub g2sd: f64,
    /// Per-element core-core exponent `alp6` (only some elements; MNDO uses
    /// pairwise `alpb`/`xfac` for parametrized pairs).
    pub alpha: f64,
    /// Per-element core additive term `poc` (`poc_6`, MOPAC `pocord`): overrides
    /// the core monopole `po(9)` when nonzero. Only a few elements (e.g. Sc, Fe,
    /// Ni) define it; used in electroncore attraction and corecore repulsion.
    pub poc: f64,
    /// MNDO per-element core-core Gaussian corrections `(K, L, M)` from
    /// `gues61/62/63` (nonzero triples only; at most 2 for MNDO).
    pub gauss: Vec<(f64, f64, f64)>,
    /// Number of valence AOs: 1 (s), 4 (s,p) or 9 (s,p,d).
    pub n_orb: usize,
    /// Core charge (valence-electron count, MOPAC `tore`).
    pub core_charge: f64,
    /// Initial shell occupancies (`ios`/`iop`/`iod`), used by the SAD guess
    /// and the isolated-atom energy.
    pub occ_s: f64,
    pub occ_p: f64,
    pub occ_d: f64,
    /// MOPAC `main_group` flag: one-center integrals from `Gss...Hsp` directly.
    pub main_group: bool,
    /// MOPAC `ndelec`: d electrons folded into the core in `Eisol` bookkeeping.
    pub ndelec: i32,
    /// Experimental atomic heat of formation (eV).
    pub eheat_ev: f64,
    /// Isolated-atom electronic energy (eV).
    pub e_isol: f64,
    /// Atomic mass (amu).
    pub mass: f64,
    // Derived NDDO multipole terms (Bohr).
    pub dd: f64,
    pub qq: f64,
    pub rho0: f64,
    pub rho1: f64,
    pub rho2: f64,
    /// MNDO-d charge separations `ddp(1..=6)` (index 0 unused; 2=sp dipole,
    /// 3=pp quadrupole, 4=sd quadrupole, 5=pd dipole, 6=dd quadrupole), Bohr.
    pub ddp: [f64; 7],
    /// MNDO-d additive Klopman terms `po(1..=9)` (1=ss, 2=sp, 3=pp, 4=sd, 5=pd,
    /// 6=dd, 7=pp-monopole, 8=dd-monopole, 9=core), Bohr.
    pub po: [f64; 10],
    /// One-center spd two-electron integrals (d elements only).
    pub onecenter: Option<crate::onecenter::OneCenterSpd>,
}

impl NddoElement {
    pub fn has_p(&self) -> bool {
        self.n_orb >= 4
    }
    pub fn has_d(&self) -> bool {
        self.n_orb >= 9
    }
}

/// Pairwise MNDO core-core parameters (`alpb` in A1-like units, unitless `xfac`).
#[derive(Clone, Copy, Debug)]
pub struct PairParams {
    pub alpha: f64,
    pub x: f64,
}

#[derive(Clone, Debug, Default)]
pub struct NddoParameters {
    /// Model carried by this parameter set.
    pub method: Method,
    pub elements: HashMap<u8, NddoElement>,
    /// Symmetric key `(max(zi,zj), min(zi,zj))`.
    pub pair: HashMap<(u8, u8), PairParams>,
}

impl NddoParameters {
    /// Load the OpenMOPAC-v23.2.5 MNDO parameter set.
    pub fn mndo() -> Result<Self> {
        Self::from_tables(
            Method::Mndo,
            "MNDO",
            data_tables::MNDO_PARAM_CSV,
            data_tables::MNDO_PAIR_CSV,
        )
    }

    /// Load the OpenMOPAC-v23.2.5 MNDO/d parameter set.
    pub fn mndod() -> Result<Self> {
        Self::from_tables(
            Method::MndoD,
            "MNDO/d",
            data_tables::MNDOD_PARAM_CSV,
            data_tables::MNDOD_PAIR_CSV,
        )
    }

    /// Load parameters for an NDDO method. ZINDO/S uses a distinct INDO engine.
    pub fn for_method(method: Method) -> Result<Self> {
        match method {
            Method::Mndo => Self::mndo(),
            Method::MndoD => Self::mndod(),
            other => Err(XndoError::InvalidInput(format!(
                "{} is not handled by the MNDO/NDDO parameter container: {}",
                other,
                other.implementation_note()
            ))),
        }
    }

    /// Return a process-wide parsed table for repeated calculations.
    pub fn cached_for_method(method: Method) -> Result<&'static Self> {
        use std::sync::OnceLock;
        static MNDO: OnceLock<NddoParameters> = OnceLock::new();
        static MNDOD: OnceLock<NddoParameters> = OnceLock::new();
        match method {
            Method::Mndo => Ok(MNDO
                .get_or_init(|| Self::mndo().expect("embedded MNDO parameter tables are valid"))),
            Method::MndoD => Ok(MNDOD.get_or_init(|| {
                Self::mndod().expect("embedded MNDO/d parameter tables are valid")
            })),
            other => Err(XndoError::InvalidInput(format!(
                "{} is not handled by the MNDO/NDDO parameter container: {}",
                other,
                other.implementation_note()
            ))),
        }
    }

    pub fn element(&self, z: u8) -> Result<&NddoElement> {
        self.elements.get(&z).ok_or(XndoError::MissingElement(z))
    }

    /// Pairwise `alpb`/`xfac`; error if MOPAC has no parameters for the pair.
    pub fn pair(&self, zi: u8, zj: u8) -> Result<&PairParams> {
        let key = (zi.max(zj), zi.min(zj));
        self.pair.get(&key).ok_or_else(|| {
            XndoError::MissingParameter(format!(
                "{} diatomic parameters for Z={zi},Z={zj}",
                self.method
            ))
        })
    }

    fn from_tables(
        method: Method,
        model_name: &str,
        param_csv: &str,
        pair_csv: &str,
    ) -> Result<Self> {
        let edata = data_tables::element_data();

        // --- per-element parameters ---
        let t = CsvTable::parse(param_csv).ok_or_else(|| {
            XndoError::InvalidInput(format!("empty {model_name} parameter table"))
        })?;
        let col = |n: &str| -> Result<usize> {
            t.col(n).ok_or_else(|| {
                XndoError::MissingParameter(format!("{model_name} parameter column {n}"))
            })
        };
        let c_z = col("z")?;
        let names = [
            "uss", "upp", "udd", "zs", "zp", "zd", "betas", "betap", "betad", "gss", "gsp", "gpp",
            "gp2", "hsp", "zsn", "zpn", "zdn", "f0sd", "g2sd", "alp", "poc",
        ];
        let mut ci = HashMap::new();
        for n in names {
            ci.insert(n, col(n)?);
        }
        let gcols = [
            (col("g1_k")?, col("g1_l")?, col("g1_m")?),
            (col("g2_k")?, col("g2_l")?, col("g2_m")?),
            (col("g3_k")?, col("g3_l")?, col("g3_m")?),
            (col("g4_k")?, col("g4_l")?, col("g4_m")?),
        ];

        let mut elements = HashMap::new();
        for row in &t.rows {
            let z = t.f64_at(row, c_z) as u8;
            // OpenMOPAC records with Z >= 100 are not chemical elements and
            // are intentionally not loaded.
            if z >= 100 {
                continue;
            }
            let get = |n: &str| t.f64_at(row, ci[n]);
            let (u_ss, beta_s) = (get("uss"), get("betas"));
            if u_ss == 0.0 && beta_s == 0.0 {
                continue;
            }
            if z as usize >= edata.len() || z == 0 {
                continue;
            }
            let ed = edata[z as usize];
            let (zeta_p, zeta_d) = (get("zp"), get("zd"));
            let n_orb = if zeta_d > 0.0 {
                9
            } else if zeta_p > 0.0 {
                4
            } else {
                1
            };
            let mut gauss = Vec::new();
            for (kk, ll, mm) in gcols {
                let (k, l, m) = (t.f64_at(row, kk), t.f64_at(row, ll), t.f64_at(row, mm));
                if k != 0.0 || l != 0.0 {
                    gauss.push((k, l, m));
                }
            }

            let n = ed.npq_s;
            let zeta_s = get("zs");
            let (mut g_ss, mut g_sp, mut g_pp, mut g_p2) =
                (get("gss"), get("gsp"), get("gpp"), get("gp2"));
            // MOPAC calpar floors: zeta_p >= 0.3, hsp >= 1e-7, hpp >= 0.1.
            let mut h_sp = get("hsp").max(1.0e-7);
            // Transition metals (main_group = false): the one-center
            // Coulomb/exchange integrals are NOT the fitted values but are
            // recomputed from the internal exponents `zsn`/`zpn` via the
            // SlaterCondon radial integrals (MOPAC `sp_two_electron`). This is
            // essential  using the fitted values makes the d shell far too
            // deep and the SCF collapses to an unphysical charge-transfer state.
            let (zsn, zpn) = (get("zsn"), get("zpn"));
            if !ed.main_group && zsn > 1.0e-4 && zpn > 1.0e-4 {
                use crate::onecenter::slater_rsc;
                let ns = ed.npq_s as i32;
                let np = ed.npq_p as i32;
                g_ss = slater_rsc(0, ns, zsn, ns, zsn, ns, zsn, ns, zsn);
                g_sp = slater_rsc(0, ns, zsn, ns, zsn, np, zpn, np, zpn);
                h_sp = (slater_rsc(1, ns, zsn, np, zpn, ns, zsn, np, zpn) / 3.0).max(1.0e-7);
                let r033 = slater_rsc(0, np, zpn, np, zpn, np, zpn, np, zpn);
                let r233 = slater_rsc(2, np, zpn, np, zpn, np, zpn, np, zpn);
                g_pp = r033 + 0.16 * r233;
                g_p2 = r033 - 0.08 * r233;
            }

            let (dd, qq, rho1, rho2) = if n_orb >= 4 && zeta_p > 0.0 {
                let zp = zeta_p.max(0.3);
                // MOPAC uses the periodic-table row PQN (`nspqn`) for the multipole charge
                // separations, not the per-element Slater `n` (see `multipole_pqn`).
                let qn = multipole_pqn(z);
                let dd = dd_charge_sep(qn, zeta_s, zp);
                let qq = qq_charge_sep(qn, zp);
                let hpp = (0.5 * (g_pp - g_p2)).max(0.1);
                let rho1 = additive_rho1(h_sp, dd);
                let rho2 = additive_rho2(hpp, qq);
                (dd, qq, rho1, rho2)
            } else {
                (0.0, 0.0, 0.0, 0.0)
            };
            let rho0 = if g_ss > 0.0 {
                0.5 * HARTREE_TO_EV / g_ss
            } else {
                0.0
            };

            // Isolated-atom electronic energy: average-of-configuration
            // coefficients in the MOPAC shell occupancies (`calpar.F90:75-113`).
            // Exact for main-group elements; the d-block adds `eiscor` terms in
            // the d milestone. Hydrogen is special-cased (eisol = uss).
            let (ios, iop) = (ed.occ_s, ed.occ_p);
            let e_isol = if z == 1 {
                u_ss
            } else {
                u_ss * ios
                    + get("upp") * iop
                    + get("udd") * ed.occ_d
                    + g_ss * gssc(ios)
                    + g_pp * gppc(iop)
                    + g_sp * gspc(ios, iop)
                    + g_p2 * gp2c(iop)
                    + h_sp * hspc(ios, iop)
            };

            let mut elem = NddoElement {
                z,
                n,
                n_s: ed.npq_s,
                n_p: ed.npq_p,
                n_d: ed.npq_d,
                u_ss,
                u_pp: get("upp"),
                u_dd: get("udd"),
                zeta_s,
                zeta_p,
                zeta_d,
                beta_s,
                beta_p: get("betap"),
                beta_d: get("betad"),
                g_ss,
                g_sp,
                g_pp,
                g_p2,
                h_sp,
                zsn: get("zsn"),
                zpn: get("zpn"),
                zdn: get("zdn"),
                f0sd: get("f0sd"),
                g2sd: get("g2sd"),
                alpha: get("alp"),
                poc: get("poc"),
                gauss,
                n_orb,
                core_charge: ed.tore,
                occ_s: ed.occ_s,
                occ_p: ed.occ_p,
                occ_d: ed.occ_d,
                main_group: ed.main_group,
                ndelec: ed.ndelec,
                eheat_ev: ed.eheat_kcal / EV_TO_KCAL,
                e_isol,
                mass: ed.mass,
                dd,
                qq,
                rho0,
                rho1,
                rho2,
                ddp: [0.0; 7],
                po: [0.0; 10],
                onecenter: None,
            };
            // Derived MNDO-d multipole terms. For a d element the one-center
            // spd integrals are built first (they feed the d-multipole solver
            // and add the d-shell isolated-atom energy); the sp additive terms
            // then follow MOPAC `inid`'s main-group overwrite.
            derive_multipoles(&mut elem);
            elements.insert(z, elem);
        }
        if elements.is_empty() {
            return Err(XndoError::InvalidInput(format!(
                "no {model_name} elements parsed from parameter table"
            )));
        }

        // --- pairwise alpb/xfac ---
        let tp = CsvTable::parse(pair_csv).ok_or_else(|| {
            XndoError::InvalidInput(format!("empty {model_name} pair-parameter table"))
        })?;
        let (c_zi, c_zj, c_a, c_x) = (
            tp.col("zi")
                .ok_or_else(|| XndoError::MissingParameter("zi".into()))?,
            tp.col("zj")
                .ok_or_else(|| XndoError::MissingParameter("zj".into()))?,
            tp.col("alpb")
                .ok_or_else(|| XndoError::MissingParameter("alpb".into()))?,
            tp.col("xfac")
                .ok_or_else(|| XndoError::MissingParameter("xfac".into()))?,
        );
        let mut pair = HashMap::new();
        for row in &tp.rows {
            let zi = tp.f64_at(row, c_zi) as u8;
            let zj = tp.f64_at(row, c_zj) as u8;
            let key = (zi.max(zj), zi.min(zj));
            pair.insert(
                key,
                PairParams {
                    alpha: tp.f64_at(row, c_a),
                    x: tp.f64_at(row, c_x),
                },
            );
        }

        Ok(Self {
            method,
            elements,
            pair,
        })
    }
}

// Average-of-configuration coefficients for the isolated-atom electronic
// energy, in the integer shell occupancies `ios`/`iop` (MOPAC `calpar.F90:75-97`,
// transcribed exactly so the transition-metal path stays correct):
//   gssc = max(ios-1, 0)
//   gspc = ios*iop
//   l    = min(iop, 6-iop)
//   gp2c = (iop*(iop-1))/2 [integer] + (l*(l-1))/4
//   gppc = -(l*(l-1))/4
//   hspc = -iop*ios/2

fn gssc(ios: f64) -> f64 {
    (ios - 1.0).max(0.0)
}
fn gspc(ios: f64, iop: f64) -> f64 {
    ios * iop
}
fn hspc(ios: f64, iop: f64) -> f64 {
    -iop * ios * 0.5
}
fn gp2c(iop: f64) -> f64 {
    let k = iop as i64;
    let l = k.min(6 - k);
    ((k * (k - 1)) / 2) as f64 + (l * (l - 1)) as f64 / 4.0
}
fn gppc(iop: f64) -> f64 {
    let k = iop as i64;
    let l = k.min(6 - k);
    -((l * (l - 1)) as f64) / 4.0
}

/// Multipole principal quantum number `nspqn` used by the two-center multipole charge
/// separations (MOPAC `calpar.F90:44`: `data nspqn/2*1, 8*2, 8*3, 18*4, 18*5, 32*6, 21*0/`).
///
/// This is the **periodic-table row** number, NOT the per-element Slater principal quantum
/// number `n`: MNDO assigns some elements a Slater `n` one higher than their row (e.g. Ne uses a
/// 3s Slater orbital, Ar a 4s), but the `dd`/`qq` multipole formulas always use the row number.
/// For most main-group elements the two coincide (C/N/O = 2), so only the outliers (chiefly the
/// noble gases) were affected.
pub fn multipole_pqn(z: u8) -> f64 {
    match z {
        1..=2 => 1.0,
        3..=10 => 2.0,
        11..=18 => 3.0,
        19..=36 => 4.0,
        37..=54 => 5.0,
        55..=86 => 6.0,
        _ => 0.0,
    }
}

/// Dipole charge separation `dd` (Bohr). MOPAC `calpar.F90`. `qn` is the multipole principal
/// quantum number ([`multipole_pqn`]).
pub fn dd_charge_sep(qn: f64, zs: f64, zp: f64) -> f64 {
    (2.0 * qn + 1.0) * (4.0 * zs * zp).powf(qn + 0.5)
        / (zs + zp).powf(2.0 * qn + 2.0)
        / 3.0_f64.sqrt()
}

/// Quadrupole charge separation `qq` (Bohr). MOPAC `calpar.F90`. `qn` is the multipole principal
/// quantum number ([`multipole_pqn`]).
pub fn qq_charge_sep(qn: f64, zp: f64) -> f64 {
    ((4.0 * qn * qn + 6.0 * qn + 2.0) / 20.0).sqrt() / zp
}

/// Additive term `rho1` reproducing the one-center dipole integral `H_sp` (Bohr).
///
/// Solves `H_sp(au) = 12 d  12 / (4 D12 + 1/d2)` for `d`, returning `rho1 = 0.5/d`.
///
/// **Exactly 5 secant iterations**, mirroring MOPAC `calpar.F90:159-171` (`jmax = 5`)  the
/// solver is *not* run to convergence. For most elements it has effectively converged by 5
/// iterations, but for elements whose target integral was floored (`hpp  0.1` when `gpp < gp2`,
/// e.g. every noble gas) the 5-iteration truncation differs from the converged root, and matching
/// MOPAC's truncation is required for bit-agreement.
pub fn additive_rho1(hsp_ev: f64, dd: f64) -> f64 {
    // MOPAC floors hsp to 1e-7 eV before the secant (calpar.F90:111).
    let hsp = hsp_ev.max(1.0e-7) / HARTREE_TO_EV;
    let g = |d: f64| 0.5 * d - 0.5 / (4.0 * dd * dd + 1.0 / (d * d)).sqrt();
    let gdd1 = (hsp / (dd * dd)).powf(1.0 / 3.0);
    let (mut d1, mut d2) = (gdd1, gdd1 + 0.04);
    for _ in 0..5 {
        let df = d2 - d1;
        let (h1, h2) = (g(d1), g(d2));
        if (h2 - h1).abs() < 1.0e-25 {
            break;
        }
        let d3 = d1 + df * (hsp - h1) / (h2 - h1);
        d1 = d2;
        d2 = d3;
    }
    0.5 / d2
}

/// Additive term `rho2` reproducing the one-center quadrupole integral `H_pp` (Bohr).
///
/// Solves `H_pp(au) = 14 q  12/(4 D22 + 1/q2) + 14/(8 D22 + 1/q2)` for `q`, returning
/// `rho2 = 0.5/q`. **Exactly 5 secant iterations** (MOPAC `calpar.F90:172-185`, `jmax = 5`);
/// see [`additive_rho1`].
pub fn additive_rho2(hpp_ev: f64, qq: f64) -> f64 {
    let hpp = hpp_ev / HARTREE_TO_EV;
    let g = |q: f64| {
        0.25 * q - 0.5 / (4.0 * qq * qq + 1.0 / (q * q)).sqrt()
            + 0.25 / (8.0 * qq * qq + 1.0 / (q * q)).sqrt()
    };
    // Start point: p4 = 2^4 = 16, matching MOPAC `gqq = (p4*hpp/(ev*48*qq^4))^0.2`.
    let gqq = (16.0 * hpp / (48.0 * qq.powi(4))).powf(0.2);
    let (mut q1, mut q2) = (gqq, gqq + 0.04);
    for _ in 0..5 {
        let qf = q2 - q1;
        let (h1, h2) = (g(q1), g(q2));
        if (h2 - h1).abs() < 1.0e-25 {
            break;
        }
        let q3 = q1 + qf * (hpp - h1) / (h2 - h1);
        q1 = q2;
        q2 = q3;
    }
    0.5 / q2
}

/// Fill the MNDO-d charge separations `ddp` and additive terms `po` on an
/// element (MOPAC `aijm`/`ddpo`/`poij` + the main-group overwrite in `inid`).
/// For a d element this also builds and attaches [`crate::onecenter::OneCenterSpd`]
/// and folds its `eisol_d` into `e_isol`.
fn derive_multipoles(elem: &mut NddoElement) {
    // One-center spd integrals (d elements): needed for the d additive-term
    // targets (`repd`) and the isolated-atom d-shell energy.
    if elem.has_d() {
        let oc = crate::onecenter::OneCenterSpd::build(elem);
        elem.e_isol += oc.eisol_d;
        // aij (unnormalized-exponent multipole normalizations, MOPAC `aijm`).
        let aij = aijm(elem);
        // po(1)/ss monopole.
        if elem.g_ss > 0.1 {
            elem.po[1] = poij(0, 1.0, elem.g_ss);
        }
        // sp dipole, pp quadrupole.
        let d_sp = aij[2] / 12.0_f64.sqrt();
        elem.ddp[2] = d_sp;
        elem.po[2] = poij(1, d_sp, elem.h_sp);
        elem.po[7] = elem.po[1];
        let d_pp = (aij[3] * 0.1).sqrt();
        elem.ddp[3] = d_pp;
        elem.po[3] = poij(2, d_pp, 0.5 * (elem.g_pp - elem.g_p2));
        // d multipoles.
        let da = (1.0_f64 / 60.0).sqrt();
        let d_sd = (aij[4] * da).sqrt();
        elem.ddp[4] = d_sd;
        elem.po[4] = poij(2, d_sd, oc.repd[19]);
        let d_pd = aij[5] / 20.0_f64.sqrt();
        elem.ddp[5] = d_pd;
        elem.po[5] = poij(1, d_pd, oc.repd[23] - 1.8 * oc.repd[35]);
        let fg_dd = 0.2 * (oc.repd[29] + 2.0 * oc.repd[30] + 2.0 * oc.repd[31]);
        elem.po[8] = if fg_dd > 1e-5 {
            poij(0, 1.0, fg_dd)
        } else {
            1e5
        };
        let d_dd = (aij[6] / 14.0).sqrt();
        elem.ddp[6] = d_dd;
        elem.po[6] = poij(2, d_dd, oc.repd[44] - (20.0 / 35.0) * oc.repd[52]);
        elem.po[9] = elem.po[1];
        elem.onecenter = Some(oc);
    } else if elem.n_orb >= 4 {
        // Main-group sp element: MOPAC `inid` overwrites po(1,2,3,7) and ddp(2,3)
        // from the calpar am/ad/aq/dd/qq path.
        elem.po[1] = elem.rho0;
        elem.po[2] = elem.rho1;
        elem.po[3] = elem.rho2;
        elem.po[7] = elem.rho0;
        elem.po[9] = elem.rho0;
        elem.ddp[2] = elem.dd;
        elem.ddp[3] = elem.qq * 2.0_f64.sqrt();
    } else {
        // Hydrogen / s-only: monopole core term only.
        elem.po[1] = elem.rho0;
        elem.po[7] = elem.rho0;
        elem.po[9] = elem.rho0;
    }
    // Core additive-term override (MOPAC `inid`): `po(9) = pocord` when the
    // element defines one (Sc, Fe, Ni, ...); otherwise `po(9) = po(1)`. `po(9)`
    // enters BOTH the electroncore attraction and the corecore repulsion,
    // where the two nearly cancel. Our two-center d path seeds the s/p
    // electron-core block from the `rho0`-based s/p path, so `po(9)` currently
    // reaches only the d-core terms; applying `poc` there alone breaks the
    // cancellation and worsens the (already physical) Sc/Fe/Ni energies. Until
    // the s/p electron-core seed is made `po(9)`-consistent, keep `po(9)=po(1)`
    // (self-consistent with the core-core), leaving a ~1.5 kcal/mol residual on
    // those three elements. `elem.poc` is retained for that future rework.
    let _ = elem.poc;
}

/// MOPAC `aijm`/`aijl`: multipole normalization factors from the *valence*
/// (unnormalized) Slater exponents. Returns `aij[1..=6]` (index 0 unused).
fn aijm(elem: &NddoElement) -> [f64; 7] {
    let mut aij = [0.0f64; 7];
    let (z1, z2, z3) = (elem.zeta_s, elem.zeta_p, elem.zeta_d);
    let nsp = elem.n_s as i32;
    if elem.z < 3 || z1 * z2 < 0.01 {
        return aij;
    }
    aij[2] = aijl(z1, z2, nsp, nsp, 1);
    aij[3] = aijl(z2, z2, nsp, nsp, 2);
    if elem.has_d() {
        let nd = elem.n_d as i32;
        aij[4] = aijl(z1, z3, nsp, nd, 2);
        aij[5] = aijl(z2, z3, nsp, nd, 1);
        aij[6] = aijl(z3, z3, nd, nd, 2);
    }
    aij
}

/// MOPAC `aijl` (mndod.F90:2146). `fx(i) = (i-1)!` (1-based factorial).
fn aijl(z1: f64, z2: f64, n1: i32, n2: i32, l: i32) -> f64 {
    fn fac(n: i32) -> f64 {
        // n! as a plain factorial (n >= 0).
        (1..=n).map(|k| k as f64).product::<f64>().max(1.0)
    }
    let zz = z1 + z2 + 1e-20;
    // MOPAC uses fx(n1+n2+l+1) = (n1+n2+l)! and fx(2*n1+1)=(2*n1)!, fx(2*n2+1)=(2*n2)!.
    fac(n1 + n2 + l) / (fac(2 * n1) * fac(2 * n2)).sqrt()
        * (2.0 * z1 / zz).powi(n1)
        * (2.0 * z1 / zz).sqrt()
        * (2.0 * z2 / zz).powi(n2)
        * (2.0 * z2 / zz).sqrt()
        * 2.0_f64.powi(l)
        / zz.powi(l)
}

/// MOPAC `poij` (mndod.F90:165): additive Klopman term (Bohr) reproducing the
/// one-center multipole integral `fg` (eV). Golden-section minimization for
/// `l = 1, 2`; closed form for the monopole `l = 0`.
fn poij(l: i32, d: f64, fg: f64) -> f64 {
    if l == 0 {
        return 0.5 * HARTREE_TO_EV / fg;
    }
    let ev4 = HARTREE_TO_EV / 4.0;
    let ev8 = HARTREE_TO_EV / 8.0;
    let dsq = d * d;
    let (mut a1, mut a2) = (0.1, 5.0);
    let (mut f1, mut f2) = (0.0, 0.0);
    for _ in 0..100 {
        let delta = a2 - a1;
        if delta < 1e-8 {
            break;
        }
        let y1 = a1 + delta * 0.382;
        let y2 = a1 + delta * 0.618;
        match l {
            1 => {
                f1 = (ev4 * (1.0 / y1 - 1.0 / (y1 * y1 + dsq).sqrt()) - fg).powi(2);
                f2 = (ev4 * (1.0 / y2 - 1.0 / (y2 * y2 + dsq).sqrt()) - fg).powi(2);
            }
            _ => {
                f1 = (ev8
                    * (1.0 / y1 - 2.0 / (y1 * y1 + dsq * 0.5).sqrt()
                        + 1.0 / (y1 * y1 + dsq).sqrt())
                    - fg)
                    .powi(2);
                f2 = (ev8
                    * (1.0 / y2 - 2.0 / (y2 * y2 + dsq * 0.5).sqrt()
                        + 1.0 / (y2 * y2 + dsq).sqrt())
                    - fg)
                    .powi(2);
            }
        }
        if f1 < f2 {
            a2 = y2;
        } else {
            a1 = y1;
        }
    }
    if f1 >= f2 {
        a2
    } else {
        a1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_core_elements() {
        let p = NddoParameters::mndo().unwrap();
        for z in [1u8, 6, 7, 8] {
            let e = p.element(z).unwrap();
            assert!(e.rho0 > 0.0, "rho0 must be positive for Z={z}");
        }
        // Hydrogen reference values from OpenMOPAC v23.2.5 MNDO.
        let h = p.element(1).unwrap();
        assert!((h.u_ss - (-11.906276)).abs() < 1e-9);
        assert!((h.beta_s - (-6.989064)).abs() < 1e-9);
        assert!((h.zeta_s - 1.331967).abs() < 1e-9);
        assert!((h.g_ss - 12.848).abs() < 1e-9);
        assert_eq!(h.n_orb, 1);
        assert_eq!(h.dd, 0.0);
        // Carbon.
        let c = p.element(6).unwrap();
        assert!((c.u_ss - (-52.279745)).abs() < 1e-9);
        assert!(c.dd > 0.0 && c.qq > 0.0 && c.rho1 > 0.0 && c.rho2 > 0.0);
        // H-H pair parameters.
        let nah = p.pair(11, 1).unwrap();
        assert!((nah.alpha - 2.509533).abs() < 1e-9);
        assert!((nah.x - 10.000001).abs() < 1e-9);
    }

    #[test]
    fn mndo_sulfur_uses_sp_orbitals() {
        let p = NddoParameters::mndo().unwrap();
        // MNDO sulfur uses the ordinary four-orbital s/p basis.
        let s = p.element(16).unwrap();
        assert!(!s.has_d(), "MNDO sulfur must use only s/p orbitals");
        assert_eq!(s.n_orb, 4);
        assert!(s.onecenter.is_none());
        // All additive terms and separations finite and positive.
        for i in 1..=9 {
            assert!(s.po[i].is_finite(), "po[{i}] not finite");
        }
        assert!(s.po[1] > 0.0 && s.po[2] > 0.0 && s.po[3] > 0.0);
        assert!(s.ddp[2] > 0.0 && s.ddp[3] > 0.0);
        // Core rho defaults to the ss additive term.
        assert!((s.po[9] - s.po[1]).abs() < 1e-12);
        // Isolated-atom energy is finite (d-shell eiscor folded in).
        assert!(s.e_isol.is_finite());
    }

    #[test]
    fn sp_element_multipoles_match_calpar_path() {
        let p = NddoParameters::mndo().unwrap();
        let o = p.element(8).unwrap();
        // MOPAC inid overwrite: po(1)=rho0, po(2)=rho1, po(3)=rho2, ddp(2)=dd,
        // ddp(3)=qq*sqrt(2).
        assert!((o.po[1] - o.rho0).abs() < 1e-12);
        assert!((o.po[2] - o.rho1).abs() < 1e-12);
        assert!((o.po[3] - o.rho2).abs() < 1e-12);
        assert!((o.ddp[2] - o.dd).abs() < 1e-12);
        assert!((o.ddp[3] - o.qq * 2.0_f64.sqrt()).abs() < 1e-12);
    }
}
