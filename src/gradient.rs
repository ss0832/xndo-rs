// SPDX-License-Identifier: GPL-3.0-or-later

//! Nuclear gradients of the MNDO total energy.
//!
//! Because NDDO works in an orthonormal AO basis, the SCF energy is stationary with respect to
//! the density, so the nuclear gradient is the derivative of the energy expression at the
//! **fixed converged density** (there is no Pulay/overlap-constraint term). Three routines are
//! provided:
//!
//! * [`closed_form_gradient`]  the primary, **fully closed-form** gradient (forward-mode
//!   dual-number AD of every integral kernel; radial *and* angular overlap analytic for
//!   `n <= 3`). No SCF re-runs and no finite differences. This is what the optimizer uses.
//! * [`analytic_gradient`] - the same Hellmann-Feynman gradient with the electronic term taken
//!   by fixed-density central differences (the core-core term stays closed-form). Kept for the
//!   open-shell path and as a cross-check.
//! * [`numerical_gradient`] - a full-SCF central-difference gradient, kept as an independent
//!   correctness reference (each Cartesian component re-runs the SCF twice).

use crate::basis::Basis;
use crate::dual::Dual;
use crate::error::Result;
use crate::fock::build_fock;
use crate::hamiltonian::build_core;
use crate::linalg::Matrix;
use crate::math::Vec3;
use crate::params::NddoParameters;
use crate::repulsion::core_core_energy;
use crate::scf::{run_nddo_with_parameters, NddoOptions, NddoResult};
use crate::system::Molecule;

/// Electronic energy (eV) at a **fixed density** matrix (no SCF, no core-core term).
pub fn electronic_energy_at_fixed_density(
    molecule: &Molecule,
    params: &NddoParameters,
    density: &Matrix,
) -> Result<f64> {
    let basis = Basis::build(molecule, params)?;
    let core = build_core(molecule, &basis, params)?;
    let f = build_fock(molecule, &basis, params, &core, density)?;
    let electronic = 0.5 * (density.frobenius_dot(&core.h_core) + density.frobenius_dot(&f));
    Ok(electronic)
}

/// Total MNDO energy (eV) at a **fixed density** (electronic + core-core).
pub fn energy_at_fixed_density(
    molecule: &Molecule,
    params: &NddoParameters,
    density: &Matrix,
) -> Result<f64> {
    Ok(
        electronic_energy_at_fixed_density(molecule, params, density)?
            + core_core_energy(molecule, params)?,
    )
}

/// Hellmann-Feynman nuclear gradient. `step` is the displacement in Bohr (default 5e-4).
pub fn analytic_gradient(
    molecule: &Molecule,
    params: &NddoParameters,
    options: &NddoOptions,
    step: f64,
) -> Result<GradientResult> {
    use rayon::prelude::*;

    let scf = run_nddo_with_parameters(molecule, params, options)?;
    let energy_ev = scf.total_ev;
    let nat = molecule.atoms.len();
    let density = scf.density.clone();

    // Core-core repulsion: exact closed-form derivative.
    let mut gradient = crate::repulsion::core_core_gradient(molecule, params)?;

    // Electronic term: Hellmann-Feynman (fixed converged density) central difference of the
    // electronic energy only; the 3N components are independent, so run them on Rayon.
    let comps: Vec<(usize, usize)> = (0..nat).flat_map(|a| (0..3).map(move |k| (a, k))).collect();
    let electronic: Vec<(usize, usize, f64)> = comps
        .par_iter()
        .map(|&(a, k)| -> Result<(usize, usize, f64)> {
            let mut plus = molecule.clone();
            let mut minus = molecule.clone();
            displace(&mut plus.atoms[a].position, k, step);
            displace(&mut minus.atoms[a].position, k, -step);
            let ep = electronic_energy_at_fixed_density(&plus, params, &density)?;
            let em = electronic_energy_at_fixed_density(&minus, params, &density)?;
            Ok((a, k, (ep - em) / (2.0 * step)))
        })
        .collect::<Result<Vec<_>>>()?;
    for (a, k, g) in electronic {
        match k {
            0 => gradient[a].x += g,
            1 => gradient[a].y += g,
            _ => gradient[a].z += g,
        }
    }

    let forces: Vec<Vec3> = gradient.iter().map(|g| *g * -1.0).collect();
    let max_gradient = gradient
        .iter()
        .flat_map(|g| g.to_array())
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    Ok(GradientResult {
        scf,
        energy_ev,
        gradient,
        forces,
        max_gradient,
    })
}

#[derive(Clone, Debug)]
pub struct GradientResult {
    /// Converged SCF result at the input geometry.
    pub scf: NddoResult,
    /// Total energy (eV).
    pub energy_ev: f64,
    /// Gradient dE/dR in eV/Bohr (atomic-unit length).
    pub gradient: Vec<Vec3>,
    /// Forces = -gradient (eV/Bohr).
    pub forces: Vec<Vec3>,
    /// Largest gradient component magnitude (eV/Bohr).
    pub max_gradient: f64,
}

/// Finite-difference nuclear gradient. `step` is the displacement in Bohr (default 5e-4).
pub fn numerical_gradient(
    molecule: &Molecule,
    params: &NddoParameters,
    options: &NddoOptions,
    step: f64,
) -> Result<GradientResult> {
    let scf = run_nddo_with_parameters(molecule, params, options)?;
    let energy_ev = scf.total_ev;
    let nat = molecule.atoms.len();
    let mut gradient = vec![Vec3::zero(); nat];

    let energy_at = |m: &Molecule| -> Result<f64> {
        Ok(run_nddo_with_parameters(m, params, options)?.total_ev)
    };

    for (a, atom_gradient) in gradient.iter_mut().enumerate() {
        for k in 0..3 {
            let mut plus = molecule.clone();
            let mut minus = molecule.clone();
            displace(&mut plus.atoms[a].position, k, step);
            displace(&mut minus.atoms[a].position, k, -step);
            let ep = energy_at(&plus)?;
            let em = energy_at(&minus)?;
            let g = (ep - em) / (2.0 * step);
            set_component(atom_gradient, k, g);
        }
    }

    let forces: Vec<Vec3> = gradient.iter().map(|g| *g * -1.0).collect();
    let max_gradient = gradient
        .iter()
        .flat_map(|g| g.to_array())
        .fold(0.0_f64, |m, v| m.max(v.abs()));

    Ok(GradientResult {
        scf,
        energy_ev,
        gradient,
        forces,
        max_gradient,
    })
}

/// Fully closed-form (dual-number) Hellmann-Feynman gradient. The two-electron and
/// core-attraction integral derivatives, the overlap (radial *and* angular, for valence shells
/// `n <= 3`), and the core-core term are all exact forward-mode AD: no SCF re-runs and no
/// finite differences. Heavy-element and d-shell overlaps differentiate through
/// the general Slater quadrature with forward-mode AD.
/// Fully closed-form Hellmann-Feynman gradient, robust to axis-aligned d-containing geometries.
///
/// Thin wrapper over `closed_form_gradient_core`. The d-containing two-center rotation is singular
/// for a bond on the global z-axis, where the degenerate branch zeroes the rotation derivatives
/// so a d-pair sitting on `+z` in an **asymmetric** environment would get a wrong transverse
/// force (0.2 eV/Bohr in testing). Because the energy is rotationally invariant, we detect such a
/// pair, evaluate the gradient in a generic frame, and rotate it back (`g_i = R0T g'_i`; see
/// [`crate::frame`]). Geometries with no near-axis d-containing pair are bit-identical to the core.
pub fn closed_form_gradient(
    molecule: &Molecule,
    params: &NddoParameters,
    options: &NddoOptions,
) -> Result<GradientResult> {
    if !crate::frame::needs_generic_frame(molecule, params)? {
        return closed_form_gradient_core(molecule, params, options);
    }
    match crate::frame::pick_generic_frame(molecule, params)? {
        Some(r0) => {
            let rotated = crate::frame::rotate_molecule(molecule, &r0);
            let mut res = closed_form_gradient_core(&rotated, params, options)?;
            // Rotate the per-atom gradient/forces back to the input frame.
            res.gradient = crate::frame::back_rotate_gradient(&res.gradient, &r0);
            res.forces = res.gradient.iter().map(|g| *g * -1.0).collect();
            res.max_gradient = res
                .gradient
                .iter()
                .flat_map(|g| g.to_array())
                .fold(0.0_f64, |m, v| m.max(v.abs()));
            // Replace the rotated-frame SCF with one converged at the input geometry so callers
            // reading `res.scf` (density, MOs) get the original frame. Integral *values* are
            // exact at any orientation (only their derivatives needed the rotated frame), so this
            // SCF is correct; the extra solve only runs for the rare on-axis d-containing case.
            res.scf = run_nddo_with_parameters(molecule, params, options)?;
            res.energy_ev = res.scf.total_ev;
            Ok(res)
        }
        // No candidate frame cleared every d-pair off the z-axis (extremely unlikely):
        // the FD-of-energy gradient is correct at any orientation (energy values are exact).
        None => numerical_gradient(molecule, params, options, 5.0e-4),
    }
}

fn closed_form_gradient_core(
    molecule: &Molecule,
    params: &NddoParameters,
    options: &NddoOptions,
) -> Result<GradientResult> {
    let scf = run_nddo_with_parameters(molecule, params, options)?;
    if scf.unrestricted {
        // Open-shell: spin-resolved closed-form fixed-density (Hellmann-Feynman) gradient.
        let energy_ev = scf.total_ev;
        let gradient = fixed_density_gradient_uhf(molecule, params, &scf)?;
        let forces: Vec<Vec3> = gradient.iter().map(|g| *g * -1.0).collect();
        let max_gradient = gradient
            .iter()
            .flat_map(|g| g.to_array())
            .fold(0.0_f64, |m, v| m.max(v.abs()));
        return Ok(GradientResult {
            scf,
            energy_ev,
            gradient,
            forces,
            max_gradient,
        });
    }
    let energy_ev = scf.total_ev;
    let gradient = fixed_density_gradient(molecule, params, &scf.density)?;
    let forces: Vec<Vec3> = gradient.iter().map(|g| *g * -1.0).collect();
    let max_gradient = gradient
        .iter()
        .flat_map(|g| g.to_array())
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    Ok(GradientResult {
        scf,
        energy_ev,
        gradient,
        forces,
        max_gradient,
    })
}

/// Total closed-form gradient (core-core + electronic) at an **arbitrary fixed density** `p`
/// (no SCF solve). Finite-differencing this over the nuclei at fixed `p` gives the skeleton
/// (fixed-density) second derivative used by the analytic Hessian.
pub fn fixed_density_gradient(
    molecule: &Molecule,
    params: &NddoParameters,
    p: &Matrix,
) -> Result<Vec<Vec3>> {
    let basis = Basis::build(molecule, params)?;
    let mut gradient = crate::repulsion::core_core_gradient(molecule, params)?;
    let elec = electronic_gradient_fixed_density(molecule, params, &basis, p)?;
    for (g, e) in gradient.iter_mut().zip(&elec) {
        *g += *e;
    }
    Ok(gradient)
}

/// Resonance  for an orbital (0 = s, 1..4 = p, 4..9 = d).
#[inline]
pub(crate) fn resonance_beta(elem: &crate::params::NddoElement, orb: u8) -> f64 {
    match orb {
        0 => elem.beta_s,
        1..=3 => elem.beta_p,
        _ => elem.beta_d,
    }
}

/// Per-pair dual-number two-electron integrals + overlap (99), dispatching to
/// the spd path when either atom has d orbitals. Returns the ordered atom
/// indices `(a, b)` (`te.e1b`/overlap rows belong to `a`, `e2a`/cols to `b`).
pub(crate) type PairDual = (
    usize,
    usize,
    crate::integrals::PairTwoElecG<Dual>,
    [[Dual; 9]; 9],
);
pub(crate) fn pair_dual(
    molecule: &Molecule,
    params: &NddoParameters,
    u: usize,
    v: usize,
) -> Result<PairDual> {
    use crate::integrals::pair_two_electron_dual;
    use crate::integrals_d::pair_two_electron_spd_dual;
    use crate::overlap::{diatom_overlap_dual, diatom_overlap_spd_dual, embed4_dual};
    let eu = params.element(molecule.atoms[u].z)?;
    let ev = params.element(molecule.atoms[v].z)?;
    if eu.n_orb == 0 || ev.n_orb == 0 {
        let (a, b) = if eu.n_orb > 0 { (u, v) } else { (v, u) };
        let ea = params.element(molecule.atoms[a].z)?;
        let eb = params.element(molecule.atoms[b].z)?;
        let displacement = molecule.atoms[b].position - molecule.atoms[a].position;
        let te = crate::integrals::pair_with_point_core_g::<Dual>(
            ea,
            eb,
            [
                Dual::var(displacement.x, 0),
                Dual::var(displacement.y, 1),
                Dual::var(displacement.z, 2),
            ],
        );
        Ok((a, b, te, [[Dual::constant(0.0); 9]; 9]))
    } else if eu.has_d() || ev.has_d() {
        let d = molecule.atoms[v].position - molecule.atoms[u].position;
        let te = pair_two_electron_spd_dual(eu, ev, d);
        let s = diatom_overlap_spd_dual(eu, ev, d);
        Ok((u, v, te, s))
    } else {
        let (a, b) = if eu.has_p() || !ev.has_p() {
            (u, v)
        } else {
            (v, u)
        };
        let ea = params.element(molecule.atoms[a].z)?;
        let eb = params.element(molecule.atoms[b].z)?;
        let (pa, pb) = (molecule.atoms[a].position, molecule.atoms[b].position);
        let te = pair_two_electron_dual(ea, eb, pb - pa);
        let s = embed4_dual(diatom_overlap_dual(ea, pa, eb, pb)?);
        Ok((a, b, te, s))
    }
}

pub fn electronic_gradient_fixed_density(
    molecule: &Molecule,
    params: &NddoParameters,
    basis: &Basis,
    p: &Matrix,
) -> Result<Vec<Vec3>> {
    let nat = molecule.atoms.len();
    let mut gradient = vec![Vec3::zero(); nat];
    for u in 0..nat {
        for v in (u + 1)..nat {
            let (a, b, te, s) = pair_dual(molecule, params, u, v)?;
            let ea = params.element(molecule.atoms[a].z)?;
            let eb = params.element(molecule.atoms[b].z)?;
            let (oa, ob) = (basis.atom_offset[a], basis.atom_offset[b]);
            let (na, nb) = (basis.atom_norb[a], basis.atom_norb[b]);
            let mut f = [0.0_f64; 3];
            for i in 0..na {
                let bi = resonance_beta(ea, basis.aos[oa + i].orb);
                for j in 0..nb {
                    let bj = resonance_beta(eb, basis.aos[ob + j].orb);
                    let coef = p[(oa + i, ob + j)] * (bi + bj);
                    for (ax, fx) in f.iter_mut().enumerate() {
                        *fx += coef * s[i][j].d[ax];
                    }
                }
            }
            for i in 0..na {
                for j in 0..na {
                    let coef = p[(oa + i, oa + j)];
                    for (ax, fx) in f.iter_mut().enumerate() {
                        *fx += coef * te.e1b[i][j].d[ax];
                    }
                }
            }
            for k in 0..nb {
                for l in 0..nb {
                    let coef = p[(ob + k, ob + l)];
                    for (ax, fx) in f.iter_mut().enumerate() {
                        *fx += coef * te.e2a[k][l].d[ax];
                    }
                }
            }
            for mu in 0..na {
                for nu in 0..na {
                    for la in 0..nb {
                        for si in 0..nb {
                            let dw = te.two_e(mu, nu, la, si).d;
                            let coul = p[(oa + mu, oa + nu)] * p[(ob + la, ob + si)];
                            let exch = -0.5 * p[(oa + mu, ob + la)] * p[(oa + nu, ob + si)];
                            let coef = coul + exch;
                            for (ax, fx) in f.iter_mut().enumerate() {
                                *fx += coef * dw[ax];
                            }
                        }
                    }
                }
            }
            gradient[b] += Vec3::new(f[0], f[1], f[2]);
            gradient[a] -= Vec3::new(f[0], f[1], f[2]);
        }
    }
    Ok(gradient)
}

/// Total closed-form UHF gradient (core-core + spin-resolved electronic) at the converged
/// open-shell density. `P_alpha = (P_tot + S)/2` and `P_beta = (P_tot - S)/2` are
/// reconstructed from total density and spin density `S = P_alpha - P_beta`.
/// This is Hellmann-Feynman differentiation in the orthonormal NDDO basis.
pub fn fixed_density_gradient_uhf(
    molecule: &Molecule,
    params: &NddoParameters,
    scf: &NddoResult,
) -> Result<Vec<Vec3>> {
    let basis = Basis::build(molecule, params)?;
    let pt = &scf.density;
    let spin = scf.spin_density.as_ref().ok_or_else(|| {
        crate::error::XndoError::InvalidInput("UHF gradient requires a spin density".into())
    })?;
    let mut pa = pt.clone();
    let mut pb = pt.clone();
    {
        let n = pt.as_slice().len();
        let (pas, pbs) = (pa.as_mut_slice(), pb.as_mut_slice());
        let (pts, ss) = (pt.as_slice(), spin.as_slice());
        for i in 0..n {
            pas[i] = 0.5 * (pts[i] + ss[i]);
            pbs[i] = 0.5 * (pts[i] - ss[i]);
        }
    }
    let mut gradient = crate::repulsion::core_core_gradient(molecule, params)?;
    let elec = electronic_gradient_fixed_density_spin(molecule, params, &basis, pt, &pa, &pb)?;
    for (g, e) in gradient.iter_mut().zip(&elec) {
        *g += *e;
    }
    Ok(gradient)
}

/// Spin-resolved electronic part of the closed-form gradient at fixed densities: resonance,
/// electroncore attraction, and Coulomb use the **total** density `P_tot`; exchange uses the
/// **same-spin** densities `P`, `P` (`[P_ P_ + P_ P_](|)`). Reduces to the
/// RHF form when `P = P = P_tot/2`.
pub fn electronic_gradient_fixed_density_spin(
    molecule: &Molecule,
    params: &NddoParameters,
    basis: &Basis,
    pt: &Matrix,
    pa: &Matrix,
    pb: &Matrix,
) -> Result<Vec<Vec3>> {
    let nat = molecule.atoms.len();
    let mut gradient = vec![Vec3::zero(); nat];
    for u in 0..nat {
        for v in (u + 1)..nat {
            let (a, b, te, s) = pair_dual(molecule, params, u, v)?;
            let ea = params.element(molecule.atoms[a].z)?;
            let eb = params.element(molecule.atoms[b].z)?;
            let (oa, ob) = (basis.atom_offset[a], basis.atom_offset[b]);
            let (na, nb) = (basis.atom_norb[a], basis.atom_norb[b]);
            let mut f = [0.0_f64; 3];
            // Resonance S (total density).
            for i in 0..na {
                let bi = resonance_beta(ea, basis.aos[oa + i].orb);
                for j in 0..nb {
                    let bj = resonance_beta(eb, basis.aos[ob + j].orb);
                    let coef = pt[(oa + i, ob + j)] * (bi + bj);
                    for (ax, fx) in f.iter_mut().enumerate() {
                        *fx += coef * s[i][j].d[ax];
                    }
                }
            }
            // Electroncore attraction (total density).
            for i in 0..na {
                for j in 0..na {
                    let coef = pt[(oa + i, oa + j)];
                    for (ax, fx) in f.iter_mut().enumerate() {
                        *fx += coef * te.e1b[i][j].d[ax];
                    }
                }
            }
            for k in 0..nb {
                for l in 0..nb {
                    let coef = pt[(ob + k, ob + l)];
                    for (ax, fx) in f.iter_mut().enumerate() {
                        *fx += coef * te.e2a[k][l].d[ax];
                    }
                }
            }
            // Two-electron: Coulomb from P_tot, exchange from same-spin P/P.
            for mu in 0..na {
                for nu in 0..na {
                    for la in 0..nb {
                        for si in 0..nb {
                            let dw = te.two_e(mu, nu, la, si).d;
                            let coul = pt[(oa + mu, oa + nu)] * pt[(ob + la, ob + si)];
                            let exch = -(pa[(oa + mu, ob + la)] * pa[(oa + nu, ob + si)]
                                + pb[(oa + mu, ob + la)] * pb[(oa + nu, ob + si)]);
                            let coef = coul + exch;
                            for (ax, fx) in f.iter_mut().enumerate() {
                                *fx += coef * dw[ax];
                            }
                        }
                    }
                }
            }
            gradient[b] += Vec3::new(f[0], f[1], f[2]);
            gradient[a] -= Vec3::new(f[0], f[1], f[2]);
        }
    }
    Ok(gradient)
}

#[inline]
fn displace(p: &mut Vec3, k: usize, d: f64) {
    match k {
        0 => p.x += d,
        1 => p.y += d,
        _ => p.z += d,
    }
}

#[inline]
fn set_component(v: &mut Vec3, k: usize, val: f64) {
    match k {
        0 => v.x = val,
        1 => v.y = val,
        _ => v.z = val,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn analytic_matches_full_scf_gradient() {
        // The Hellmann-Feynman (fixed-density) gradient must match the full-SCF finite
        // difference on a molecule displaced away from equilibrium (nonzero forces).
        let mol = Molecule::from_xyz_str(
            "3\nwater\nO 0.0 0.0 0.0\nH 1.02 0.0 0.0\nH -0.28 0.96 0.0\n",
            0.0,
        )
        .unwrap();
        let params = NddoParameters::mndo().unwrap();
        let opts = NddoOptions::default();
        let a = analytic_gradient(&mol, &params, &opts, 1.0e-4).unwrap();
        let n = numerical_gradient(&mol, &params, &opts, 1.0e-4).unwrap();
        let mut max_delta = 0.0_f64;
        for (ga, gn) in a.gradient.iter().zip(&n.gradient) {
            for k in 0..3 {
                max_delta = max_delta.max((ga.get(k) - gn.get(k)).abs());
            }
        }
        eprintln!("analytic-vs-numerical gradient max delta = {max_delta:.3e} eV/Bohr");
        assert!(max_delta < 1.0e-4, "gradient mismatch {max_delta:.3e}");
        // Forces must be nonzero for this distorted geometry.
        assert!(a.max_gradient > 1.0e-2);
    }

    #[test]
    fn water_gradient_matches_mopac_reference() {
        // OpenMOPAC v23.2.5: `MNDO PRECISE 1SCF GRADIENTS` at this exact XYZ
        // geometry. MOPAC reports kcal mol^-1 Angstrom^-1; the kernel uses
        // eV Bohr^-1. Keeping the oracle values in the test makes the
        // gradient criterion reproducible without requiring MOPAC at test time.
        let mol = Molecule::from_xyz_str(
            "3\nwater\nO 0.0 0.0 0.0\nH 1.05 0.0 0.0\nH -0.30 1.02 0.0\n",
            0.0,
        )
        .unwrap();
        let oracle_kcal_per_ang = [
            [-71.683_451_236_528_7, -108.175_252_702_434_2, 0.0],
            [112.372_811_821_249_6, -6.704_146_187_807_9, 0.0],
            [-40.689_360_520_280_1, 114.879_398_957_316_5, 0.0],
        ];
        let scale = crate::constants::KCAL_TO_EV * crate::constants::BOHR_TO_ANGSTROM;
        let params = NddoParameters::mndo().unwrap();
        let got = closed_form_gradient(&mol, &params, &NddoOptions::default()).unwrap();
        let mut max_delta = 0.0_f64;
        for (g, reference) in got.gradient.iter().zip(oracle_kcal_per_ang) {
            for (axis, value) in reference.into_iter().enumerate() {
                max_delta = max_delta.max((g.get(axis) - value * scale).abs());
            }
        }
        // 0.02 kcal mol^-1 Angstrom^-1 is the M3 target from the plan.
        let target = 0.02 * scale;
        assert!(
            max_delta < target,
            "MOPAC gradient mismatch {max_delta:.3e}"
        );
    }

    #[test]
    fn closed_form_matches_numerical_gradient() {
        // The fully closed-form (dual-number) gradient must match the full-SCF finite
        // difference on a molecule with s and p atoms displaced from equilibrium.
        let mol = Molecule::from_xyz_str(
            "4\nformaldehyde\nC 0.0 0.0 0.0\nO 0.03 0.0 1.25\nH 0.95 0.02 -0.55\nH -0.94 -0.03 -0.52\n",
            0.0,
        )
        .unwrap();
        let params = NddoParameters::mndo().unwrap();
        let opts = NddoOptions::default();
        let cf = closed_form_gradient(&mol, &params, &opts).unwrap();
        let n = numerical_gradient(&mol, &params, &opts, 1.0e-4).unwrap();
        let mut max_delta = 0.0_f64;
        for (gc, gn) in cf.gradient.iter().zip(&n.gradient) {
            for k in 0..3 {
                max_delta = max_delta.max((gc.get(k) - gn.get(k)).abs());
            }
        }
        eprintln!("closed-form-vs-numerical gradient max delta = {max_delta:.3e} eV/Bohr");
        assert!(
            max_delta < 5.0e-5,
            "closed-form gradient mismatch {max_delta:.3e}"
        );
    }

    #[test]
    fn analytic_d_orbital_gradient_matches_energy_differences() {
        // Br has a MNDO valence d shell. The analytic spd-AD gradient must
        // agree with the independent full-SCF central difference.
        let mol = Molecule::from_xyz_str(
            "5\nCH3Br\nC 0.0 0.0 0.0\nBr 0.0 0.0 -2.10\nH 1.03 0.0 0.40\nH -0.515 0.892 0.40\nH -0.515 -0.892 0.40\n",
            0.0,
        )
        .unwrap();
        let params = NddoParameters::mndo().unwrap();
        let opts = NddoOptions::default();
        let cf = closed_form_gradient(&mol, &params, &opts).unwrap();
        let num = numerical_gradient(&mol, &params, &opts, 5.0e-4).unwrap();
        let mut max_delta = 0.0f64;
        for (a, b) in cf.gradient.iter().zip(&num.gradient) {
            for k in 0..3 {
                max_delta = max_delta.max((a.get(k) - b.get(k)).abs());
            }
        }
        assert!(
            max_delta < 1e-6,
            "analytic d gradient mismatch {max_delta:.3e}"
        );
    }

    #[test]
    fn closed_form_gradient_uhf_radical() {
        // Methyl radical (doublet, UHF), distorted from planar: the spin-resolved closed-form
        // gradient must match the full-SCF finite difference (no fixed-density FD fallback).
        let mol = Molecule::from_xyz_str(
            "4\nmethyl\nC 0.0 0.0 0.05\nH 1.12 0.0 0.0\nH -0.55 0.95 0.0\nH -0.55 -0.95 0.0\n",
            0.0,
        )
        .unwrap();
        let params = NddoParameters::mndo().unwrap();
        let opts = NddoOptions {
            multiplicity: 2,
            ..NddoOptions::default()
        };
        let cf = closed_form_gradient(&mol, &params, &opts).unwrap();
        let n = numerical_gradient(&mol, &params, &opts, 1.0e-4).unwrap();
        assert!(cf.scf.unrestricted);
        let mut max_delta = 0.0_f64;
        for (gc, gn) in cf.gradient.iter().zip(&n.gradient) {
            for k in 0..3 {
                max_delta = max_delta.max((gc.get(k) - gn.get(k)).abs());
            }
        }
        eprintln!("UHF closed-form-vs-numerical gradient max delta = {max_delta:.3e}");
        assert!(max_delta < 5.0e-5, "UHF gradient mismatch {max_delta:.3e}");
        assert!(cf.max_gradient > 1.0e-2);
    }
}
