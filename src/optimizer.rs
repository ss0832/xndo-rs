// SPDX-License-Identifier: GPL-3.0-or-later

//! Rust-native L-BFGS geometry optimization.
//!
//! [`optimize`] uses the unified analytic-gradient dispatcher and therefore
//! supports every executable ground-state method. Positions are Bohr
//! internally; gradients are eV/Bohr.

use crate::error::{Result, XndoError};
use crate::math::Vec3;
use crate::method::Method;
use crate::scf::NddoOptions;
use crate::system::Molecule;
use crate::xndo::{run_gradient, run_method, CalculationResult};

#[derive(Clone, Debug)]
pub struct OptOptions {
    pub max_iter: usize,
    /// Convergence on the max gradient component (eV/Bohr).
    pub gtol: f64,
    /// L-BFGS history length.
    pub history: usize,
}

impl Default for OptOptions {
    fn default() -> Self {
        Self {
            max_iter: 200,
            gtol: 1.0e-3,
            history: 8,
        }
    }
}

#[derive(Clone, Debug)]
pub struct OptStep {
    pub energy_ev: f64,
    pub heat_of_formation_kcal: Option<f64>,
    pub max_gradient: f64,
    pub positions: Vec<Vec3>,
}

#[derive(Clone, Debug)]
pub struct OptResult {
    pub method: Method,
    pub molecule: Molecule,
    pub energy_ev: f64,
    pub heat_of_formation_kcal: Option<f64>,
    pub converged: bool,
    pub iterations: usize,
    pub trajectory: Vec<OptStep>,
}

/// Optimize the ground-state geometry for any executable native method.
///
/// The search direction uses L-BFGS with an Armijo backtracking line search.
/// Trial geometries require only a single-point energy; analytic gradients are
/// evaluated at accepted steps. Registered methods without a numerical engine
/// are rejected rather than substituted with a nearby Hamiltonian.
pub fn optimize(
    molecule: &Molecule,
    method: Method,
    scf_options: &NddoOptions,
    opt: &OptOptions,
) -> Result<OptResult> {
    if opt.max_iter == 0 {
        return Err(XndoError::InvalidInput(
            "optimization max_iter must be at least 1".to_owned(),
        ));
    }
    if !opt.gtol.is_finite() || opt.gtol <= 0.0 {
        return Err(XndoError::InvalidInput(
            "optimization gtol must be finite and greater than 0 eV/Bohr".to_owned(),
        ));
    }
    if opt.history == 0 {
        return Err(XndoError::InvalidInput(
            "optimization history must be at least 1".to_owned(),
        ));
    }
    if !method.supports_gradient() {
        return Err(XndoError::InvalidInput(method.execution_error()));
    }

    let nat = molecule.atoms.len();
    let ndof = 3 * nat;
    let mut mol = molecule.clone();
    let mut x = flatten(&mol);
    let gradient = run_gradient(&mol, method, scf_options)?;
    let mut g = flatten_grad(&gradient.gradient);
    let mut energy = gradient.energy_ev;
    let mut heat = gradient.heat_of_formation_kcal;
    let mut max_grad = max_gradient_component(&gradient.gradient);

    let mut s_hist: Vec<Vec<f64>> = Vec::new();
    let mut y_hist: Vec<Vec<f64>> = Vec::new();
    let mut rho_hist: Vec<f64> = Vec::new();
    let mut trajectory = vec![OptStep {
        energy_ev: energy,
        heat_of_formation_kcal: heat,
        max_gradient: max_grad,
        positions: unflatten(&x),
    }];
    let mut converged = max_grad < opt.gtol;
    let mut iterations = 0;

    for iter in 0..opt.max_iter {
        if converged {
            break;
        }
        iterations = iter + 1;

        let mut q = g.clone();
        let history_len = s_hist.len();
        let mut alpha = vec![0.0; history_len];
        for i in (0..history_len).rev() {
            let value = rho_hist[i] * dot(&s_hist[i], &q);
            alpha[i] = value;
            axpy(&mut q, -value, &y_hist[i]);
        }
        let scale = if history_len > 0 {
            let sy = dot(&s_hist[history_len - 1], &y_hist[history_len - 1]);
            let yy = dot(&y_hist[history_len - 1], &y_hist[history_len - 1]);
            if yy > 0.0 {
                sy / yy
            } else {
                1.0
            }
        } else {
            0.1 / max_grad.max(1.0e-6)
        };
        for value in &mut q {
            *value *= scale;
        }
        for i in 0..history_len {
            let beta = rho_hist[i] * dot(&y_hist[i], &q);
            axpy(&mut q, alpha[i] - beta, &s_hist[i]);
        }
        let mut direction: Vec<f64> = q.iter().map(|value| -value).collect();
        if dot(&direction, &g) > 0.0 {
            direction = g.iter().map(|value| -value).collect();
        }

        let directional_derivative = dot(&g, &direction);
        let mut step = 1.0;
        let mut x_new;
        let mut accepted = false;
        loop {
            x_new = x.clone();
            axpy(&mut x_new, step, &direction);
            set_positions(&mut mol, &x_new);
            if let Ok(result) = run_method(&mol, method, scf_options) {
                if calculation_energy(&result) <= energy + 1.0e-4 * step * directional_derivative {
                    accepted = true;
                    break;
                }
            }
            step *= 0.5;
            if step < 1.0e-8 {
                break;
            }
        }
        if !accepted {
            set_positions(&mut mol, &x);
            break;
        }

        let gradient_new = run_gradient(&mol, method, scf_options)?;
        let g_new = flatten_grad(&gradient_new.gradient);
        let s: Vec<f64> = (0..ndof).map(|i| x_new[i] - x[i]).collect();
        let y: Vec<f64> = (0..ndof).map(|i| g_new[i] - g[i]).collect();
        let sy = dot(&s, &y);
        if sy > 1.0e-10 {
            s_hist.push(s);
            y_hist.push(y);
            rho_hist.push(1.0 / sy);
            if s_hist.len() > opt.history {
                s_hist.remove(0);
                y_hist.remove(0);
                rho_hist.remove(0);
            }
        }

        x = x_new;
        g = g_new;
        energy = gradient_new.energy_ev;
        heat = gradient_new.heat_of_formation_kcal;
        max_grad = max_gradient_component(&gradient_new.gradient);
        converged = max_grad < opt.gtol;
        trajectory.push(OptStep {
            energy_ev: energy,
            heat_of_formation_kcal: heat,
            max_gradient: max_grad,
            positions: unflatten(&x),
        });
    }

    set_positions(&mut mol, &x);
    Ok(OptResult {
        method,
        molecule: mol,
        energy_ev: energy,
        heat_of_formation_kcal: heat,
        converged,
        iterations,
        trajectory,
    })
}

fn calculation_energy(result: &CalculationResult) -> f64 {
    match result {
        CalculationResult::CndoIndo(value) => value.total_ev,
        CalculationResult::Nddo(value) => value.total_ev,
        CalculationResult::Mindo3(value) => value.total_ev,
        CalculationResult::ZindoS(value) => value.total_ev,
    }
}

fn max_gradient_component(gradient: &[Vec3]) -> f64 {
    gradient
        .iter()
        .flat_map(|value| value.to_array())
        .fold(0.0_f64, |largest, value| largest.max(value.abs()))
}

fn flatten(mol: &Molecule) -> Vec<f64> {
    let mut v = Vec::with_capacity(3 * mol.atoms.len());
    for a in &mol.atoms {
        v.push(a.position.x);
        v.push(a.position.y);
        v.push(a.position.z);
    }
    v
}
fn flatten_grad(g: &[Vec3]) -> Vec<f64> {
    let mut v = Vec::with_capacity(3 * g.len());
    for gi in g {
        v.push(gi.x);
        v.push(gi.y);
        v.push(gi.z);
    }
    v
}
fn unflatten(x: &[f64]) -> Vec<Vec3> {
    x.chunks(3).map(|c| Vec3::new(c[0], c[1], c[2])).collect()
}
fn set_positions(mol: &mut Molecule, x: &[f64]) {
    for (i, a) in mol.atoms.iter_mut().enumerate() {
        a.position = Vec3::new(x[3 * i], x[3 * i + 1], x[3 * i + 2]);
    }
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn axpy(y: &mut [f64], a: f64, x: &[f64]) {
    for (yi, xi) in y.iter_mut().zip(x) {
        *yi += a * xi;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn optimizes_water() {
        // OpenMOPAC v23.2.5 MNDO optimization of the same distorted water.
        let xyz = "3\nwater\nO 0.0 0.0 0.0\nH 1.05 0.0 0.0\nH -0.30 1.02 0.0\n";
        let mol = Molecule::from_xyz_str(xyz, 0.0).unwrap();
        let res = optimize(
            &mol,
            Method::Mndo,
            &NddoOptions::default(),
            &OptOptions::default(),
        )
        .unwrap();
        eprintln!(
            "opt H2O: converged={} iters={} dHf={:.3} kcal/mol maxgrad={:.2e}",
            res.converged,
            res.iterations,
            res.heat_of_formation_kcal.unwrap(),
            res.trajectory.last().unwrap().max_gradient
        );
        assert!(res.converged);
        assert!((res.heat_of_formation_kcal.unwrap() - (-60.9470958081165)).abs() < 0.005);
        assert!(res.trajectory.last().unwrap().energy_ev < res.trajectory[0].energy_ev);
    }

    #[test]
    fn unified_optimizer_runs_every_native_method() {
        let xyz = "3\nwater\nO 0.0 0.0 0.0\nH 1.05 0.0 0.0\nH -0.30 1.02 0.0\n";
        let molecule = Molecule::from_xyz_str(xyz, 0.0).unwrap();
        let options = NddoOptions::default();
        let opt = OptOptions {
            max_iter: 100,
            gtol: 2.0e-3,
            ..OptOptions::default()
        };
        for method in [
            Method::Cndo2,
            Method::Indo,
            Method::Mndo,
            Method::MndoD,
            Method::Mindo3,
            Method::ZindoS,
        ] {
            let result = optimize(&molecule, method, &options, &opt).unwrap();
            let first = result.trajectory.first().unwrap();
            let last = result.trajectory.last().unwrap();
            assert!(result.converged, "{method}: {last:?}");
            assert!(last.energy_ev < first.energy_ev, "{method}");
            assert!(last.max_gradient < first.max_gradient, "{method}");
        }
    }
}
