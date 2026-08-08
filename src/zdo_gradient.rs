// SPDX-License-Identifier: GPL-3.0-or-later

//! Shared contractions for analytic gradients of orthogonal ZDO Hamiltonians.

use crate::dual::Scalar;
use crate::error::{Result, XndoError};
use crate::linalg::{solve_linear, solve_linear_matrix, Matrix};

pub(crate) fn atom_population(density: &Matrix, offset: usize, norb: usize) -> f64 {
    (0..norb).map(|i| density[(offset + i, offset + i)]).sum()
}

pub(crate) fn spin_densities(total: &Matrix, spin: &Matrix) -> (Matrix, Matrix) {
    debug_assert_eq!(total.rows, spin.rows);
    debug_assert_eq!(total.cols, spin.cols);
    let mut alpha = Matrix::zeros(total.rows, total.cols);
    let mut beta = Matrix::zeros(total.rows, total.cols);
    for i in 0..total.as_slice().len() {
        alpha.as_mut_slice()[i] = 0.5 * (total.as_slice()[i] + spin.as_slice()[i]);
        beta.as_mut_slice()[i] = 0.5 * (total.as_slice()[i] - spin.as_slice()[i]);
    }
    (alpha, beta)
}

/// Assemble the geometry-dependent energy of one atom pair. All scalar
/// coefficients are evaluated at the converged density; only `gamma`,
/// `resonance`, and `core` carry nuclear derivatives.
#[allow(clippy::too_many_arguments)]
pub(crate) fn pair_energy<S: Scalar>(
    gamma: S,
    core: S,
    core_a: f64,
    core_b: f64,
    population_a: f64,
    population_b: f64,
    exchange: f64,
    resonance: S,
) -> S {
    let gamma_weight =
        -core_b * population_a - core_a * population_b + population_a * population_b - exchange;
    gamma * gamma_weight + resonance * 2.0 + core
}

pub(crate) fn exchange_weight_rhf(
    density: &Matrix,
    offset_a: usize,
    norb_a: usize,
    offset_b: usize,
    norb_b: usize,
) -> f64 {
    let mut value = 0.0;
    for i in 0..norb_a {
        for j in 0..norb_b {
            let p = density[(offset_a + i, offset_b + j)];
            value += 0.5 * p * p;
        }
    }
    value
}

pub(crate) fn exchange_weight_uhf(
    alpha: &Matrix,
    beta: &Matrix,
    offset_a: usize,
    norb_a: usize,
    offset_b: usize,
    norb_b: usize,
) -> f64 {
    let mut value = 0.0;
    for i in 0..norb_a {
        for j in 0..norb_b {
            let pa = alpha[(offset_a + i, offset_b + j)];
            let pb = beta[(offset_a + i, offset_b + j)];
            value += pa * pa + pb * pb;
        }
    }
    value
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

fn transform_to_mo(c: &Matrix, ao: &Matrix) -> Matrix {
    c.transpose_matmul_seq(ao).matmul_seq(c)
}

fn rotation_density(c: &Matrix, n_occ: usize, amplitudes: &[f64], occupation: f64) -> Matrix {
    let nmo = c.cols;
    let nvirt = nmo - n_occ;
    debug_assert_eq!(amplitudes.len(), n_occ * nvirt);
    let mut density = Matrix::zeros(c.rows, c.rows);
    for i in 0..n_occ {
        for av in 0..nvirt {
            let a = n_occ + av;
            let u = amplitudes[i * nvirt + av] * occupation;
            if u == 0.0 {
                continue;
            }
            for mu in 0..c.rows {
                let left = c[(mu, a)] * u;
                let occupied_mu = c[(mu, i)] * u;
                for nu in 0..c.rows {
                    density[(mu, nu)] += left * c[(nu, i)] + occupied_mu * c[(nu, a)];
                }
            }
        }
    }
    density
}

fn coefficient_derivative(c: &Matrix, eps: &[f64], fock_mo: &Matrix) -> Result<Matrix> {
    let nmo = c.cols;
    let mut derivative = Matrix::zeros(c.rows, nmo);
    for p in 0..nmo {
        for q in 0..nmo {
            if p == q {
                continue;
            }
            let gap = eps[p] - eps[q];
            if gap.abs() < 1.0e-10 {
                if fock_mo[(q, p)].abs() > 1.0e-7 {
                    return Err(XndoError::LinearAlgebra(format!(
                        "degenerate orbital response ({p},{q}) with nonzero coupling {:.3e}",
                        fock_mo[(q, p)]
                    )));
                }
                continue;
            }
            let u = fock_mo[(q, p)] / gap;
            for mu in 0..c.rows {
                derivative[(mu, p)] += c[(mu, q)] * u;
            }
        }
    }
    Ok(derivative)
}

#[derive(Clone, Debug)]
pub(crate) struct RhfResponse {
    pub density: Matrix,
    pub fock_mo: Matrix,
    pub mo_coeff: Matrix,
    pub mo_energies: Vec<f64>,
}

/// Solve the static RHF CPHF equations for several Cartesian perturbations.
/// The orbital-response matrix is factored only once.
pub(crate) fn solve_rhf_responses<F>(
    c: &Matrix,
    eps: &[f64],
    n_occ: usize,
    skeleton_focks: &[Matrix],
    response_fock: F,
) -> Result<Vec<RhfResponse>>
where
    F: Fn(&Matrix) -> Result<Matrix>,
{
    let nvirt = c.cols.saturating_sub(n_occ);
    let nov = n_occ * nvirt;
    if nov == 0 {
        return skeleton_focks
            .iter()
            .map(|skeleton_fock| {
                let fock_mo = transform_to_mo(c, skeleton_fock);
                Ok(RhfResponse {
                    density: Matrix::zeros(c.rows, c.rows),
                    mo_coeff: coefficient_derivative(c, eps, &fock_mo)?,
                    mo_energies: (0..c.cols).map(|p| fock_mo[(p, p)]).collect(),
                    fock_mo,
                })
            })
            .collect();
    }
    let skeleton_mos: Vec<Matrix> = skeleton_focks
        .iter()
        .map(|skeleton| transform_to_mo(c, skeleton))
        .collect();
    let mut rhs = Matrix::zeros(nov, skeleton_focks.len());
    let mut system = Matrix::identity(nov);
    for (coordinate, skeleton_mo) in skeleton_mos.iter().enumerate() {
        for i in 0..n_occ {
            for av in 0..nvirt {
                let row = i * nvirt + av;
                let a = n_occ + av;
                let gap = eps[i] - eps[a];
                if gap.abs() < 1.0e-10 {
                    return Err(XndoError::LinearAlgebra(format!(
                        "vanishing RHF occupied-virtual gap ({i},{a})"
                    )));
                }
                rhs[(row, coordinate)] = skeleton_mo[(a, i)] / gap;
            }
        }
    }
    for col in 0..nov {
        let mut unit = vec![0.0; nov];
        unit[col] = 1.0;
        let trial_density = rotation_density(c, n_occ, &unit, 2.0);
        let induced_mo = transform_to_mo(c, &response_fock(&trial_density)?);
        for i in 0..n_occ {
            for av in 0..nvirt {
                let row = i * nvirt + av;
                let a = n_occ + av;
                system[(row, col)] -= induced_mo[(a, i)] / (eps[i] - eps[a]);
            }
        }
    }
    let amplitudes = solve_linear_matrix(&system, &rhs)?;
    skeleton_focks
        .iter()
        .enumerate()
        .map(|(coordinate, skeleton_fock)| {
            let column: Vec<f64> = (0..nov)
                .map(|index| amplitudes[(index, coordinate)])
                .collect();
            let density = rotation_density(c, n_occ, &column, 2.0);
            let fock = add_matrix(skeleton_fock, &response_fock(&density)?);
            let fock_mo = transform_to_mo(c, &fock);
            Ok(RhfResponse {
                density,
                mo_coeff: coefficient_derivative(c, eps, &fock_mo)?,
                mo_energies: (0..c.cols).map(|p| fock_mo[(p, p)]).collect(),
                fock_mo,
            })
        })
        .collect()
}

#[derive(Clone, Debug)]
pub(crate) struct UhfResponse {
    pub density_alpha: Matrix,
    pub density_beta: Matrix,
    pub fock_alpha_mo: Matrix,
    pub fock_beta_mo: Matrix,
    pub mo_coeff_alpha: Matrix,
    pub mo_coeff_beta: Matrix,
    pub mo_energies_alpha: Vec<f64>,
    pub mo_energies_beta: Vec<f64>,
}

#[derive(Clone, Debug)]
pub(crate) struct UhfSecondResponse {
    pub mo_coeff_alpha: Matrix,
    pub mo_coeff_beta: Matrix,
    pub mo_energies_alpha: Vec<f64>,
    pub mo_energies_beta: Vec<f64>,
}

#[derive(Clone, Debug)]
pub(crate) struct RhfSecondResponse {
    pub mo_coeff: Matrix,
    pub mo_energies: Vec<f64>,
}

fn paired_response_diis(
    history_alpha: &[Matrix],
    history_beta: &[Matrix],
    errors_alpha: &[Matrix],
    errors_beta: &[Matrix],
) -> Option<(Matrix, Matrix)> {
    let count = history_alpha.len();
    if count < 2
        || history_beta.len() != count
        || errors_alpha.len() != count
        || errors_beta.len() != count
    {
        return None;
    }
    let scale = (0..count)
        .map(|i| {
            (errors_alpha[i].frobenius_dot(&errors_alpha[i])
                + errors_beta[i].frobenius_dot(&errors_beta[i]))
            .sqrt()
            .max(1.0e-300)
        })
        .collect::<Vec<_>>();
    let dimension = count + 1;
    let mut system = Matrix::zeros(dimension, dimension);
    for i in 0..count {
        for j in 0..count {
            let mut value = (errors_alpha[i].frobenius_dot(&errors_alpha[j])
                + errors_beta[i].frobenius_dot(&errors_beta[j]))
                / (scale[i] * scale[j]);
            if i == j {
                value += 1.0e-10;
            }
            system[(i, j)] = value;
        }
        system[(i, count)] = -1.0;
        system[(count, i)] = -1.0;
    }
    let mut rhs = vec![0.0; dimension];
    rhs[count] = -1.0;
    let raw = solve_linear(&system, &rhs).ok()?;
    let mut coefficients = (0..count).map(|i| raw[i] / scale[i]).collect::<Vec<_>>();
    let sum = coefficients.iter().sum::<f64>();
    if !sum.is_finite() || sum.abs() < 1.0e-300 {
        return None;
    }
    for coefficient in &mut coefficients {
        *coefficient /= sum;
    }
    let mut alpha = Matrix::zeros(history_alpha[0].rows, history_alpha[0].cols);
    let mut beta = Matrix::zeros(history_beta[0].rows, history_beta[0].cols);
    for (index, coefficient) in coefficients.into_iter().enumerate() {
        for element in 0..alpha.as_slice().len() {
            alpha.as_mut_slice()[element] += coefficient * history_alpha[index].as_slice()[element];
            beta.as_mut_slice()[element] += coefficient * history_beta[index].as_slice()[element];
        }
    }
    Some((alpha, beta))
}

fn mixed_density(
    c: &Matrix,
    cx: &Matrix,
    cy: &Matrix,
    cxy: &Matrix,
    n_occ: usize,
    occupation: f64,
) -> Matrix {
    let mut density = Matrix::zeros(c.rows, c.rows);
    for i in 0..n_occ {
        for mu in 0..c.rows {
            for nu in 0..c.rows {
                density[(mu, nu)] += occupation
                    * (cxy[(mu, i)] * c[(nu, i)]
                        + cx[(mu, i)] * cy[(nu, i)]
                        + cy[(mu, i)] * cx[(nu, i)]
                        + c[(mu, i)] * cxy[(nu, i)]);
            }
        }
    }
    density
}

/// Solve the mixed second-order RHF response by analytic fixed-point
/// iteration. The map is the differentiated CPHF fixed point; no displaced
/// geometry or finite difference is involved.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_rhf_second_response<F>(
    c: &Matrix,
    eps: &[f64],
    n_occ: usize,
    first_x: &RhfResponse,
    first_y: &RhfResponse,
    known_fock_xy: &Matrix,
    response_fock: F,
) -> Result<RhfSecondResponse>
where
    F: Fn(&Matrix) -> Result<Matrix>,
{
    let ux = c.transpose_matmul_seq(&first_x.mo_coeff);
    let uy = c.transpose_matmul_seq(&first_y.mo_coeff);
    let mut density = mixed_density(
        c,
        &first_x.mo_coeff,
        &first_y.mo_coeff,
        &Matrix::zeros(c.rows, c.cols),
        n_occ,
        2.0,
    );
    for iteration in 0..500 {
        let final_fock = add_matrix(known_fock_xy, &response_fock(&density)?);
        let final_fock_mo = transform_to_mo(c, &final_fock);
        let mut final_eps = vec![0.0; c.cols];
        let mut v = Matrix::zeros(c.cols, c.cols);
        for p in 0..c.cols {
            for q in 0..c.cols {
                if p == q {
                    continue;
                }
                let gap = eps[p] - eps[q];
                if gap.abs() < 1.0e-10 {
                    continue;
                }
                let mut numerator = final_fock_mo[(q, p)]
                    - first_x.mo_energies[p] * uy[(q, p)]
                    - first_y.mo_energies[p] * ux[(q, p)];
                for r in 0..c.cols {
                    numerator +=
                        first_x.fock_mo[(q, r)] * uy[(r, p)] + first_y.fock_mo[(q, r)] * ux[(r, p)];
                }
                v[(q, p)] = numerator / gap;
            }
            let mut normalization = 0.0;
            for q in 0..c.cols {
                normalization += ux[(q, p)] * uy[(q, p)];
            }
            v[(p, p)] = -normalization;
            let mut energy = final_fock_mo[(p, p)];
            for q in 0..c.cols {
                energy +=
                    first_x.fock_mo[(p, q)] * uy[(q, p)] + first_y.fock_mo[(p, q)] * ux[(q, p)];
            }
            final_eps[p] = energy;
        }
        let final_cxy = c.matmul_seq(&v);
        let next = mixed_density(
            c,
            &first_x.mo_coeff,
            &first_y.mo_coeff,
            &final_cxy,
            n_occ,
            2.0,
        );
        let error = next.rms_difference(&density);
        if error < 1.0e-11 {
            return Ok(RhfSecondResponse {
                mo_coeff: final_cxy,
                mo_energies: final_eps,
            });
        }
        for index in 0..density.as_slice().len() {
            density.as_mut_slice()[index] =
                0.45 * next.as_slice()[index] + 0.55 * density.as_slice()[index];
        }
        if iteration == 499 {
            return Err(XndoError::LinearAlgebra(format!(
                "second-order RHF response did not converge (density RMS {error:.3e})"
            )));
        }
    }
    unreachable!()
}

/// Solve the coupled alpha/beta static UHF response equations.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_uhf_responses<F>(
    ca: &Matrix,
    cb: &Matrix,
    eps_a: &[f64],
    eps_b: &[f64],
    n_alpha: usize,
    n_beta: usize,
    skeleton_alpha: &[Matrix],
    skeleton_beta: &[Matrix],
    response_fock: F,
) -> Result<Vec<UhfResponse>>
where
    F: Fn(&Matrix, &Matrix) -> Result<(Matrix, Matrix)>,
{
    if skeleton_alpha.len() != skeleton_beta.len() {
        return Err(XndoError::LinearAlgebra(
            "UHF response skeleton count mismatch".to_string(),
        ));
    }
    let nva = ca.cols - n_alpha;
    let nvb = cb.cols - n_beta;
    let nova = n_alpha * nva;
    let novb = n_beta * nvb;
    let nov = nova + novb;
    for i in 0..n_alpha {
        for a in n_alpha..ca.cols {
            if (eps_a[i] - eps_a[a]).abs() < 1.0e-9 {
                return Err(XndoError::LinearAlgebra(format!(
                    "degenerate alpha occupied-virtual UHF response ({i},{a}); a state-specific Hessian is not unique"
                )));
            }
        }
    }
    for i in 0..n_beta {
        for a in n_beta..cb.cols {
            if (eps_b[i] - eps_b[a]).abs() < 1.0e-9 {
                return Err(XndoError::LinearAlgebra(format!(
                    "degenerate beta occupied-virtual UHF response ({i},{a}); a state-specific Hessian is not unique"
                )));
            }
        }
    }
    let skel_a_mo: Vec<Matrix> = skeleton_alpha
        .iter()
        .map(|skeleton| transform_to_mo(ca, skeleton))
        .collect();
    let skel_b_mo: Vec<Matrix> = skeleton_beta
        .iter()
        .map(|skeleton| transform_to_mo(cb, skeleton))
        .collect();
    let mut rhs = Matrix::zeros(nov, skeleton_alpha.len());
    let mut system = Matrix::identity(nov);
    for coordinate in 0..skeleton_alpha.len() {
        for i in 0..n_alpha {
            for av in 0..nva {
                let a = n_alpha + av;
                rhs[(i * nva + av, coordinate)] =
                    skel_a_mo[coordinate][(a, i)] / (eps_a[i] - eps_a[a]);
            }
        }
        for i in 0..n_beta {
            for av in 0..nvb {
                let a = n_beta + av;
                rhs[(nova + i * nvb + av, coordinate)] =
                    skel_b_mo[coordinate][(a, i)] / (eps_b[i] - eps_b[a]);
            }
        }
    }
    for col in 0..nov {
        let mut ua = vec![0.0; nova];
        let mut ub = vec![0.0; novb];
        if col < nova {
            ua[col] = 1.0;
        } else {
            ub[col - nova] = 1.0;
        }
        let dpa = rotation_density(ca, n_alpha, &ua, 1.0);
        let dpb = rotation_density(cb, n_beta, &ub, 1.0);
        let (dfa, dfb) = response_fock(&dpa, &dpb)?;
        let dfa_mo = transform_to_mo(ca, &dfa);
        let dfb_mo = transform_to_mo(cb, &dfb);
        for i in 0..n_alpha {
            for av in 0..nva {
                let row = i * nva + av;
                let a = n_alpha + av;
                system[(row, col)] -= dfa_mo[(a, i)] / (eps_a[i] - eps_a[a]);
            }
        }
        for i in 0..n_beta {
            for av in 0..nvb {
                let row = nova + i * nvb + av;
                let a = n_beta + av;
                system[(row, col)] -= dfb_mo[(a, i)] / (eps_b[i] - eps_b[a]);
            }
        }
    }
    let amplitudes = solve_linear_matrix(&system, &rhs)?;
    (0..skeleton_alpha.len())
        .map(|coordinate| {
            let alpha: Vec<f64> = (0..nova)
                .map(|index| amplitudes[(index, coordinate)])
                .collect();
            let beta: Vec<f64> = (0..novb)
                .map(|index| amplitudes[(nova + index, coordinate)])
                .collect();
            let density_alpha = rotation_density(ca, n_alpha, &alpha, 1.0);
            let density_beta = rotation_density(cb, n_beta, &beta, 1.0);
            let (induced_alpha, induced_beta) = response_fock(&density_alpha, &density_beta)?;
            let fock_alpha = add_matrix(&skeleton_alpha[coordinate], &induced_alpha);
            let fock_beta = add_matrix(&skeleton_beta[coordinate], &induced_beta);
            let fock_alpha_mo = transform_to_mo(ca, &fock_alpha);
            let fock_beta_mo = transform_to_mo(cb, &fock_beta);
            Ok(UhfResponse {
                density_alpha,
                density_beta,
                mo_coeff_alpha: coefficient_derivative(ca, eps_a, &fock_alpha_mo)?,
                mo_coeff_beta: coefficient_derivative(cb, eps_b, &fock_beta_mo)?,
                mo_energies_alpha: (0..ca.cols)
                    .map(|orbital| fock_alpha_mo[(orbital, orbital)])
                    .collect(),
                mo_energies_beta: (0..cb.cols)
                    .map(|orbital| fock_beta_mo[(orbital, orbital)])
                    .collect(),
                fock_alpha_mo,
                fock_beta_mo,
            })
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn second_orbital_response(
    c: &Matrix,
    eps: &[f64],
    cx: &Matrix,
    cy: &Matrix,
    eps_x: &[f64],
    eps_y: &[f64],
    fock_x_mo: &Matrix,
    fock_y_mo: &Matrix,
    fock_xy_mo: &Matrix,
) -> (Matrix, Vec<f64>) {
    let ux = c.transpose_matmul_seq(cx);
    let uy = c.transpose_matmul_seq(cy);
    let mut rotation = Matrix::zeros(c.cols, c.cols);
    let mut energies = vec![0.0; c.cols];
    for p in 0..c.cols {
        for q in 0..c.cols {
            if p == q {
                continue;
            }
            let gap = eps[p] - eps[q];
            if gap.abs() < 1.0e-10 {
                continue;
            }
            let mut numerator = fock_xy_mo[(q, p)] - eps_x[p] * uy[(q, p)] - eps_y[p] * ux[(q, p)];
            for r in 0..c.cols {
                numerator += fock_x_mo[(q, r)] * uy[(r, p)] + fock_y_mo[(q, r)] * ux[(r, p)];
            }
            rotation[(q, p)] = numerator / gap;
        }
        rotation[(p, p)] = -(0..c.cols).map(|q| ux[(q, p)] * uy[(q, p)]).sum::<f64>();
        energies[p] = fock_xy_mo[(p, p)]
            + (0..c.cols)
                .map(|q| fock_x_mo[(p, q)] * uy[(q, p)] + fock_y_mo[(p, q)] * ux[(q, p)])
                .sum::<f64>();
    }
    (c.matmul_seq(&rotation), energies)
}

/// Solve the coupled mixed second-order UHF response without displaced
/// geometries. The alpha and beta density maps are iterated together.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_uhf_second_response<F>(
    ca: &Matrix,
    cb: &Matrix,
    eps_a: &[f64],
    eps_b: &[f64],
    n_alpha: usize,
    n_beta: usize,
    first_x: &UhfResponse,
    first_y: &UhfResponse,
    known_alpha_xy: &Matrix,
    known_beta_xy: &Matrix,
    response_fock: F,
) -> Result<UhfSecondResponse>
where
    F: Fn(&Matrix, &Matrix) -> Result<(Matrix, Matrix)>,
{
    let zero_a = Matrix::zeros(ca.rows, ca.cols);
    let zero_b = Matrix::zeros(cb.rows, cb.cols);
    let mut density_alpha = mixed_density(
        ca,
        &first_x.mo_coeff_alpha,
        &first_y.mo_coeff_alpha,
        &zero_a,
        n_alpha,
        1.0,
    );
    let mut density_beta = mixed_density(
        cb,
        &first_x.mo_coeff_beta,
        &first_y.mo_coeff_beta,
        &zero_b,
        n_beta,
        1.0,
    );
    let mut history_alpha = Vec::new();
    let mut history_beta = Vec::new();
    let mut errors_alpha = Vec::new();
    let mut errors_beta = Vec::new();
    const MAX_DIIS: usize = 8;
    for iteration in 0..500 {
        let (induced_alpha, induced_beta) = response_fock(&density_alpha, &density_beta)?;
        let fock_alpha = add_matrix(known_alpha_xy, &induced_alpha);
        let fock_beta = add_matrix(known_beta_xy, &induced_beta);
        let fock_alpha_mo = transform_to_mo(ca, &fock_alpha);
        let fock_beta_mo = transform_to_mo(cb, &fock_beta);
        let (coeff_alpha, energies_alpha) = second_orbital_response(
            ca,
            eps_a,
            &first_x.mo_coeff_alpha,
            &first_y.mo_coeff_alpha,
            &first_x.mo_energies_alpha,
            &first_y.mo_energies_alpha,
            &first_x.fock_alpha_mo,
            &first_y.fock_alpha_mo,
            &fock_alpha_mo,
        );
        let (coeff_beta, energies_beta) = second_orbital_response(
            cb,
            eps_b,
            &first_x.mo_coeff_beta,
            &first_y.mo_coeff_beta,
            &first_x.mo_energies_beta,
            &first_y.mo_energies_beta,
            &first_x.fock_beta_mo,
            &first_y.fock_beta_mo,
            &fock_beta_mo,
        );
        let next_alpha = mixed_density(
            ca,
            &first_x.mo_coeff_alpha,
            &first_y.mo_coeff_alpha,
            &coeff_alpha,
            n_alpha,
            1.0,
        );
        let next_beta = mixed_density(
            cb,
            &first_x.mo_coeff_beta,
            &first_y.mo_coeff_beta,
            &coeff_beta,
            n_beta,
            1.0,
        );
        let error = next_alpha
            .rms_difference(&density_alpha)
            .max(next_beta.rms_difference(&density_beta));
        if error < 1.0e-11 {
            return Ok(UhfSecondResponse {
                mo_coeff_alpha: coeff_alpha,
                mo_coeff_beta: coeff_beta,
                mo_energies_alpha: energies_alpha,
                mo_energies_beta: energies_beta,
            });
        }
        let mut error_alpha = next_alpha.clone();
        let mut error_beta = next_beta.clone();
        for index in 0..error_alpha.as_slice().len() {
            error_alpha.as_mut_slice()[index] -= density_alpha.as_slice()[index];
            error_beta.as_mut_slice()[index] -= density_beta.as_slice()[index];
        }
        history_alpha.push(next_alpha.clone());
        history_beta.push(next_beta.clone());
        errors_alpha.push(error_alpha);
        errors_beta.push(error_beta);
        if history_alpha.len() > MAX_DIIS {
            history_alpha.remove(0);
            history_beta.remove(0);
            errors_alpha.remove(0);
            errors_beta.remove(0);
        }
        if let Some((alpha, beta)) =
            paired_response_diis(&history_alpha, &history_beta, &errors_alpha, &errors_beta)
        {
            density_alpha = alpha;
            density_beta = beta;
        } else {
            for index in 0..density_alpha.as_slice().len() {
                density_alpha.as_mut_slice()[index] =
                    0.45 * next_alpha.as_slice()[index] + 0.55 * density_alpha.as_slice()[index];
                density_beta.as_mut_slice()[index] =
                    0.45 * next_beta.as_slice()[index] + 0.55 * density_beta.as_slice()[index];
            }
        }
        if iteration == 499 {
            return Err(XndoError::LinearAlgebra(format!(
                "second-order UHF response did not converge (density RMS {error:.3e})"
            )));
        }
    }
    unreachable!()
}
