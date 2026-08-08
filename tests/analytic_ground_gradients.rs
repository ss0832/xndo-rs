// SPDX-License-Identifier: GPL-3.0-or-later

use xndo_rs::{
    run_gradient, run_method, Atom, CalculationResult, CndoIndoOptions, Method, Mindo3Options,
    Molecule, NddoOptions, Reference, Vec3, ZindoOptions, ZindoParameters,
};

fn molecule(atoms: &[(u8, [f64; 3])], multiplicity: usize) -> Molecule {
    Molecule::new(
        atoms
            .iter()
            .map(|&(z, r)| Atom {
                z,
                position: Vec3::new(r[0], r[1], r[2]),
            })
            .collect(),
    )
    .with_multiplicity(multiplicity)
}

fn options(multiplicity: usize, reference: Reference) -> NddoOptions {
    NddoOptions {
        multiplicity,
        reference,
        max_scf: 800,
        e_tol: 1.0e-10,
        p_tol: 1.0e-9,
        damping: 0.20,
        ..NddoOptions::default()
    }
}

fn energy(mol: &Molecule, method: Method, options: &NddoOptions) -> f64 {
    match run_method(mol, method, options).unwrap() {
        CalculationResult::CndoIndo(r) => r.total_ev,
        CalculationResult::Nddo(r) => r.total_ev,
        CalculationResult::Mindo3(r) => r.total_ev,
        CalculationResult::ZindoS(r) => r.total_ev,
    }
}

fn component(v: Vec3, axis: usize) -> f64 {
    match axis {
        0 => v.x,
        1 => v.y,
        _ => v.z,
    }
}

fn displace(mol: &mut Molecule, atom: usize, axis: usize, delta: f64) {
    match axis {
        0 => mol.atoms[atom].position.x += delta,
        1 => mol.atoms[atom].position.y += delta,
        _ => mol.atoms[atom].position.z += delta,
    }
}

fn assert_full_gradient(mol: &Molecule, method: Method, options: &NddoOptions, tolerance: f64) {
    let analytic = run_gradient(mol, method, options).unwrap();
    let step = 2.0e-4;
    let mut max_error = 0.0_f64;
    for atom in 0..mol.atoms.len() {
        for axis in 0..3 {
            let mut plus = mol.clone();
            let mut minus = mol.clone();
            displace(&mut plus, atom, axis, step);
            displace(&mut minus, atom, axis, -step);
            let finite_difference =
                (energy(&plus, method, options) - energy(&minus, method, options)) / (2.0 * step);
            let error = (component(analytic.gradient[atom], axis) - finite_difference).abs();
            max_error = max_error.max(error);
        }
    }
    let translation = analytic
        .gradient
        .iter()
        .copied()
        .fold(Vec3::zero(), |sum, value| sum + value)
        .norm();
    eprintln!(
        "{} {:?}: max analytic/FD error={:.3e} eV/Bohr, translation={:.3e}",
        method, options.reference, max_error, translation
    );
    assert!(
        max_error < tolerance,
        "{method} gradient error {max_error:.6e}"
    );
    assert!(
        translation < 1.0e-9,
        "{method} translation residual {translation:.6e}"
    );
}

fn cndo_options(options: &NddoOptions) -> CndoIndoOptions {
    CndoIndoOptions {
        charge: options.charge,
        multiplicity: options.multiplicity,
        reference: options.reference,
        max_scf: options.max_scf,
        e_tol_ev: options.e_tol,
        p_tol: options.p_tol,
        damping: options.damping,
    }
}

fn assert_cndo_hessian(mol: &Molecule, method: Method, options: &NddoOptions, tolerance: f64) {
    let (_, _, analytic) =
        xndo_rs::cndo_indo::analytic_ground_hessian(mol, method, &cndo_options(options)).unwrap();
    let step = 2.0e-4;
    let ncoord = 3 * mol.atoms.len();
    let mut max_error = 0.0_f64;
    for column in 0..ncoord {
        let mut plus = mol.clone();
        let mut minus = mol.clone();
        displace(&mut plus, column / 3, column % 3, step);
        displace(&mut minus, column / 3, column % 3, -step);
        let gp = run_gradient(&plus, method, options).unwrap();
        let gm = run_gradient(&minus, method, options).unwrap();
        for row in 0..ncoord {
            let value = (component(gp.gradient[row / 3], row % 3)
                - component(gm.gradient[row / 3], row % 3))
                / (2.0 * step);
            max_error = max_error.max((analytic[(row, column)] - value).abs());
        }
    }
    eprintln!("{method} Hessian max analytic/gradient-FD error={max_error:.3e} eV/Bohr2");
    assert!(
        max_error < tolerance,
        "{method} Hessian error {max_error:.6e}"
    );
}

fn assert_mindo_hessian(mol: &Molecule, options: &NddoOptions, tolerance: f64) {
    let mindo_options = Mindo3Options {
        charge: options.charge,
        multiplicity: options.multiplicity,
        reference: options.reference,
        max_scf: options.max_scf,
        e_tol_ev: options.e_tol,
        p_tol: options.p_tol,
        damping: options.damping,
    };
    let (_, _, analytic) = xndo_rs::mindo3::analytic_ground_hessian(mol, &mindo_options).unwrap();
    let step = 2.0e-4;
    let ncoord = 3 * mol.atoms.len();
    let mut max_error = 0.0_f64;
    for column in 0..ncoord {
        let mut plus = mol.clone();
        let mut minus = mol.clone();
        displace(&mut plus, column / 3, column % 3, step);
        displace(&mut minus, column / 3, column % 3, -step);
        let gp = run_gradient(&plus, Method::Mindo3, options).unwrap();
        let gm = run_gradient(&minus, Method::Mindo3, options).unwrap();
        for row in 0..ncoord {
            let value = (component(gp.gradient[row / 3], row % 3)
                - component(gm.gradient[row / 3], row % 3))
                / (2.0 * step);
            max_error = max_error.max((analytic[(row, column)] - value).abs());
        }
    }
    eprintln!("MINDO/3 Hessian max analytic/gradient-FD error={max_error:.3e} eV/Bohr2");
    assert!(
        max_error < tolerance,
        "MINDO/3 Hessian error {max_error:.6e}"
    );
}

fn assert_zindo_hessian(mol: &Molecule, options: &NddoOptions, tolerance: f64) {
    let params = ZindoParameters::cached().unwrap();
    let zindo_options = ZindoOptions {
        charge: options.charge,
        multiplicity: options.multiplicity,
        max_scf: options.max_scf,
        e_tol_ev: options.e_tol,
        p_tol: options.p_tol,
        damping: options.damping,
        ..ZindoOptions::default()
    };
    let (_, _, analytic) =
        xndo_rs::zindo::analytic_ground_hessian(mol, params, &zindo_options).unwrap();
    let step = 2.0e-4;
    let ncoord = 3 * mol.atoms.len();
    let mut max_error = 0.0_f64;
    let mut max_index = (0usize, 0usize);
    for column in 0..ncoord {
        let mut plus = mol.clone();
        let mut minus = mol.clone();
        displace(&mut plus, column / 3, column % 3, step);
        displace(&mut minus, column / 3, column % 3, -step);
        let gp = run_gradient(&plus, Method::ZindoS, options).unwrap();
        let gm = run_gradient(&minus, Method::ZindoS, options).unwrap();
        for row in 0..ncoord {
            let value = (component(gp.gradient[row / 3], row % 3)
                - component(gm.gradient[row / 3], row % 3))
                / (2.0 * step);
            let error = (analytic[(row, column)] - value).abs();
            if error > max_error {
                max_error = error;
                max_index = (row, column);
            }
        }
    }
    eprintln!(
        "ZINDO/S Hessian max analytic/gradient-FD error={max_error:.3e} eV/Bohr2 at {:?}",
        max_index
    );
    assert!(
        max_error < tolerance,
        "ZINDO/S Hessian error {max_error:.6e}"
    );
}

#[test]
fn native_rhf_gradients_match_reconverged_finite_differences() {
    let water = molecule(
        &[
            (8, [0.12, -0.08, 0.04]),
            (1, [1.58, 0.21, -0.13]),
            (1, [-0.37, 1.41, 0.29]),
        ],
        1,
    );
    let options = options(1, Reference::Rhf);
    for method in [Method::Cndo2, Method::Indo, Method::Mindo3, Method::ZindoS] {
        assert_full_gradient(&water, method, &options, 2.0e-4);
    }
}

#[test]
fn native_uhf_gradients_match_reconverged_finite_differences() {
    let methyl = molecule(
        &[
            (6, [0.07, -0.02, 0.05]),
            (1, [1.72, 0.18, -0.11]),
            (1, [-0.71, 1.43, 0.22]),
            (1, [-0.63, -1.36, -0.31]),
        ],
        2,
    );
    let options = options(2, Reference::Uhf);
    for method in [Method::Cndo2, Method::Indo, Method::Mindo3, Method::ZindoS] {
        assert_full_gradient(&methyl, method, &options, 3.0e-4);
    }
}

#[test]
fn cndo_indo_rhf_and_uhf_hessians_match_gradient_differences() {
    let water = molecule(
        &[
            (8, [0.12, -0.08, 0.04]),
            (1, [1.58, 0.21, -0.13]),
            (1, [-0.37, 1.41, 0.29]),
        ],
        1,
    );
    let rhf = options(1, Reference::Rhf);
    for method in [Method::Cndo2, Method::Indo] {
        assert_cndo_hessian(&water, method, &rhf, 3.0e-4);
    }

    let methyl = molecule(
        &[
            (6, [0.07, -0.02, 0.05]),
            (1, [1.72, 0.18, -0.11]),
            (1, [-0.71, 1.43, 0.22]),
            (1, [-0.63, -1.36, -0.31]),
        ],
        2,
    );
    let uhf = options(2, Reference::Uhf);
    for method in [Method::Cndo2, Method::Indo] {
        assert_cndo_hessian(&methyl, method, &uhf, 5.0e-4);
    }
}

#[test]
fn mindo3_rhf_and_uhf_hessians_match_gradient_differences() {
    let water = molecule(
        &[
            (8, [0.12, -0.08, 0.04]),
            (1, [1.58, 0.21, -0.13]),
            (1, [-0.37, 1.41, 0.29]),
        ],
        1,
    );
    assert_mindo_hessian(&water, &options(1, Reference::Rhf), 3.0e-4);

    let methyl = molecule(
        &[
            (6, [0.07, -0.02, 0.05]),
            (1, [1.72, 0.18, -0.11]),
            (1, [-0.71, 1.43, 0.22]),
            (1, [-0.63, -1.36, -0.31]),
        ],
        2,
    );
    assert_mindo_hessian(&methyl, &options(2, Reference::Uhf), 5.0e-4);
}

#[test]
fn zindo_s_rhf_and_uhf_hessians_match_gradient_differences() {
    let water = molecule(
        &[
            (8, [0.12, -0.08, 0.04]),
            (1, [1.58, 0.21, -0.13]),
            (1, [-0.37, 1.41, 0.29]),
        ],
        1,
    );
    assert_zindo_hessian(&water, &options(1, Reference::Rhf), 3.0e-4);

    let methyl = molecule(
        &[
            (6, [0.07, -0.02, 0.05]),
            (1, [1.72, 0.18, -0.11]),
            (1, [-0.71, 1.43, 0.22]),
            (1, [-0.63, -1.36, -0.31]),
        ],
        2,
    );
    assert_zindo_hessian(&methyl, &options(2, Reference::Uhf), 6.0e-4);
}

#[test]
fn zindo_s_uhf_derivatives_cover_nh2_and_ho2_radicals() {
    let radicals = [
        (
            "NH2",
            molecule(
                &[
                    (7, [0.10, -0.05, 0.03]),
                    (1, [1.79, 0.23, -0.12]),
                    (1, [-0.52, 1.64, 0.21]),
                ],
                2,
            ),
        ),
        (
            "HO2",
            molecule(
                &[
                    (8, [0.10, -0.10, 0.05]),
                    (8, [2.45, 0.20, -0.15]),
                    (1, [3.31, 1.55, 0.20]),
                ],
                2,
            ),
        ),
    ];
    let uhf = options(2, Reference::Uhf);
    for (name, radical) in radicals {
        eprintln!("additional ZINDO/S UHF radical: {name}");
        assert_full_gradient(&radical, Method::ZindoS, &uhf, 4.0e-4);
        assert_zindo_hessian(&radical, &uhf, 7.0e-4);
    }
}

#[test]
fn oh_radical_is_retained_and_degenerate_hessian_is_explicit() {
    let oh = molecule(&[(8, [0.11, -0.07, 0.03]), (1, [1.91, 0.24, -0.16])], 2);
    let uhf = options(2, Reference::Uhf);
    assert_full_gradient(&oh, Method::ZindoS, &uhf, 4.0e-4);
    let zindo_options = ZindoOptions {
        multiplicity: 2,
        reference: Reference::Uhf,
        max_scf: uhf.max_scf,
        e_tol_ev: uhf.e_tol,
        p_tol: uhf.p_tol,
        damping: uhf.damping,
        ..ZindoOptions::default()
    };
    let error = xndo_rs::zindo::analytic_ground_hessian(
        &oh,
        ZindoParameters::cached().unwrap(),
        &zindo_options,
    )
    .unwrap_err();
    assert!(
        error.to_string().contains("not unique"),
        "unexpected OH response result: {error}"
    );
}
