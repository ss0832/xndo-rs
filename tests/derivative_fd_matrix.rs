// SPDX-License-Identifier: GPL-3.0-or-later

//! Full finite-difference validation matrix for the derivative-capable xndo-rs methods.
//!
//! These tests deliberately compare *every Cartesian component* of the production
//! analytic gradient with full-SCF central differences of the total energy, and
//! *every Hessian matrix element* with central differences of the production
//! analytic gradient.  They are ignored in the default developer loop because
//! each Hessian test requires many SCF solves; CI explicitly executes them with
//! `cargo test --test derivative_fd_matrix -- --ignored --nocapture`.
//!
//! Open-shell MNDO/d systems containing d AOs use the spin-coupled analytic
//! UCPHF/spd-AD path and are checked separately with an SH radical.

use xndo_rs::{
    analytic_hessian, closed_form_gradient, numerical_gradient, numerical_hessian, Method,
    Molecule, NddoOptions, NddoParameters, Reference,
};

#[derive(Clone, Copy)]
struct Case {
    label: &'static str,
    method: Method,
    xyz: &'static str,
    charge: f64,
    multiplicity: usize,
    reference: Reference,
    grad_max_tol: f64,
    grad_rms_tol: f64,
    hess_max_tol: f64,
    hess_rms_tol: f64,
}

fn opts(case: Case) -> NddoOptions {
    NddoOptions {
        charge: case.charge,
        multiplicity: case.multiplicity,
        reference: case.reference,
        // Tighten the ordinary SCF before the energy finite difference as well;
        // analytic_hessian() tightens its own copy internally.
        e_tol: 1.0e-10,
        p_tol: 1.0e-9,
        max_scf: 600,
        ..NddoOptions::default()
    }
}

fn molecule(case: Case) -> Molecule {
    Molecule::from_xyz_str(case.xyz, case.charge)
        .unwrap_or_else(|e| panic!("{}: invalid XYZ: {e}", case.label))
        .with_multiplicity(case.multiplicity)
}

fn gradient_metrics(a: &[xndo_rs::Vec3], b: &[xndo_rs::Vec3]) -> (f64, f64) {
    assert_eq!(a.len(), b.len());
    let mut max_abs = 0.0_f64;
    let mut ss = 0.0_f64;
    let mut n = 0usize;
    for (ga, gb) in a.iter().zip(b) {
        for k in 0..3 {
            let d = ga.get(k) - gb.get(k);
            max_abs = max_abs.max(d.abs());
            ss += d * d;
            n += 1;
        }
    }
    (max_abs, (ss / n.max(1) as f64).sqrt())
}

fn matrix_metrics(a: &xndo_rs::Matrix, b: &xndo_rs::Matrix) -> (f64, f64, f64) {
    assert_eq!((a.rows, a.cols), (b.rows, b.cols));
    let mut max_abs = 0.0_f64;
    let mut max_asym = 0.0_f64;
    let mut ss = 0.0_f64;
    let mut n = 0usize;
    for i in 0..a.rows {
        for j in 0..a.cols {
            let d = a[(i, j)] - b[(i, j)];
            max_abs = max_abs.max(d.abs());
            max_asym = max_asym.max((a[(i, j)] - a[(j, i)]).abs());
            ss += d * d;
            n += 1;
        }
    }
    (max_abs, (ss / n.max(1) as f64).sqrt(), max_asym)
}

fn validate_case(case: Case) {
    let mol = molecule(case);
    let params = NddoParameters::for_method(case.method).unwrap();
    let options = opts(case);

    // Independent first derivative gate: analytic/fixed-density gradient versus
    // central finite difference of fully reconverged total energies.
    let ga = closed_form_gradient(&mol, &params, &options).unwrap();
    let gn = numerical_gradient(&mol, &params, &options, 5.0e-4).unwrap();
    let (gmax, grms) = gradient_metrics(&ga.gradient, &gn.gradient);
    eprintln!(
        "{} gradient: max|analytic-energyFD|={:.6e}, RMS={:.6e} eV/Bohr",
        case.label, gmax, grms
    );
    assert!(
        gmax <= case.grad_max_tol && grms <= case.grad_rms_tol,
        "{} gradient FD gate failed: max={:.6e} (tol {:.3e}), rms={:.6e} (tol {:.3e})",
        case.label,
        gmax,
        case.grad_max_tol,
        grms,
        case.grad_rms_tol
    );

    // Independent second derivative gate: analytic CPHF/AD Hessian versus
    // central finite difference of analytic gradients. Compare the *entire*
    // 3N x 3N matrix, not selected elements or frequencies only.
    let ha = analytic_hessian(&mol, &params, &options, 1.0e-3).unwrap();
    let hn = numerical_hessian(&mol, &params, &options, 1.0e-3).unwrap();
    let (hmax, hrms, hasym) = matrix_metrics(&ha, &hn);
    eprintln!(
        "{} Hessian: max|analytic-gradFD|={:.6e}, RMS={:.6e}, max asym={:.6e} eV/Bohr^2",
        case.label, hmax, hrms, hasym
    );
    assert!(
        hasym <= 1.0e-8,
        "{} analytic Hessian is asymmetric: {:.6e}",
        case.label,
        hasym
    );
    assert!(
        hmax <= case.hess_max_tol && hrms <= case.hess_rms_tol,
        "{} Hessian FD gate failed: max={:.6e} (tol {:.3e}), rms={:.6e} (tol {:.3e})",
        case.label,
        hmax,
        case.hess_max_tol,
        hrms,
        case.hess_rms_tol
    );
}

const CASES: &[Case] = &[
    Case {
        label: "MNDO/RHF/H2O",
        method: Method::Mndo,
        xyz: "3\nwater\nO 0.000000 0.000000 0.000000\nH 0.958400 0.000000 0.000000\nH -0.240000 0.927800 0.000000\n",
        charge: 0.0, multiplicity: 1, reference: Reference::Rhf,
        grad_max_tol: 2.0e-4, grad_rms_tol: 8.0e-5,
        hess_max_tol: 2.0e-3, hess_rms_tol: 6.0e-4,
    },
    Case {
        label: "MNDO/RHF/NH3",
        method: Method::Mndo,
        xyz: "4\nammonia\nN 0.000000 0.000000 0.120000\nH 0.940000 0.000000 -0.280000\nH -0.470000 0.814000 -0.280000\nH -0.470000 -0.814000 -0.280000\n",
        charge: 0.0, multiplicity: 1, reference: Reference::Rhf,
        grad_max_tol: 3.0e-4, grad_rms_tol: 1.0e-4,
        hess_max_tol: 3.0e-3, hess_rms_tol: 8.0e-4,
    },
    Case {
        label: "MNDO/RHF/H2CO",
        method: Method::Mndo,
        xyz: "4\nformaldehyde\nC 0.000000 0.000000 0.000000\nO 1.205000 0.040000 0.030000\nH -0.610000 0.920000 -0.040000\nH -0.570000 -0.950000 0.070000\n",
        charge: 0.0, multiplicity: 1, reference: Reference::Rhf,
        grad_max_tol: 4.0e-4, grad_rms_tol: 1.5e-4,
        hess_max_tol: 4.0e-3, hess_rms_tol: 1.0e-3,
    },
    Case {
        label: "MNDO/UHF/CH3",
        method: Method::Mndo,
        xyz: "4\nmethyl radical\nC 0.000000 0.000000 0.050000\nH 1.070000 0.000000 -0.020000\nH -0.535000 0.927000 0.010000\nH -0.515000 -0.942000 -0.040000\n",
        charge: 0.0, multiplicity: 2, reference: Reference::Uhf,
        grad_max_tol: 8.0e-4, grad_rms_tol: 3.0e-4,
        hess_max_tol: 5.0e-3, hess_rms_tol: 1.5e-3,
    },
    Case {
        label: "MNDO/UHF/OH",
        method: Method::Mndo,
        xyz: "2\nhydroxyl radical\nO 0.000000 0.000000 0.000000\nH 0.965000 0.070000 -0.030000\n",
        charge: 0.0, multiplicity: 2, reference: Reference::Uhf,
        grad_max_tol: 8.0e-4, grad_rms_tol: 3.0e-4,
        hess_max_tol: 5.0e-3, hess_rms_tol: 1.5e-3,
    },
    // MNDO itself also contains explicit 9-AO d-shell parameterizations for
    // transition metals. Cheap off-axis M-H release cases exercise a broad,
    // SCF-stable representative set rather than only main-group s/p molecules.
    // Charges are chosen only to give an even
    // valence-electron count for the RHF derivative gate; these are numerical
    // derivative regressions, not claims about preferred chemical states.
    Case { label: "MNDO/RHF/ScH", method: Method::Mndo, xyz: "2\nScH derivative gate\nSc 0.020 -0.030 0.010\nH 1.420 0.610 -0.430\n", charge: 0.0, multiplicity: 1, reference: Reference::Rhf, grad_max_tol: 3.0e-3, grad_rms_tol: 1.2e-3, hess_max_tol: 3.0e-2, hess_rms_tol: 1.0e-2 },
    Case { label: "MNDO/RHF/TiH+", method: Method::Mndo, xyz: "2\nTiH+ derivative gate\nTi 0.020 -0.030 0.010\nH 1.430 0.620 -0.420\n", charge: 1.0, multiplicity: 1, reference: Reference::Rhf, grad_max_tol: 3.0e-3, grad_rms_tol: 1.2e-3, hess_max_tol: 3.0e-2, hess_rms_tol: 1.0e-2 },
    Case { label: "MNDO/RHF/FeH-", method: Method::Mndo, xyz: "2\nFeH- derivative gate\nFe 0.020 -0.030 0.010\nH 1.620 0.170 0.110\n", charge: -1.0, multiplicity: 1, reference: Reference::Rhf, grad_max_tol: 3.0e-3, grad_rms_tol: 1.2e-3, hess_max_tol: 3.0e-2, hess_rms_tol: 1.0e-2 },
    Case { label: "MNDO/RHF/CoH", method: Method::Mndo, xyz: "2\nCoH derivative gate\nCo 0.020 -0.030 0.010\nH 1.470 0.660 -0.380\n", charge: 0.0, multiplicity: 1, reference: Reference::Rhf, grad_max_tol: 3.0e-3, grad_rms_tol: 1.2e-3, hess_max_tol: 3.0e-2, hess_rms_tol: 1.0e-2 },
    Case { label: "MNDO/RHF/CuH", method: Method::Mndo, xyz: "2\nCuH derivative gate\nCu 0.020 -0.030 0.010\nH 1.490 0.680 -0.360\n", charge: 0.0, multiplicity: 1, reference: Reference::Rhf, grad_max_tol: 3.0e-3, grad_rms_tol: 1.2e-3, hess_max_tol: 3.0e-2, hess_rms_tol: 1.0e-2 },
    Case { label: "MNDO/RHF/ZrH+", method: Method::Mndo, xyz: "2\nZrH+ derivative gate\nZr 0.020 -0.030 0.010\nH 1.570 0.700 -0.350\n", charge: 1.0, multiplicity: 1, reference: Reference::Rhf, grad_max_tol: 4.0e-3, grad_rms_tol: 1.5e-3, hess_max_tol: 4.0e-2, hess_rms_tol: 1.3e-2 },
    Case { label: "MNDO/RHF/MoH3+", method: Method::Mndo, xyz: "2\nMoH3+ derivative gate\nMo 0.020 -0.030 0.010\nH 1.580 0.710 -0.340\n", charge: 3.0, multiplicity: 1, reference: Reference::Rhf, grad_max_tol: 4.0e-3, grad_rms_tol: 1.5e-3, hess_max_tol: 4.0e-2, hess_rms_tol: 1.3e-2 },
    Case { label: "MNDO/RHF/AgH", method: Method::Mndo, xyz: "2\nAgH derivative gate\nAg 0.020 -0.030 0.010\nH 1.600 0.730 -0.320\n", charge: 0.0, multiplicity: 1, reference: Reference::Rhf, grad_max_tol: 4.0e-3, grad_rms_tol: 1.5e-3, hess_max_tol: 4.0e-2, hess_rms_tol: 1.3e-2 },
    Case { label: "MNDO/RHF/PtH+", method: Method::Mndo, xyz: "2\nPtH+ derivative gate\nPt 0.020 -0.030 0.010\nH 1.610 0.740 -0.310\n", charge: 1.0, multiplicity: 1, reference: Reference::Rhf, grad_max_tol: 4.0e-3, grad_rms_tol: 1.5e-3, hess_max_tol: 4.0e-2, hess_rms_tol: 1.3e-2 },
    // Closed-shell MNDO/d cases deliberately cover every element in the
    // built-in MNDO/d table that carries an explicit d exponent: Al, Si, P,
    // S, Cl, Br and I.  This guards against validating only one convenient
    // d-bearing molecule while leaving element-specific spd parameter paths
    // unexercised.
    Case {
        label: "MNDO-d/RHF/AlH3",
        method: Method::MndoD,
        xyz: "4\naluminum trihydride\nAl 0.020000 -0.030000 0.010000\nH 1.550000 0.120000 0.080000\nH -0.850000 1.310000 -0.090000\nH -0.760000 -1.360000 0.140000\n",
        charge: 0.0, multiplicity: 1, reference: Reference::Rhf,
        grad_max_tol: 1.5e-3, grad_rms_tol: 6.0e-4,
        hess_max_tol: 1.0e-2, hess_rms_tol: 3.0e-3,
    },
    Case {
        label: "MNDO-d/RHF/SiH4",
        method: Method::MndoD,
        xyz: "5\nsilane\nSi 0.030000 -0.020000 0.010000\nH 0.850000 0.870000 0.820000\nH -0.870000 -0.820000 0.890000\nH -0.810000 0.900000 -0.850000\nH 0.890000 -0.850000 -0.830000\n",
        charge: 0.0, multiplicity: 1, reference: Reference::Rhf,
        grad_max_tol: 8.0e-4, grad_rms_tol: 3.0e-4,
        hess_max_tol: 5.0e-3, hess_rms_tol: 1.5e-3,
    },
    Case {
        label: "MNDO-d/RHF/PH3",
        method: Method::MndoD,
        xyz: "4\nphosphine\nP 0.020000 -0.010000 0.180000\nH 1.390000 0.100000 -0.360000\nH -0.770000 1.160000 -0.310000\nH -0.680000 -1.250000 -0.410000\n",
        charge: 0.0, multiplicity: 1, reference: Reference::Rhf,
        grad_max_tol: 1.5e-3, grad_rms_tol: 6.0e-4,
        hess_max_tol: 1.0e-2, hess_rms_tol: 3.0e-3,
    },
    Case {
        label: "MNDO-d/RHF/H2S",
        method: Method::MndoD,
        xyz: "3\nhydrogen sulfide\nS 0.010000 -0.020000 0.030000\nH 1.330000 0.130000 -0.060000\nH -0.390000 1.270000 0.160000\n",
        charge: 0.0, multiplicity: 1, reference: Reference::Rhf,
        grad_max_tol: 1.0e-3, grad_rms_tol: 4.0e-4,
        hess_max_tol: 7.0e-3, hess_rms_tol: 2.0e-3,
    },
    Case {
        label: "MNDO-d/RHF/HCl",
        method: Method::MndoD,
        xyz: "2\nhydrogen chloride off-axis\nH -0.210000 0.080000 -0.100000\nCl 1.080000 0.470000 0.380000\n",
        charge: 0.0, multiplicity: 1, reference: Reference::Rhf,
        grad_max_tol: 1.5e-3, grad_rms_tol: 6.0e-4,
        hess_max_tol: 1.0e-2, hess_rms_tol: 3.0e-3,
    },
    Case {
        label: "MNDO-d/RHF/HBr",
        method: Method::MndoD,
        xyz: "2\nhydrogen bromide off-axis\nH -0.180000 0.090000 -0.070000\nBr 1.190000 0.510000 0.420000\n",
        charge: 0.0, multiplicity: 1, reference: Reference::Rhf,
        grad_max_tol: 2.0e-3, grad_rms_tol: 8.0e-4,
        hess_max_tol: 1.2e-2, hess_rms_tol: 4.0e-3,
    },
    Case {
        label: "MNDO-d/RHF/HI",
        method: Method::MndoD,
        xyz: "2\nhydrogen iodide off-axis\nH -0.160000 0.110000 -0.060000\nI 1.330000 0.570000 0.460000\n",
        charge: 0.0, multiplicity: 1, reference: Reference::Rhf,
        grad_max_tol: 2.5e-3, grad_rms_tol: 1.0e-3,
        hess_max_tol: 1.5e-2, hess_rms_tol: 5.0e-3,
    },
    // This uses the MNDO/d parameter set without a d-bearing atom and complements
    // the separate SH test of the analytic open-shell spd/UCPHF path.
    Case {
        label: "MNDO-d/UHF/CH3-s-p",
        method: Method::MndoD,
        xyz: "4\nmethyl radical\nC 0.000000 0.000000 0.050000\nH 1.070000 0.000000 -0.020000\nH -0.535000 0.927000 0.010000\nH -0.515000 -0.942000 -0.040000\n",
        charge: 0.0, multiplicity: 2, reference: Reference::Uhf,
        grad_max_tol: 8.0e-4, grad_rms_tol: 3.0e-4,
        hess_max_tol: 5.0e-3, hess_rms_tol: 1.5e-3,
    },
];

#[test]
fn mndo_release_matrix_covers_the_declared_representative_d_elements() {
    use std::collections::BTreeSet;

    let params = NddoParameters::for_method(Method::Mndo).unwrap();
    let expected: BTreeSet<u8> = [21, 22, 26, 27, 29, 40, 42, 47, 78].into_iter().collect();
    let mut covered = BTreeSet::new();
    for &case in CASES
        .iter()
        .filter(|c| c.method == Method::Mndo && c.reference == Reference::Rhf)
    {
        let mol = molecule(case);
        for atom in &mol.atoms {
            if atom.z <= 97 && params.element(atom.z).unwrap().has_d() {
                covered.insert(atom.z);
            }
        }
    }
    assert_eq!(
        covered, expected,
        "MNDO FD release matrix must cover the declared representative d-element set"
    );
}

#[test]
fn mndod_release_matrix_covers_every_explicit_d_element() {
    use std::collections::BTreeSet;

    let params = NddoParameters::for_method(Method::MndoD).unwrap();
    let expected: BTreeSet<u8> = params
        .elements
        .iter()
        .filter_map(|(&z, e)| if e.has_d() { Some(z) } else { None })
        .collect();
    let mut covered = BTreeSet::new();
    for &case in CASES
        .iter()
        .filter(|c| c.method == Method::MndoD && c.reference == Reference::Rhf)
    {
        let mol = molecule(case);
        for atom in &mol.atoms {
            if params.element(atom.z).unwrap().has_d() {
                covered.insert(atom.z);
            }
        }
    }
    assert_eq!(
        covered, expected,
        "MNDO/d FD release matrix must cover every explicit d-bearing element"
    );
}

#[test]
#[ignore = "full finite-difference derivative matrix; run explicitly in CI/release validation"]
fn full_mndo_mndod_rhf_uhf_derivative_matrix() {
    let filter = std::env::var("XNDO_FD_CASE").ok();
    let mut matched = 0usize;
    for &case in CASES.iter().filter(|case| {
        filter
            .as_ref()
            .is_none_or(|needle| case.label.contains(needle))
    }) {
        matched += 1;
        validate_case(case);
    }
    assert!(matched > 0, "XNDO_FD_CASE did not match a validation case");
}

#[test]
#[ignore = "full open-shell d-AO analytic-Hessian finite-difference validation"]
fn mndod_open_shell_d_analytic_hessian_matches_gradient_differences() {
    let case = Case {
        label: "MNDO-d/UHF/SH-d-analytic",
        method: Method::MndoD,
        xyz: "2\nSH radical\nS 0.010000 -0.020000 0.030000\nH 1.340000 0.120000 -0.070000\n",
        charge: 0.0,
        multiplicity: 2,
        reference: Reference::Uhf,
        grad_max_tol: 1.5e-3,
        grad_rms_tol: 6.0e-4,
        hess_max_tol: 1.0e-4,
        hess_rms_tol: 4.0e-5,
    };
    validate_case(case);
}
