// SPDX-License-Identifier: GPL-3.0-or-later

use xndo_rs::{
    zindo_s_cis_gradients, zindo_s_cis_hessians, zindo_s_cis_spin, zindo_s_ucis,
    zindo_s_ucis_gradients, zindo_s_ucis_hessians, Atom, CisSpin, Molecule, Reference, Vec3,
    ZindoOptions, ZindoParameters,
};

fn water() -> Molecule {
    Molecule::new(vec![
        Atom {
            z: 8,
            position: Vec3::new(0.12, -0.08, 0.04),
        },
        Atom {
            z: 1,
            position: Vec3::new(1.58, 0.21, -0.13),
        },
        Atom {
            z: 1,
            position: Vec3::new(-0.37, 1.41, 0.29),
        },
    ])
}

fn methyl() -> Molecule {
    Molecule::new(vec![
        Atom {
            z: 6,
            position: Vec3::new(0.07, -0.02, 0.05),
        },
        Atom {
            z: 1,
            position: Vec3::new(1.72, 0.18, -0.11),
        },
        Atom {
            z: 1,
            position: Vec3::new(-0.71, 1.43, 0.22),
        },
        Atom {
            z: 1,
            position: Vec3::new(-0.63, -1.36, -0.31),
        },
    ])
    .with_multiplicity(2)
}

fn nh2_radical() -> Molecule {
    Molecule::new(vec![
        Atom {
            z: 7,
            position: Vec3::new(0.10, -0.05, 0.03),
        },
        Atom {
            z: 1,
            position: Vec3::new(1.79, 0.23, -0.12),
        },
        Atom {
            z: 1,
            position: Vec3::new(-0.52, 1.64, 0.21),
        },
    ])
    .with_multiplicity(2)
}

fn ho2_radical() -> Molecule {
    Molecule::new(vec![
        Atom {
            z: 8,
            position: Vec3::new(0.10, -0.10, 0.05),
        },
        Atom {
            z: 8,
            position: Vec3::new(2.45, 0.20, -0.15),
        },
        Atom {
            z: 1,
            position: Vec3::new(3.31, 1.55, 0.20),
        },
    ])
    .with_multiplicity(2)
}

fn oh_radical() -> Molecule {
    Molecule::new(vec![
        Atom {
            z: 8,
            position: Vec3::new(0.11, -0.07, 0.03),
        },
        Atom {
            z: 1,
            position: Vec3::new(1.91, 0.24, -0.16),
        },
    ])
    .with_multiplicity(2)
}

fn ucis_options() -> ZindoOptions {
    ZindoOptions {
        multiplicity: 2,
        reference: Reference::Uhf,
        n_states: 1,
        max_scf: 800,
        e_tol_ev: 1.0e-10,
        p_tol: 1.0e-9,
        damping: 0.20,
        ..ZindoOptions::default()
    }
}

#[test]
fn unrestricted_cis_state_gradient_matches_reconverged_energy_differences() {
    let mol = methyl();
    let params = ZindoParameters::cached().unwrap();
    let options = ucis_options();
    let analytic = zindo_s_ucis_gradients(&mol, params, &options).unwrap();
    assert_eq!(analytic.states.len(), 1);
    let step = 2.0e-4;
    let mut max_error = 0.0_f64;
    for coordinate in 0..12 {
        let mut plus = mol.clone();
        let mut minus = mol.clone();
        displace(&mut plus, coordinate, step);
        displace(&mut minus, coordinate, -step);
        let ep = zindo_s_ucis(&plus, params, &options).unwrap();
        let em = zindo_s_ucis(&minus, params, &options).unwrap();
        let finite_difference = (ep.states[0].state_total_energy_ev
            - em.states[0].state_total_energy_ev)
            / (2.0 * step);
        let predicted = component(
            analytic.states[0].state_gradient[coordinate / 3],
            coordinate % 3,
        );
        max_error = max_error.max((predicted - finite_difference).abs());
    }
    eprintln!("ZINDO/S UCIS state-gradient max error={max_error:.3e} eV/Bohr");
    assert!(max_error < 5.0e-4, "UCIS gradient error {max_error:.6e}");
}

#[test]
fn unrestricted_cis_gradients_cover_nh2_and_ho2_radicals() {
    let params = ZindoParameters::cached().unwrap();
    let options = ucis_options();
    let step = 2.0e-4;
    for (name, mol) in [("NH2", nh2_radical()), ("HO2", ho2_radical())] {
        let analytic = zindo_s_ucis_gradients(&mol, params, &options).unwrap();
        assert_eq!(analytic.states.len(), 1);
        let mut max_error = 0.0_f64;
        for coordinate in 0..(3 * mol.atoms.len()) {
            let mut plus = mol.clone();
            let mut minus = mol.clone();
            displace(&mut plus, coordinate, step);
            displace(&mut minus, coordinate, -step);
            let ep = zindo_s_ucis(&plus, params, &options).unwrap();
            let em = zindo_s_ucis(&minus, params, &options).unwrap();
            let finite_difference = (ep.states[0].state_total_energy_ev
                - em.states[0].state_total_energy_ev)
                / (2.0 * step);
            let predicted = component(
                analytic.states[0].state_gradient[coordinate / 3],
                coordinate % 3,
            );
            max_error = max_error.max((predicted - finite_difference).abs());
        }
        eprintln!("ZINDO/S {name} UCIS state-gradient max error={max_error:.3e} eV/Bohr");
        assert!(
            max_error < 6.0e-4,
            "{name} UCIS gradient error {max_error:.6e}"
        );
    }
}

#[test]
fn oh_radical_ucis_spectrum_is_retained_with_explicit_degenerate_response() {
    let mol = oh_radical();
    let params = ZindoParameters::cached().unwrap();
    let options = ucis_options();
    let spectrum = zindo_s_ucis(&mol, params, &options).unwrap();
    assert!(!spectrum.states.is_empty());
    assert_eq!(spectrum.spin, CisSpin::Unrestricted);
    let gradient_error = zindo_s_ucis_gradients(&mol, params, &options).unwrap_err();
    assert!(
        gradient_error.to_string().contains("not unique"),
        "unexpected OH UCIS gradient result: {gradient_error}"
    );
    let error = zindo_s_ucis_hessians(&mol, params, &options).unwrap_err();
    assert!(
        error.to_string().contains("not unique"),
        "unexpected OH UCIS response result: {error}"
    );
}

#[test]
fn near_degenerate_planar_methyl_state_derivatives_are_explicit() {
    let mol = Molecule::from_xyz_str(
        "4\nplanar methyl radical\nC 0 0 0\nH 1.0790 0 0\nH -0.5395 0.9344 0\nH -0.5395 -0.9344 0\n",
        0.0,
    )
    .unwrap()
    .with_multiplicity(2);
    let params = ZindoParameters::cached().unwrap();
    let options = ucis_options();
    assert!(!zindo_s_ucis(&mol, params, &options)
        .unwrap()
        .states
        .is_empty());
    let gradient_error = zindo_s_ucis_gradients(&mol, params, &options).unwrap_err();
    assert!(gradient_error.to_string().contains("not unique"));
    let hessian_error = zindo_s_ucis_hessians(&mol, params, &options).unwrap_err();
    assert!(hessian_error.to_string().contains("not unique"));
}

#[test]
fn unrestricted_cis_state_hessian_matches_analytic_gradient_differences() {
    let mol = methyl();
    let params = ZindoParameters::cached().unwrap();
    let options = ucis_options();
    let analytic = zindo_s_ucis_hessians(&mol, params, &options).unwrap();
    assert_eq!(analytic.states.len(), 1);
    let step = 2.0e-4;
    let mut max_error = 0.0_f64;
    for column in 0..12 {
        let mut plus = mol.clone();
        let mut minus = mol.clone();
        displace(&mut plus, column, step);
        displace(&mut minus, column, -step);
        let gp = zindo_s_ucis_gradients(&plus, params, &options).unwrap();
        let gm = zindo_s_ucis_gradients(&minus, params, &options).unwrap();
        for row in 0..12 {
            let finite_difference = (component(gp.states[0].state_gradient[row / 3], row % 3)
                - component(gm.states[0].state_gradient[row / 3], row % 3))
                / (2.0 * step);
            max_error = max_error
                .max((analytic.states[0].state_hessian[(row, column)] - finite_difference).abs());
        }
    }
    eprintln!("ZINDO/S UCIS state-Hessian max error={max_error:.3e} eV/Bohr2");
    assert!(max_error < 8.0e-4, "UCIS Hessian error {max_error:.6e}");
}

fn displace(mol: &mut Molecule, coordinate: usize, delta: f64) {
    let atom = coordinate / 3;
    match coordinate % 3 {
        0 => mol.atoms[atom].position.x += delta,
        1 => mol.atoms[atom].position.y += delta,
        _ => mol.atoms[atom].position.z += delta,
    }
}

fn component(value: Vec3, axis: usize) -> f64 {
    match axis {
        0 => value.x,
        1 => value.y,
        _ => value.z,
    }
}

#[test]
fn singlet_and_triplet_cis_state_gradients_match_reconverged_energy_differences() {
    let mol = water();
    let params = ZindoParameters::cached().unwrap();
    let options = ZindoOptions {
        n_states: 3,
        max_scf: 800,
        e_tol_ev: 1.0e-10,
        p_tol: 1.0e-9,
        damping: 0.20,
        ..ZindoOptions::default()
    };
    let step = 2.0e-4;
    for spin in [CisSpin::Singlet, CisSpin::Triplet] {
        let analytic = zindo_s_cis_gradients(&mol, params, &options, spin).unwrap();
        assert!(!analytic.states.is_empty());
        let nstates = analytic.states.len().min(2);
        let mut max_error = 0.0_f64;
        for coordinate in 0..9 {
            let mut plus = mol.clone();
            let mut minus = mol.clone();
            displace(&mut plus, coordinate, step);
            displace(&mut minus, coordinate, -step);
            let ep = zindo_s_cis_spin(&plus, params, &options, spin).unwrap();
            let em = zindo_s_cis_spin(&minus, params, &options, spin).unwrap();
            for state in 0..nstates {
                let finite_difference = (ep.states[state].state_total_energy_ev
                    - em.states[state].state_total_energy_ev)
                    / (2.0 * step);
                let predicted = component(
                    analytic.states[state].state_gradient[coordinate / 3],
                    coordinate % 3,
                );
                max_error = max_error.max((predicted - finite_difference).abs());
            }
        }
        eprintln!(
            "ZINDO/S {} CIS state-gradient max error={max_error:.3e} eV/Bohr",
            spin.as_str()
        );
        assert!(max_error < 3.0e-4, "CIS gradient error {max_error:.6e}");
    }
}

#[test]
fn singlet_and_triplet_cis_state_hessians_match_analytic_gradient_differences() {
    let mol = water();
    let params = ZindoParameters::cached().unwrap();
    let options = ZindoOptions {
        n_states: 1,
        max_scf: 800,
        e_tol_ev: 1.0e-10,
        p_tol: 1.0e-9,
        damping: 0.20,
        ..ZindoOptions::default()
    };
    let step = 2.0e-4;
    for spin in [CisSpin::Singlet, CisSpin::Triplet] {
        let analytic = zindo_s_cis_hessians(&mol, params, &options, spin).unwrap();
        assert_eq!(analytic.states.len(), 1);
        let mut max_error = 0.0_f64;
        for column in 0..9 {
            let mut plus = mol.clone();
            let mut minus = mol.clone();
            displace(&mut plus, column, step);
            displace(&mut minus, column, -step);
            let gp = zindo_s_cis_gradients(&plus, params, &options, spin).unwrap();
            let gm = zindo_s_cis_gradients(&minus, params, &options, spin).unwrap();
            for row in 0..9 {
                let finite_difference = (component(gp.states[0].state_gradient[row / 3], row % 3)
                    - component(gm.states[0].state_gradient[row / 3], row % 3))
                    / (2.0 * step);
                max_error = max_error.max(
                    (analytic.states[0].state_hessian[(row, column)] - finite_difference).abs(),
                );
            }
        }
        eprintln!(
            "ZINDO/S {} CIS state-Hessian max error={max_error:.3e} eV/Bohr2",
            spin.as_str()
        );
        assert!(max_error < 5.0e-4, "CIS Hessian error {max_error:.6e}");
    }
}
