// SPDX-License-Identifier: GPL-3.0-or-later

use xndo_rs::{
    CalculationResult, Method, MethodStatus, Molecule, NddoOptions, NddoParameters, ZindoOptions,
    ZindoParameters,
};

#[test]
fn method_aliases_are_stable() {
    assert_eq!(Method::parse("MNDO"), Some(Method::Mndo));
    assert_eq!(Method::parse("MNDO/d"), Some(Method::MndoD));
    assert_eq!(Method::parse("MNDOD"), Some(Method::MndoD));
    assert_eq!(Method::parse("INDO"), Some(Method::Indo));
    assert_eq!(Method::parse("INDO/S"), Some(Method::ZindoS));
    assert_eq!(Method::parse("ZINDO/S"), Some(Method::ZindoS));
}

#[test]
fn mndo_and_mndod_parameter_sets_load() {
    let mndo = NddoParameters::mndo().unwrap();
    let mndod = NddoParameters::mndod().unwrap();
    assert_eq!(mndo.method, Method::Mndo);
    assert_eq!(mndod.method, Method::MndoD);
    assert!((mndo.element(1).unwrap().u_ss + 11.906276).abs() < 1.0e-9);
    assert!(mndod.element(14).unwrap().has_d());
    assert!(!mndod.element(6).unwrap().has_d());
}

#[test]
fn mndo_water_smoke() {
    let mol = Molecule::from_xyz_str("3\nwater\nO 0 0 0\nH 0.9584 0 0\nH -0.2400 0.9278 0\n", 0.0)
        .unwrap();
    let p = NddoParameters::mndo().unwrap();
    let r = xndo_rs::run_nddo_with_parameters(&mol, &p, &NddoOptions::default()).unwrap();
    assert!(r.total_ev.is_finite());
    assert!(r.converged);
}

#[test]
fn mndod_sih4_smoke() {
    let mol = Molecule::from_xyz_str(
        "5\nsilane\nSi 0 0 0\nH 0.85 0.85 0.85\nH -0.85 -0.85 0.85\nH -0.85 0.85 -0.85\nH 0.85 -0.85 -0.85\n",
        0.0,
    ).unwrap();
    let p = NddoParameters::mndod().unwrap();
    let r = xndo_rs::run_nddo_with_parameters(&mol, &p, &NddoOptions::default()).unwrap();
    assert!(r.total_ev.is_finite());
    assert!(r.converged);
}

#[test]
fn zindo_parameter_table_and_water_cis_smoke() {
    let p = ZindoParameters::standard().unwrap();
    assert_eq!(p.elements.len(), 38);
    let c = p.element(6).unwrap();
    assert!((c.fg[14] - 6.897842 / 3.0).abs() < 1.0e-12);
    assert!((c.fg[15] - 4.509913 / 25.0).abs() < 1.0e-12);

    let mol = Molecule::from_xyz_str(
        "3\nwater\nO 0 0 0\nH 0.7586 0 0.5043\nH -0.7586 0 0.5043\n",
        0.0,
    )
    .unwrap();
    let o = ZindoOptions {
        n_states: 3,
        ..ZindoOptions::default()
    };
    let r = xndo_rs::zindo_s_cis(&mol, &p, &o).unwrap();
    assert!(r.ground.total_ev.is_finite());
    assert!(r.states.iter().all(|s| s.energy_ev.is_finite()));
}

#[test]
fn full_legacy_registry_is_stable() {
    assert_eq!(Method::ALL.len(), 18);
    for name in [
        "cndo1", "cndo2", "cndo3", "indo", "indo1", "indo2", "indo3", "zindo1", "zindo2", "zindo3",
        "zindo/s", "mindo1", "mindo2", "mindo3", "sindo", "msindo", "mndo", "mndod",
    ] {
        assert!(
            Method::parse(name).is_some(),
            "missing registry method {name}"
        );
    }
    assert_eq!(Method::Cndo2.status(), MethodStatus::Native);
    assert_eq!(Method::Indo.status(), MethodStatus::Native);
    assert_eq!(Method::Mindo3.status(), MethodStatus::Native);
    assert_eq!(Method::Cndo3.status(), MethodStatus::NonCanonical);
    assert!(Method::Cndo2.supports_energy());
    assert!(Method::Indo.supports_energy());
}

#[test]
fn cndo2_and_indo_single_points_are_independently_wired() {
    use xndo_rs::{run_cndo_indo, CndoIndoOptions};
    let water =
        Molecule::from_xyz_str("3\nwater\nO 0 0 0\nH 0.9584 0 0\nH -0.2400 0.9278 0\n", 0.0)
            .unwrap();
    for method in [Method::Cndo2, Method::Indo] {
        let r = run_cndo_indo(&water, method, &CndoIndoOptions::default()).unwrap();
        assert!(r.converged);
        assert!(r.total_ev.is_finite());
        assert_eq!(r.method, method);
        assert_eq!(r.charges.len(), 3);
        assert!(r.charges.iter().sum::<f64>().abs() < 1.0e-7);
    }
}

#[test]
fn cndo2_and_indo_high_level_dispatch_are_native() {
    use xndo_rs::{run_method, CalculationResult};
    let h2 = Molecule::from_xyz_str("2\nh2\nH 0 0 -0.37\nH 0 0 0.37\n", 0.0).unwrap();
    for method in [Method::Cndo2, Method::Indo] {
        match run_method(&h2, method, &NddoOptions::default()).unwrap() {
            CalculationResult::CndoIndo(r) => assert_eq!(r.method, method),
            other => panic!("unexpected result variant for {method}: {other:?}"),
        }
    }
}

#[test]
fn indo_does_not_accept_sulfur_without_complete_indo_row() {
    use xndo_rs::{run_cndo_indo, CndoIndoOptions};
    let h2s = Molecule::from_xyz_str("3\nh2s\nS 0 0 0\nH 1.34 0 0\nH -0.45 1.26 0\n", 0.0).unwrap();
    assert!(run_cndo_indo(&h2s, Method::Cndo2, &CndoIndoOptions::default()).is_ok());
    let err = run_cndo_indo(&h2s, Method::Indo, &CndoIndoOptions::default()).unwrap_err();
    assert!(err.to_string().contains("H/Li/C/N/O"));
}

#[test]
fn mindo3_parameter_smoke() {
    let oxygen = xndo_rs::mindo3_element(8).unwrap();
    assert_eq!(oxygen.n_orb, 4);
    assert!((oxygen.uss + 91.73).abs() < 1.0e-12);
    let ch = xndo_rs::mindo3_pair(6, 1).unwrap();
    assert!((ch.0 - 0.315011).abs() < 1.0e-12);
}

#[test]
fn mindo3_uhf_methyl_radical_smoke_and_spin_trace() {
    use xndo_rs::{run_mindo3, Mindo3Options, Reference};
    let mol = Molecule::from_xyz_str(
        "4\nmethyl radical\nC 0 0 0.05\nH 1.07 0 -0.02\nH -0.535 0.927 0.01\nH -0.515 -0.942 -0.04\n",
        0.0,
    ).unwrap().with_multiplicity(2);
    let opts = Mindo3Options {
        multiplicity: 2,
        reference: Reference::Uhf,
        max_scf: 600,
        damping: 0.35,
        ..Mindo3Options::default()
    };
    let r = run_mindo3(&mol, &opts).unwrap();
    assert!(r.converged);
    assert!(r.unrestricted);
    assert_eq!((r.n_alpha, r.n_beta), (4, 3));
    assert!(r.mo_energies_beta_ev.is_some());
    let sd = r.spin_density.as_ref().unwrap();
    let trace: f64 = (0..sd.rows).map(|i| sd[(i, i)]).sum();
    assert!(
        (trace - 1.0).abs() < 1.0e-7,
        "UHF spin-density trace={trace}"
    );
}

#[test]
fn zindo_excited_state_properties_are_normalized_and_spin_resolved() {
    use xndo_rs::{zindo_s_cis_spin, CisSpin};
    let p = ZindoParameters::standard().unwrap();
    let mol = Molecule::from_xyz_str(
        "3\nwater\nO 0 0 0\nH 0.7586 0 0.5043\nH -0.7586 0 0.5043\n",
        0.0,
    )
    .unwrap();
    let o = ZindoOptions {
        n_states: 3,
        ..ZindoOptions::default()
    };

    for spin in [CisSpin::Singlet, CisSpin::Triplet] {
        let spec = zindo_s_cis_spin(&mol, &p, &o, spin).unwrap();
        assert!(!spec.states.is_empty());
        for st in &spec.states {
            assert_eq!(st.spin, spin);
            assert!(st.energy_ev.is_finite() && st.energy_cm1.is_finite());
            assert!(st.state_total_energy_ev.is_finite());
            assert!(st.state_total_energy_hartree.is_finite());
            assert!(
                (st.state_total_energy_ev - (spec.ground.total_ev + st.energy_ev)).abs() < 1.0e-10
            );
            assert!(st.charge_transfer_distance_angstrom.is_finite());
            assert!(st.permanent_dipole_magnitude_debye.is_finite());
            assert!(st.difference_dipole_magnitude_debye.is_finite());
            match spin {
                CisSpin::Singlet => {
                    assert_eq!(st.spin_multiplicity, 1);
                    assert_eq!(st.s2_expectation, 0.0);
                }
                CisSpin::Triplet => {
                    assert_eq!(st.spin_multiplicity, 3);
                    assert_eq!(st.s2_expectation, 2.0);
                }
                CisSpin::Unrestricted => unreachable!("this loop covers RHF-CIS sectors"),
            }
            let qsum: f64 = st.charges.iter().sum();
            assert!(qsum.abs() < 1.0e-7, "excited-state charge sum={qsum}");
            let hsum: f64 = st.hole_population.iter().sum();
            let esum: f64 = st.electron_population.iter().sum();
            assert!((hsum - 1.0).abs() < 1.0e-6, "hole population sum={hsum}");
            assert!(
                (esum - 1.0).abs() < 1.0e-6,
                "electron population sum={esum}"
            );
            if spin == CisSpin::Triplet {
                assert_eq!(st.oscillator_strength, 0.0);
                assert_eq!(st.transition_dipole_au, [0.0; 3]);
            }
        }
    }
}

#[test]
fn zindo_high_level_dispatch_accepts_uhf_reference() {
    use xndo_rs::{run_method, Reference};
    let mol = Molecule::from_xyz_str(
        "3\nwater\nO 0 0 0\nH 0.7586 0 0.5043\nH -0.7586 0 0.5043\n",
        0.0,
    )
    .unwrap();
    let opts = NddoOptions {
        reference: Reference::Uhf,
        ..NddoOptions::default()
    };
    let result = run_method(&mol, Method::ZindoS, &opts).unwrap();
    match result {
        CalculationResult::ZindoS(result) => assert!(result.unrestricted),
        _ => panic!("unexpected high-level ZINDO/S result variant"),
    }
}

#[test]
fn executable_api_and_alias_matrix_is_explicit() {
    assert_eq!(
        Method::Cndo2.api_names(),
        &[
            "single_point",
            "gradient",
            "forces",
            "optimize",
            "hessian",
            "frequencies"
        ]
    );
    assert_eq!(
        Method::Cndo2.accepted_strings(),
        &["cndo2", "cndo/2", "cndo"]
    );
    assert_eq!(Method::Indo.api_names(), Method::Cndo2.api_names());
    assert_eq!(Method::Indo.accepted_strings(), &["indo"]);
    assert_eq!(
        Method::Mndo.api_names(),
        &[
            "single_point",
            "gradient",
            "forces",
            "optimize",
            "hessian",
            "frequencies"
        ]
    );
    assert_eq!(Method::Mndo.accepted_strings(), &["mndo"]);
    assert_eq!(
        Method::MndoD.accepted_strings(),
        &["mndod", "mndo/d", "mndo-d"]
    );
    assert_eq!(Method::Mindo3.api_names(), Method::Cndo2.api_names());
    assert_eq!(
        Method::Mindo3.accepted_strings(),
        &["mindo3", "mindo/3", "mindo"]
    );
    assert_eq!(
        Method::ZindoS.api_names(),
        &[
            "single_point",
            "gradient",
            "forces",
            "optimize",
            "hessian",
            "frequencies",
            "excited_states",
            "excited_properties",
            "uv_vis_spectrum",
            "excited_state_gradients",
            "excited_state_hessians"
        ]
    );
    assert_eq!(
        Method::ZindoS.accepted_strings(),
        &["zindo/s", "zindos", "zindo", "indo/s", "indos"]
    );
    assert!(Method::Cndo2.supports_uhf());
    assert!(Method::Indo.supports_uhf());
    assert!(Method::Mndo.supports_uhf());
    assert!(Method::MndoD.supports_uhf());
    assert!(Method::Mindo3.supports_uhf());
    assert!(Method::ZindoS.supports_uhf());
}

#[test]
fn rust_api_to_method_discovery_is_available() {
    let single = xndo_rs::methods_for_api("single_point");
    assert_eq!(
        single,
        vec![
            Method::Cndo2,
            Method::Indo,
            Method::ZindoS,
            Method::Mindo3,
            Method::Mndo,
            Method::MndoD
        ]
    );
    assert_eq!(
        xndo_rs::methods_for_api("gradient"),
        vec![
            Method::Cndo2,
            Method::Indo,
            Method::ZindoS,
            Method::Mindo3,
            Method::Mndo,
            Method::MndoD
        ]
    );
    assert_eq!(
        xndo_rs::methods_for_api("excited_properties"),
        vec![Method::ZindoS]
    );
    assert_eq!(
        xndo_rs::methods_for_api("uv_vis_spectrum"),
        vec![Method::ZindoS]
    );
    assert!(xndo_rs::methods_for_api("not_an_api").is_empty());
}

#[test]
fn every_advertised_method_string_round_trips() {
    for method in Method::ALL {
        for alias in method.accepted_strings() {
            assert_eq!(
                Method::parse(alias),
                Some(method),
                "advertised alias {alias:?} does not parse as {method}"
            );
        }
    }
}
