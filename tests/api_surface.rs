// SPDX-License-Identifier: GPL-3.0-or-later

use xndo_rs::{
    methods_for_api, run_nddo, Method, Molecule, NddoCalculator, NddoOptions, NddoParameters,
    Reference, XndoError, ZindoParameters,
};

const WATER: &str =
    "3\nwater\nO 0.0000 0.0000 0.0000\nH 0.9584 0.0000 0.0000\nH -0.2400 0.9278 0.0000\n";

#[test]
fn canonical_nddo_types_are_directly_usable() {
    let params = NddoParameters::mndo().unwrap();
    let options = NddoOptions {
        reference: Reference::Rhf,
        ..NddoOptions::default()
    };
    let calculator = NddoCalculator::with_options(params, options);
    let molecule = Molecule::from_xyz_str(WATER, 0.0).unwrap();
    let result = calculator.calculate(&molecule).unwrap();
    assert!(result.converged);
    assert!(!result.unrestricted);
    assert!(result.total_ev.is_finite());
    assert_eq!(result.charges.len(), 3);
}

#[test]
fn high_level_nddo_dispatch_uses_the_requested_method() {
    let molecule = Molecule::from_xyz_str(WATER, 0.0).unwrap();
    for method in [Method::Mndo, Method::MndoD] {
        let result = run_nddo(&molecule, method, &NddoOptions::default()).unwrap();
        assert!(result.converged, "{method}");
        assert!(result.total_ev.is_finite(), "{method}");
    }
}

#[test]
fn embedded_nddo_parameters_are_cached_per_method() {
    let first = NddoParameters::cached_for_method(Method::Mndo).unwrap();
    let second = NddoParameters::cached_for_method(Method::Mndo).unwrap();
    assert!(std::ptr::eq(first, second));
    assert_eq!(first.method, Method::Mndo);

    let mndod = NddoParameters::cached_for_method(Method::MndoD).unwrap();
    assert_eq!(mndod.method, Method::MndoD);
    assert!(!std::ptr::eq(first, mndod));
}

#[test]
fn embedded_zindo_parameters_are_cached() {
    let first = ZindoParameters::cached().unwrap();
    let second = ZindoParameters::cached().unwrap();
    assert!(std::ptr::eq(first, second));
    assert!(!first.supported_sp_elements().is_empty());
}

#[test]
fn method_capabilities_match_the_public_dispatchers() {
    assert_eq!(
        methods_for_api("single_point"),
        vec![
            Method::Cndo2,
            Method::Indo,
            Method::ZindoS,
            Method::Mindo3,
            Method::Mndo,
            Method::MndoD,
        ]
    );
    assert_eq!(
        methods_for_api("gradient"),
        vec![
            Method::Cndo2,
            Method::Indo,
            Method::ZindoS,
            Method::Mindo3,
            Method::Mndo,
            Method::MndoD,
        ]
    );
    assert_eq!(methods_for_api("hessian"), methods_for_api("gradient"));
    assert_eq!(methods_for_api("optimize"), methods_for_api("gradient"));
    assert_eq!(methods_for_api("excited_states"), vec![Method::ZindoS]);
    assert_eq!(methods_for_api("uv_vis_spectrum"), vec![Method::ZindoS]);
    assert_eq!(
        methods_for_api("excited_state_gradients"),
        vec![Method::ZindoS]
    );
    assert_eq!(
        methods_for_api("excited_state_hessians"),
        vec![Method::ZindoS]
    );
}

#[test]
fn non_nddo_method_is_rejected_without_substitution() {
    let molecule = Molecule::from_xyz_str(WATER, 0.0).unwrap();
    let error = run_nddo(&molecule, Method::Mindo3, &NddoOptions::default()).unwrap_err();
    assert!(matches!(error, XndoError::InvalidInput(_)));
    assert!(error.to_string().contains("not an MNDO/NDDO method"));
}

#[test]
fn pair_cache_limit_reports_a_recoverable_error() {
    let mut xyz = String::from("180\nwater lattice\n");
    for index in 0..60 {
        let x = 3.1 * (index % 5) as f64;
        let y = 3.1 * ((index / 5) % 4) as f64;
        let z = 3.1 * (index / 20) as f64;
        xyz.push_str(&format!(
            "O {x} {y} {z}\nH {} {y} {z}\nH {} {} {z}\n",
            x + 0.96,
            x - 0.24,
            y + 0.93
        ));
    }
    let molecule = Molecule::from_xyz_str(&xyz, 0.0).unwrap();
    let options = NddoOptions {
        integral_memory_mb: 1,
        ..NddoOptions::default()
    };
    let error = run_nddo(&molecule, Method::Mndo, &options).unwrap_err();
    assert!(matches!(error, XndoError::ResourceLimit { .. }));
}
