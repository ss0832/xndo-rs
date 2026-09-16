// SPDX-License-Identifier: GPL-3.0-or-later

//! The orbital-energy / HOMO-LUMO surface, checked against the engines rather
//! than against itself.
//!
//! `src/orbitals.rs` unit-tests the arithmetic on made-up spectra. What it
//! cannot check is that each of the four engines is wired to it correctly, that
//! all six methods answer, and that the frontier pair is the one the converged
//! SCF actually has. That is what these do.

use xndo_rs::{
    run_method, run_orbitals, CalculationResult, Method, Molecule, NddoOptions, Reference,
};

const WATER: &str =
    "3\nwater\nO 0.0000 0.0000 0.0000\nH 0.9584 0.0000 0.0000\nH -0.2400 0.9278 0.0000\n";

/// The methyl radical: an open-shell doublet, and the case where reading the
/// alpha channel alone gives the wrong LUMO.
const METHYL: &str = "4\nmethyl\nC 0.0000 0.0000 0.0000\nH 1.0790 0.0000 0.0000\n\
    H -0.5395 0.9344 0.0000\nH -0.5395 -0.9344 0.0000\n";

/// Every method with a native ground-state engine. CNDO/2 and INDO only carry
/// H, Li, C, N, O (and S for CNDO/2), which water and methyl are inside.
const ALL: [Method; 6] = [
    Method::Cndo2,
    Method::Indo,
    Method::Mndo,
    Method::MndoD,
    Method::Mindo3,
    Method::ZindoS,
];

fn options(multiplicity: usize) -> NddoOptions {
    NddoOptions {
        multiplicity,
        reference: Reference::Auto,
        ..NddoOptions::default()
    }
}

#[test]
fn every_method_reports_a_frontier_pair_and_a_full_spectrum() {
    let molecule = Molecule::from_xyz_str(WATER, 0.0).unwrap();
    for method in ALL {
        let o = run_orbitals(&molecule, method, &options(1)).unwrap_or_else(|e| {
            panic!("{method}: {e}");
        });
        // Water in a minimal valence basis: 4 AOs on O, 1 on each H.
        assert_eq!(o.n_orbitals(), 6, "{method}");
        assert_eq!(o.n_alpha, 4, "{method}");
        assert_eq!(o.n_beta, 4, "{method}");
        assert!(!o.unrestricted_reference(), "{method}");

        let homo = o.homo_ev().unwrap_or_else(|| panic!("{method}: no HOMO"));
        let lumo = o.lumo_ev().unwrap_or_else(|| panic!("{method}: no LUMO"));
        assert!(homo.is_finite() && lumo.is_finite(), "{method}");
        // A closed-shell ground state at its aufbau solution.
        assert!(
            lumo > homo,
            "{method}: LUMO {lumo} is not above HOMO {homo}"
        );
        assert_eq!(o.gap_ev(), Some(lumo - homo), "{method}");
        // The frontier pair must be the spectrum's, not a separate calculation.
        assert_eq!(homo, o.alpha_ev[o.n_alpha - 1], "{method}");
        assert_eq!(lumo, o.alpha_ev[o.n_alpha], "{method}");
        // Ascending, as the eigensolver returns it: a descending spectrum would
        // make every index above mean the opposite orbital.
        assert!(
            o.alpha_ev.windows(2).all(|w| w[0] <= w[1]),
            "{method}: spectrum is not ascending: {:?}",
            o.alpha_ev
        );
        assert_eq!(
            o.occupations_alpha(),
            vec![2.0, 2.0, 2.0, 2.0, 0.0, 0.0],
            "{method}"
        );
        assert_eq!(o.occupations_beta(), None, "{method}");
    }
}

#[test]
fn an_open_shell_doublet_takes_its_lumo_from_the_beta_channel() {
    let molecule = Molecule::from_xyz_str(METHYL, 0.0)
        .unwrap()
        .with_multiplicity(2);
    for method in ALL {
        let o = run_orbitals(&molecule, method, &options(2)).unwrap_or_else(|e| {
            panic!("{method}: {e}");
        });
        assert!(o.unrestricted_reference(), "{method}");
        assert_eq!((o.n_alpha, o.n_beta), (4, 3), "{method}");

        let lumo_alpha = o.lumo_alpha_ev().unwrap();
        let lumo_beta = o.lumo_beta_ev().unwrap();
        // The beta partner of the singly occupied orbital. This is the whole
        // reason the combined reading exists: before v0.3.0 the NDDO engine
        // discarded the beta eigenvalues and reported `lumo_alpha` as the LUMO.
        assert!(
            lumo_beta < lumo_alpha,
            "{method}: beta LUMO {lumo_beta} is not below alpha LUMO {lumo_alpha}"
        );
        assert_eq!(o.lumo_ev(), Some(lumo_beta), "{method}");
        assert_eq!(o.homo_ev(), o.homo_alpha_ev(), "{method}");
        assert_eq!(
            o.occupations_alpha(),
            vec![1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
            "{method}"
        );
        assert_eq!(
            o.occupations_beta(),
            Some(vec![1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
            "{method}"
        );
    }
}

#[test]
fn planar_methyl_puts_its_unpaired_electron_on_the_out_of_plane_p() {
    // The regression behind `tests/data/ORACLE_NOTES.md` item 25. MINDO/3 used
    // to converge to a real, aufbau-satisfying SCF solution 3.07 eV above this
    // one, with the odd electron in the C-H antibonding a1* and a p_z spin
    // population of exactly zero.
    //
    // Reference values are MOPAC7 1.15's own, for this geometry, under
    // `MINDO3 1SCF PRECISE UHF DOUBLET`: total energy -169.36017 eV, ionisation
    // potential 9.65591 eV, spin populations (s, px, py, pz) = (0.07479,
    // 0.06865, 0.06865, 1.00000). MOPAC7 is the program the MINDO/3 tables come
    // from, so this is the source and not a proxy.
    //
    // The assertion that matters is the p_z one: energy alone would also pass
    // for a solution that happened to be close, and it is *where the electron
    // is* that was wrong.
    let molecule = Molecule::from_xyz_str(METHYL, 0.0)
        .unwrap()
        .with_multiplicity(2);
    let result = run_method(&molecule, Method::Mindo3, &options(2)).unwrap();
    let CalculationResult::Mindo3(r) = &result else {
        panic!("MINDO/3 did not run through its own engine");
    };

    let total = r.total_ev;
    assert!(
        (total - -169.36017).abs() < 2.0e-3,
        "total energy {total} eV is not MOPAC7's -169.36017 (the wrong solution is -166.294)"
    );
    let homo = r.homo_ev().unwrap();
    assert!(
        (-homo - 9.65591).abs() < 2.0e-3,
        "ionisation potential {} eV is not MOPAC7's 9.65591",
        -homo
    );

    // Carbon is atom 0, AO order s, px, py, pz, and the molecule is in the xy
    // plane, so p_z is AO 3 and carries the unpaired electron.
    let spin = r
        .spin_density
        .as_ref()
        .expect("an open-shell doublet has a spin density");
    for (ao, want, label) in [
        (0usize, 0.07479, "C s"),
        (1, 0.06865, "C px"),
        (2, 0.06865, "C py"),
        (3, 1.00000, "C pz"),
    ] {
        let got = spin[(ao, ao)];
        assert!(
            (got - want).abs() < 5.0e-4,
            "{label} spin population {got} is not MOPAC7's {want}"
        );
    }
}

#[test]
fn the_nddo_result_fields_agree_with_the_dispatched_accessor() {
    // `NddoResult` carries `homo_ev`/`lumo_ev` as fields and everything else
    // reaches them through `orbitals()`. Two ways to say one thing is two ways
    // to drift, so pin them together.
    let molecule = Molecule::from_xyz_str(METHYL, 0.0)
        .unwrap()
        .with_multiplicity(2);
    for method in [Method::Mndo, Method::MndoD] {
        let result = run_method(&molecule, method, &options(2)).unwrap();
        let CalculationResult::Nddo(r) = &result else {
            panic!("{method} did not run through the NDDO engine");
        };
        let o = r.orbitals();
        assert_eq!(r.homo_ev, o.homo_ev(), "{method}");
        assert_eq!(r.lumo_ev, o.lumo_ev(), "{method}");
        assert_eq!(r.homo_lumo_gap_ev(), o.gap_ev(), "{method}");
        assert_eq!(r.n_alpha, r.n_occ, "{method}");
        // The beta channel is kept, not dropped.
        let beta = r
            .mo_energies_beta
            .as_ref()
            .unwrap_or_else(|| panic!("{method}: UHF beta eigenvalues were discarded"));
        assert_eq!(beta.len(), r.mo_energies.len(), "{method}");
        let cb = r
            .mo_coeff_beta
            .as_ref()
            .unwrap_or_else(|| panic!("{method}: UHF beta coefficients were discarded"));
        assert_eq!((cb.rows, cb.cols), (r.mo_coeff.rows, r.mo_coeff.cols));
        // The beta coefficients must belong to the beta eigenvalues, not be a
        // copy of the alpha ones: an open shell breaks spin symmetry.
        assert_ne!(beta, &r.mo_energies, "{method}");
    }
}

#[test]
fn a_restricted_reference_leaves_the_beta_channel_absent_rather_than_duplicated() {
    // `None` and "a copy of alpha" are not the same claim: the first says there
    // is one spectrum, the second says there are two that happen to agree.
    let molecule = Molecule::from_xyz_str(WATER, 0.0).unwrap();
    for method in ALL {
        let o = run_orbitals(&molecule, method, &options(1)).unwrap();
        assert_eq!(o.beta_ev, None, "{method}");
        // ... but a caller asking for the beta channel still gets the right
        // numbers, because for a restricted reference they are the alpha ones.
        assert_eq!(o.beta_or_alpha_ev(), o.alpha_ev.as_slice(), "{method}");
        assert_eq!(o.homo_beta_ev(), o.homo_alpha_ev(), "{method}");
    }
}

#[test]
fn the_frontier_pair_survives_an_explicit_uhf_request_on_a_closed_shell() {
    // Forcing UHF on a singlet should reach the same fixed point, so the
    // frontier pair must not move. It also exercises the path where both
    // channels are present and equal.
    //
    // Both SCFs are tightened well past the defaults first. At the default
    // thresholds the two paths stop at slightly different residuals and the
    // orbital energies differ in the seventh decimal -- which says nothing about
    // whether they found the same solution. Tightening is what makes the
    // comparison mean something: the agreement below is then the fixed point's,
    // not the convergence threshold's.
    let molecule = Molecule::from_xyz_str(WATER, 0.0).unwrap();
    let tight = NddoOptions {
        e_tol: 1.0e-13,
        p_tol: 1.0e-11,
        max_scf: 400,
        ..NddoOptions::default()
    };
    let uhf = NddoOptions {
        reference: Reference::Uhf,
        ..tight.clone()
    };
    for method in ALL {
        let rhf = run_orbitals(&molecule, method, &tight).unwrap();
        let unr = run_orbitals(&molecule, method, &uhf).unwrap();
        assert!(unr.unrestricted_reference(), "{method}");
        let (a, b) = (rhf.homo_ev().unwrap(), unr.homo_ev().unwrap());
        assert!(
            (a - b).abs() < 1.0e-8,
            "{method}: HOMO moved from {a} to {b} between RHF and UHF on a closed shell"
        );
        let (a, b) = (rhf.lumo_ev().unwrap(), unr.lumo_ev().unwrap());
        assert!(
            (a - b).abs() < 1.0e-8,
            "{method}: LUMO moved from {a} to {b} between RHF and UHF on a closed shell"
        );
        // Occupations change meaning with the reference: 2.0 per spatial
        // orbital restricted, 1.0 per spin orbital unrestricted.
        assert_eq!(unr.occupations_alpha()[0], 1.0, "{method}");
        assert_eq!(rhf.occupations_alpha()[0], 2.0, "{method}");
    }
}
