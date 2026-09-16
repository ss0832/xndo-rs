// SPDX-License-Identifier: GPL-3.0-or-later

//! What ZINDO/S-CIS does to photoexcitation, measured rather than asserted.
//!
//! This is the study, run as a test. Every number in
//! `studies/photoexcitation/results.tsv` and in the report rendered from it is
//! produced here; nothing is transcribed by hand.
//!
//! The measurements are **self-contained**: each one compares the model against
//! something known exactly -- an asymptotic form, an exact degeneracy, a
//! symmetry -- rather than against a literature value. That is a deliberate
//! limit. With post-HF out of scope for v0.3.0 there is no way, inside this
//! crate, to separate the CIS ansatz's error from the semiempirical
//! Hamiltonian's, so this study does not claim to; what it can do exactly, it
//! does exactly.
//!
//! Set `XNDO_STUDY_WRITE=1` to regenerate `studies/photoexcitation/results.tsv`.

use std::fmt::Write as _;

use xndo_rs::zindo::{ZINDO_AU2ANG, ZINDO_AU2EV, ZINDO_TOMK};
use xndo_rs::{
    zindo_s_cis_spin, CisSpin, ExcitedState, Molecule, Reference, ZindoOptions, ZindoParameters,
};

/// `e^2 / 4 pi eps0` in eV*Angstrom: the exact Coulomb interaction of two unit
/// charges one Angstrom apart, and the coefficient a charge-transfer state's
/// `-C/R` asymptote must have.
///
/// Written from the model's own Hartree and Bohr so that this is the same
/// number ZINDO/S would produce with `TOMK = 1`, and the comparison below is
/// about the factor and nothing else.
const EXACT_COULOMB_EV_ANGSTROM: f64 = ZINDO_AU2EV * ZINDO_AU2ANG;

/// Transition-dipole magnitude (atomic units) below which a state is called
/// dark. See [`out_of_plane_fraction`] for why it is here and not lower.
const DARK_DIPOLE_AU: f64 = 1.0e-5;

fn params() -> &'static ZindoParameters {
    ZindoParameters::cached().expect("bundled ZINDO/S parameters")
}

fn options(n_states: usize) -> ZindoOptions {
    ZindoOptions {
        reference: Reference::Rhf,
        n_states,
        ..ZindoOptions::default()
    }
}

fn roots(xyz: &str, spin: CisSpin, n_states: usize) -> Vec<ExcitedState> {
    let mol = Molecule::from_xyz_str(xyz, 0.0).expect("geometry");
    zindo_s_cis_spin(&mol, params(), &options(n_states), spin)
        .expect("CIS")
        .states
}

/// A planar ring of `symbols` on a circle of `radius`, each carrying a radial
/// substituent at `r_sub` unless its entry is `None`.
///
/// Idealised: a regular polygon, not an optimised geometry. That is the right
/// choice here and a wrong one elsewhere, so it is worth being explicit. This
/// study characterises what the *model* does with a state of a given character
/// -- whether the lowest singlet is n->pi*, how far a state's hole and electron
/// separate, how much a window truncation moves a root. Those depend on the
/// symmetry and the connectivity, which a regular polygon has exactly right.
/// They would be the wrong geometries for comparing an excitation energy with
/// an experimental band maximum, which this study does not do.
fn planar_ring(symbols: &[&str], radius: f64, substituents: &[Option<(&str, f64)>]) -> String {
    let n = symbols.len();
    let mut atoms = Vec::new();
    for (i, sym) in symbols.iter().enumerate() {
        let theta = std::f64::consts::TAU * i as f64 / n as f64;
        let (s, c) = theta.sin_cos();
        atoms.push(format!(
            "{sym} {:.6} {:.6} 0.000000",
            radius * c,
            radius * s
        ));
        if let Some(Some((sub, r_sub))) = substituents.get(i) {
            let rr = radius + r_sub;
            atoms.push(format!("{sub} {:.6} {:.6} 0.000000", rr * c, rr * s));
        }
    }
    format!("{}\nring\n{}\n", atoms.len(), atoms.join("\n"))
}

/// The chromophores this study runs over: `(name, xyz, expected character, why)`.
///
/// Composition is C/N/O/H only, so every atom is inside the ZINDO/S s-p branch
/// and nothing here is confounded by a missing parameter. The set is built
/// around the distinction that matters for the model -- states where the hole
/// is a pi orbital against states where it is a lone pair -- because that is
/// where ZINDO/S's documented failure lives.
struct Chromophore {
    name: &'static str,
    xyz: String,
    /// Unit normal of the molecular plane. Every geometry in this set is planar
    /// and is constructed in a coordinate plane, so the normal is known exactly
    /// rather than fitted -- which is what lets [`out_of_plane_fraction`] be an
    /// exact statement instead of an estimate.
    normal: [f64; 3],
    /// Whether the molecule has a lone pair that can donate into pi*, i.e.
    /// whether an n->pi* state exists at all. Read off the structure, not off
    /// the calculation.
    has_lone_pair: bool,
    why: &'static str,
}

/// Out-of-plane axis for the two coordinate planes the geometries use.
const NORMAL_Y: [f64; 3] = [0.0, 1.0, 0.0];
const NORMAL_Z: [f64; 3] = [0.0, 0.0, 1.0];

fn chromophores() -> Vec<Chromophore> {
    let c = |name, xyz, normal, has_lone_pair, why| Chromophore {
        name,
        xyz,
        normal,
        has_lone_pair,
        why,
    };
    vec![
        c(
            "ethene",
            "6\nethene\nC 0.000000 0.000000 0.667000\nC 0.000000 0.000000 -0.667000\n\
             H 0.929000 0.000000 1.237000\nH -0.929000 0.000000 1.237000\n\
             H 0.929000 0.000000 -1.237000\nH -0.929000 0.000000 -1.237000\n"
                .into(),
            NORMAL_Y,
            false,
            "the smallest pi system; one double bond and nothing else",
        ),
        c(
            "butadiene",
            "10\ns-trans butadiene\nC -1.828000 -0.396000 0.000000\nC -0.633000 0.192000 0.000000\n\
             C 0.633000 -0.192000 0.000000\nC 1.828000 0.396000 0.000000\n\
             H -2.731000 0.209000 0.000000\nH -1.918000 -1.478000 0.000000\n\
             H -0.560000 1.277000 0.000000\nH 0.560000 -1.277000 0.000000\n\
             H 1.918000 1.478000 0.000000\nH 2.731000 -0.209000 0.000000\n"
                .into(),
            NORMAL_Z,
            false,
            "conjugation over two double bonds; the first step of the polyene series",
        ),
        c(
            "benzene",
            planar_ring(&["C"; 6], 1.397, &[Some(("H", 1.084)); 6]),
            NORMAL_Z,
            false,
            "aromatic, fully delocalised, degenerate frontier orbitals",
        ),
        c(
            "pyridine",
            {
                let mut subs = vec![Some(("H", 1.084)); 6];
                subs[0] = None;
                planar_ring(&["N", "C", "C", "C", "C", "C"], 1.397, &subs)
            },
            NORMAL_Z,
            true,
            "one ring nitrogen: an in-plane lone pair above the pi system",
        ),
        c(
            "pyrazine",
            {
                let mut subs = vec![Some(("H", 1.084)); 6];
                subs[0] = None;
                subs[3] = None;
                planar_ring(&["N", "C", "C", "N", "C", "C"], 1.397, &subs)
            },
            NORMAL_Z,
            true,
            "two para nitrogens: two lone pairs, and the classic n->pi* test case",
        ),
        c(
            "formaldehyde",
            "4\nformaldehyde\nC 0.000000 0.000000 0.000000\nO 0.000000 0.000000 1.203000\n\
             H 0.943000 0.000000 -0.545000\nH -0.943000 0.000000 -0.545000\n"
                .into(),
            NORMAL_Y,
            true,
            "the reference carbonyl; the oxygen lone pair is the highest occupied orbital",
        ),
        c(
            "glyoxal",
            "6\ns-trans glyoxal\nC -0.752000 0.000000 0.000000\nC 0.752000 0.000000 0.000000\n\
             O -1.383000 1.043000 0.000000\nO 1.383000 -1.043000 0.000000\n\
             H -1.283000 -0.958000 0.000000\nH 1.283000 0.958000 0.000000\n"
                .into(),
            NORMAL_Z,
            true,
            "two conjugated carbonyls; lone pairs that can mix with each other",
        ),
        c(
            "acrolein",
            "8\ns-trans acrolein\nC 1.219000 -0.400000 0.000000\nC 0.000000 0.243000 0.000000\n\
             C -1.219000 -0.400000 0.000000\nO -2.320000 0.120000 0.000000\n\
             H 2.140000 0.170000 0.000000\nH 1.280000 -1.483000 0.000000\n\
             H 0.000000 1.328000 0.000000\nH -1.180000 -1.500000 0.000000\n"
                .into(),
            NORMAL_Z,
            true,
            "a carbonyl conjugated with a double bond: n->pi* and pi->pi* close together",
        ),
        c(
            "formamide",
            "6\nformamide\nC 0.000000 0.000000 0.000000\nO 0.000000 0.000000 1.219000\n\
             N 1.194000 0.000000 -0.683000\nH -0.939000 0.000000 -0.560000\n\
             H 1.206000 0.000000 -1.690000\nH 2.056000 0.000000 -0.180000\n"
                .into(),
            NORMAL_Y,
            true,
            "the amide chromophore; a nitrogen lone pair conjugated into the carbonyl",
        ),
    ]
}

/// A donor and an acceptor, rigid, their centres `r` Angstrom apart along z.
///
/// Ammonia donates and formaldehyde accepts: the nitrogen lone pair is the
/// highest occupied orbital in the pair and the carbonyl pi* the lowest empty
/// one, so the lowest charge-transfer state is unambiguously N -> C=O. They are
/// held far enough apart that no bond forms, and only the separation changes.
fn donor_acceptor(r: f64) -> String {
    let mut s = String::from("8\nammonia + formaldehyde\n");
    // Ammonia, pyramidal, nitrogen at the origin, hydrogens pointing -z.
    let _ = writeln!(s, "N 0.000000 0.000000 0.000000");
    for (x, y) in [(0.9377, 0.0), (-0.4689, 0.8121), (-0.4689, -0.8121)] {
        let _ = writeln!(s, "H {x:.6} {y:.6} -0.3816");
    }
    // Formaldehyde, planar, carbon on the axis at +r, oxygen pointing away.
    let _ = writeln!(s, "C 0.000000 0.000000 {:.6}", r);
    let _ = writeln!(s, "O 0.000000 0.000000 {:.6}", r + 1.203);
    let _ = writeln!(s, "H 0.943000 0.000000 {:.6}", r - 0.545);
    let _ = writeln!(s, "H -0.943000 0.000000 {:.6}", r - 0.545);
    s
}

/// Least-squares fit of `y = c0 + c1*x + c2*x^2`, returning the coefficients.
///
/// Quadratic and not linear because of what the model's Coulomb integral
/// actually is. Mataga-Nishimoto is `C/(R + a)`, not `C/R`, and
///
/// ```text
///     C/(R + a) = C/R - C*a/R^2 + O(1/R^3)
/// ```
///
/// so a fit in `1/R` alone absorbs the `a` term into the slope and reports a
/// coefficient that depends on the range it was fitted over. `a` is about
/// 1.4 Angstrom here, which is not a small fraction of any separation a
/// molecule survives, so this matters: the quadratic term is what lets the
/// `1/R` coefficient be read off without the screening length contaminating it.
///
/// Solved by Cramer's rule on the 3x3 normal equations. The design matrix is a
/// Vandermonde in `1/R` over a range of two, which is well enough conditioned
/// at this size that the textbook solution is the right one.
fn quadratic_fit(points: &[(f64, f64)]) -> [f64; 3] {
    let power_sum = |k: u32| -> f64 { points.iter().map(|p| p.0.powi(k as i32)).sum() };
    let moment = |k: u32| -> f64 { points.iter().map(|p| p.1 * p.0.powi(k as i32)).sum() };
    let a = [
        [power_sum(0), power_sum(1), power_sum(2)],
        [power_sum(1), power_sum(2), power_sum(3)],
        [power_sum(2), power_sum(3), power_sum(4)],
    ];
    let b = [moment(0), moment(1), moment(2)];
    let det3 = |m: &[[f64; 3]; 3]| -> f64 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det3(&a);
    let mut out = [0.0; 3];
    for (k, slot) in out.iter_mut().enumerate() {
        let mut m = a;
        for row in 0..3 {
            m[row][k] = b[row];
        }
        *slot = det3(&m) / d;
    }
    out
}

#[test]
fn the_charge_transfer_asymptote_keeps_its_shape_and_loses_its_coefficient() {
    // The sharpest thing this model can be held to, because the answer is known
    // exactly and is not a matter of parameterisation: two separated unit
    // charges attract as `-e^2/(4 pi eps0 R)`, so a charge-transfer state's
    // energy must approach its asymptote as `-14.3996 eV*A / R`.
    //
    // ZINDO/S will not. Its two-centre Coulomb integral is Mataga-Nishimoto
    // with `TOMK = 1.2` (`zindo.rs`):
    //
    //     gamma(R) = TOMK * AU2EV / (R_bohr + TOMK * 2 * AU2EV / (g_A + g_B))
    //
    // whose large-R limit is `TOMK * AU2EV / R_bohr`, i.e. **1.2 times** the
    // exact Coulomb interaction. So the prediction is specific: the CT state
    // keeps the correct `1/R` *form* -- unlike TD-DFT, whose most documented
    // charge-transfer failure is losing the `1/R` term altogether -- and gets
    // its coefficient wrong by exactly the parameterisation factor.
    //
    // Stated before running it, which is the only way it means anything: the
    // fitted slope should be -17.28 eV*A rather than -14.40, and the ratio of
    // the two should come back as TOMK itself.
    // Far enough out that the next term in the expansion is negligible. The
    // Mataga-Nishimoto screening length here is about 1.7 Angstrom, so the
    // `1/R^3` term the quadratic fit does not carry is worth `C*a^2/R^3`:
    // 1.5e-2 eV at 15 Angstrom, which is 1.6% of the coefficient, and under
    // 1e-4 eV beyond 60. Fitting closer in and calling the answer 1.18 would
    // have been measuring the fitting range.
    let separations = [60.0_f64, 80.0, 100.0, 140.0, 200.0, 280.0, 400.0];
    let mut points = Vec::new();
    for r in separations {
        let states = roots(&donor_acceptor(r), CisSpin::Singlet, 20);
        // The charge-transfer root identifies itself: the engine reports the
        // distance between the hole and electron centroids, and for a state
        // that moves an electron from one fragment to the other that distance
        // is the fragment separation. A local excitation leaves it near zero.
        // Picking the state by the descriptor rather than by root index is what
        // keeps this honest as the ordering changes with R.
        let ct = states
            .iter()
            .filter(|s| s.charge_transfer_distance_angstrom > 0.5 * r)
            .min_by(|a, b| a.energy_ev.partial_cmp(&b.energy_ev).unwrap())
            .unwrap_or_else(|| {
                panic!(
                    "no charge-transfer root at R = {r} A; largest hole-electron \
                     separation found was {:.2} A over {} roots",
                    states
                        .iter()
                        .map(|s| s.charge_transfer_distance_angstrom)
                        .fold(0.0_f64, f64::max),
                    states.len()
                )
            });
        // The abscissa is the *measured* hole-electron separation, not the
        // nominal fragment separation. They are not the same number: the hole
        // is the nitrogen lone pair and the electron is the carbonyl pi*, whose
        // centroid sits between the C and the O, roughly 0.6 A further away
        // again. Fitting against `r` charges that offset to the 1/R
        // coefficient and recovers 1.137 instead of 1.2 -- a 5% error in the
        // headline number, from using a separation the excitation does not
        // have. The engine already computes the right one.
        println!(
            "  R = {r:5.1} A   D_CT = {:6.2} A   E = {:.6} eV",
            ct.charge_transfer_distance_angstrom, ct.energy_ev
        );
        points.push((1.0 / ct.charge_transfer_distance_angstrom, ct.energy_ev));
    }

    let [c0, slope, c2] = quadratic_fit(&points);
    let expected = -ZINDO_TOMK * EXACT_COULOMB_EV_ANGSTROM;
    println!(
        "CT asymptote: E(R) = {c0:.4} {slope:+.4}/R {c2:+.4}/R^2 eV, R in Angstrom\n  \
         1/R coefficient {slope:.4}, exact Coulomb {:.4}, ratio {:.4} (TOMK = {ZINDO_TOMK})",
        -EXACT_COULOMB_EV_ANGSTROM,
        slope / -EXACT_COULOMB_EV_ANGSTROM,
    );

    // The fit has to describe the data before its coefficients mean anything:
    // coefficients pulled through a curve the model does not have would be an
    // artefact of the fitting range.
    let worst_residual = points
        .iter()
        .map(|(x, y)| (y - (c0 + slope * x + c2 * x * x)).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        worst_residual < 2.0e-3,
        "the CT energy does not follow E_inf - C/R + D/R^2 (worst residual \
         {worst_residual:.2e} eV); coefficients fitted through the wrong form \
         would not mean what this test claims"
    );

    assert!(
        (slope - expected).abs() < 0.05,
        "CT asymptote 1/R coefficient {slope:.4} eV*A, expected {expected:.4} \
         (= -TOMK * {EXACT_COULOMB_EV_ANGSTROM:.4})"
    );

    // And the thing the study reports: the model's 1/R coefficient divided by
    // the exact one is the parameterisation factor, recovered from excitation
    // energies alone rather than read out of the source.
    let recovered = slope / -EXACT_COULOMB_EV_ANGSTROM;
    assert!(
        (recovered - ZINDO_TOMK).abs() < 0.005,
        "recovered coefficient {recovered:.4}, expected TOMK = {ZINDO_TOMK}"
    );
}

/// The fraction of a transition dipole that lies along the molecular normal.
///
/// Exact for a planar molecule, and the reason is worth stating because it is
/// what makes this a measurement rather than a label. Reflection in the
/// molecular plane is a symmetry of every molecule in this set. Under it a pi
/// orbital is odd and an in-plane orbital (a sigma bond, or a lone pair) is
/// even; of the dipole components, the two in-plane ones are even and the
/// normal one is odd. So a transition that **preserves** reflection parity --
/// pi->pi*, or sigma->sigma* -- has an even transition density and can only
/// couple to the in-plane components, while one that **changes** it --
/// n->pi*, but equally sigma->pi* and pi->sigma* -- can only couple to the
/// normal one. No transition can have both.
///
/// So this measures a parity change, which is not the same as measuring a lone
/// pair: sigma orbitals are even too, and ethene, which has no lone pair at all,
/// has an out-of-plane polarised root at 8.73 eV that is sigma <-> pi. Reading
/// "out-of-plane" as "n->pi*" would misread it, the same way reading "dark" as
/// "n->pi*" would misread benzene's two lowest singlets, which are dark and
/// pi->pi*.
///
/// Returns `None` for a state whose transition dipole vanishes, where the
/// question has no answer: a dark state is dark, and a ratio of two zeros would
/// invent a polarisation it does not have.
///
/// `DARK_DIPOLE_AU` is where "vanishes" is drawn, and it has to be drawn with
/// care because a symmetry-forbidden transition dipole is zero in exact
/// arithmetic and merely tiny in floating point. At 1e-7 it fell *inside* the
/// noise: benzene's degenerate pair at 7.0717 eV came out one "dark" and one
/// "out-of-plane", which is not a distinction the molecule has. A threshold
/// that splits a degenerate pair is measuring the arithmetic. 1e-5 sits four
/// orders below the weakest genuinely allowed transition in this set
/// (formaldehyde's, f = 0.016, |mu| = 0.1 au) and above that noise, so it
/// separates the two questions instead of cutting through one of them.
fn out_of_plane_fraction(state: &ExcitedState, normal: [f64; 3]) -> Option<f64> {
    let mu = state.transition_dipole_au;
    let total = (mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2]).sqrt();
    if total < DARK_DIPOLE_AU {
        return None;
    }
    let along = mu[0] * normal[0] + mu[1] * normal[1] + mu[2] * normal[2];
    Some((along / total).abs())
}

/// One row of `studies/photoexcitation/results.tsv`.
struct Row {
    molecule: &'static str,
    spin: &'static str,
    root: usize,
    energy_ev: f64,
    oscillator: f64,
    /// `in-plane`, `out-of-plane`, or `dark` -- what the transition dipole's
    /// orientation says, which for a planar molecule is exact. See
    /// [`out_of_plane_fraction`].
    polarisation: String,
    d_ct_angstrom: f64,
    difference_dipole_debye: f64,
    /// Why this molecule is in the set. The oracle reference files carry the
    /// same column for the same reason: it is the record of what the row is
    /// for, and the thing a reviewer checks when deciding whether the set still
    /// covers what it claims after an edit.
    why: &'static str,
}

/// Roots per molecule per spin sector. Six is enough to reach past the
/// frontier pair into states where the window edge starts to matter, and few
/// enough that the table stays readable.
const SURVEY_ROOTS: usize = 6;

fn survey() -> Vec<Row> {
    let mut rows = Vec::new();
    for chromophore in chromophores() {
        for (spin, label) in [(CisSpin::Singlet, "singlet"), (CisSpin::Triplet, "triplet")] {
            for (index, state) in roots(&chromophore.xyz, spin, SURVEY_ROOTS)
                .iter()
                .enumerate()
            {
                let polarisation = match out_of_plane_fraction(state, chromophore.normal) {
                    None => "dark".to_string(),
                    Some(f) if f > 0.5 => "out-of-plane".to_string(),
                    Some(_) => "in-plane".to_string(),
                };
                rows.push(Row {
                    molecule: chromophore.name,
                    spin: label,
                    root: index + 1,
                    energy_ev: state.energy_ev,
                    oscillator: state.oscillator_strength,
                    polarisation,
                    d_ct_angstrom: state.charge_transfer_distance_angstrom,
                    difference_dipole_debye: state.difference_dipole_magnitude_debye,
                    why: chromophore.why,
                });
            }
        }
    }
    rows
}

#[test]
fn every_transition_dipole_is_purely_in_plane_or_purely_along_the_normal() {
    // Exact, for the reason [`out_of_plane_fraction`] gives: reflection in the
    // molecular plane is a symmetry of every molecule in this set, and it
    // separates the in-plane dipole components from the normal one. A
    // transition dipole with, say, 30% of its length along the normal would
    // mean the CI vector is mixing states of different reflection symmetry,
    // which no exact treatment of this Hamiltonian can do.
    //
    // This is the floor the polarisation column of the study stands on: without
    // it, calling a state "in-plane" would be a description of a continuum
    // rather than a statement about its symmetry.
    for chromophore in chromophores() {
        for state in roots(&chromophore.xyz, CisSpin::Singlet, SURVEY_ROOTS) {
            let Some(fraction) = out_of_plane_fraction(&state, chromophore.normal) else {
                continue;
            };
            assert!(
                !(1.0e-6..=1.0 - 1.0e-6).contains(&fraction),
                "{}: a root at {:.4} eV is polarised {:.1}% along the normal, \
                 which is neither in-plane nor out-of-plane; reflection in the \
                 molecular plane forbids the mixture",
                chromophore.name,
                state.energy_ev,
                100.0 * fraction
            );
        }
    }
}

#[test]
fn the_lowest_parity_changing_root_is_lower_when_there_is_a_lone_pair_to_supply_it() {
    // What the polarisation column can and cannot settle.
    //
    // It cannot identify an n->pi* state on its own: sigma orbitals are even
    // under the plane reflection just as lone pairs are, so sigma->pi* is
    // out-of-plane polarised too. Ethene has no lone pair and still has such a
    // root, at 8.73 eV. A first draft of this test asserted the opposite and
    // ethene falsified it, which is the useful thing a test can do.
    //
    // What it can settle is the ordering, and that *is* the chemistry: a lone
    // pair is a far better donor than a sigma bond, so where one exists the
    // lowest parity-changing excitation should come from it and should lie well
    // below where a molecule without one has to find its first sigma->pi*. The
    // separation between the two groups is the measurement.
    //
    // One caveat, and it is the reason formaldehyde reads 8.8 eV below rather
    // than the 3.2 eV of its n->pi*: reflection is not the only symmetry in
    // play. Formaldehyde is C2v, where n->pi* is A2, and A2 spans no component
    // of the dipole at all -- so that state is not out-of-plane polarised, it
    // is *completely* dark, and a measure built on the transition dipole's
    // direction cannot see it. It is in the table, at 3.187 eV, identifiable by
    // its 2.55 D difference dipole. That a state can be invisible to the very
    // observable used to classify it is the honest limit of this approach, and
    // the reason the results table carries the difference dipole beside the
    // polarisation rather than instead of it.
    let mut with = f64::INFINITY;
    let mut without = f64::NEG_INFINITY;
    for chromophore in chromophores() {
        let lowest = roots(&chromophore.xyz, CisSpin::Singlet, SURVEY_ROOTS)
            .iter()
            .filter(|s| out_of_plane_fraction(s, chromophore.normal).is_some_and(|f| f > 0.5))
            .map(|s| s.energy_ev)
            .fold(f64::INFINITY, f64::min);
        if !lowest.is_finite() {
            continue;
        }
        println!(
            "  {:<14} lowest parity-changing root {lowest:7.3} eV   (lone pair: {})",
            chromophore.name, chromophore.has_lone_pair
        );
        if chromophore.has_lone_pair {
            with = with.min(lowest);
        } else {
            without = without.max(lowest);
        }
    }
    assert!(
        with.is_finite() && without.is_finite(),
        "both groups must contribute a parity-changing root for this to compare anything"
    );
    assert!(
        with < without,
        "the lowest lone-pair-donated excitation ({with:.3} eV) is not below the \
         highest sigma-donated one ({without:.3} eV); the two groups do not separate"
    );
}

#[test]
fn the_lowest_triplet_lies_below_the_lowest_singlet_everywhere() {
    // Hund's rule for excited states, and an exact statement about CIS rather
    // than an empirical trend: for the same spatial excitation the triplet
    // energy is the singlet's minus twice the exchange integral, and that
    // integral is positive definite. A molecule where this fails has not found
    // the same spatial excitation in both sectors, which is worth knowing.
    for chromophore in chromophores() {
        let singlet = roots(&chromophore.xyz, CisSpin::Singlet, 3);
        let triplet = roots(&chromophore.xyz, CisSpin::Triplet, 3);
        let (s1, t1) = (singlet[0].energy_ev, triplet[0].energy_ev);
        assert!(
            t1 < s1,
            "{}: lowest triplet {t1:.4} eV is not below the lowest singlet {s1:.4}",
            chromophore.name
        );
    }
}

#[test]
fn an_excitation_energy_does_not_depend_on_where_the_molecule_is() {
    // Translation and rotation are exact invariances, so this is a correctness
    // floor the rest of the study stands on: a number that moved when the
    // molecule did would make every comparison below meaningless. The
    // transition dipole's *magnitude* is invariant too, while its components
    // rotate, so the magnitude is what is compared.
    let mol = "4\nformaldehyde\nC 0.000000 0.000000 0.000000\nO 0.000000 0.000000 1.203000\n\
               H 0.943000 0.000000 -0.545000\nH -0.943000 0.000000 -0.545000\n";
    let moved = "4\nformaldehyde, translated and rotated\n\
                 C 5.000000 -3.000000 2.000000\nO 5.000000 -1.797000 2.000000\n\
                 H 5.000000 -3.545000 2.943000\nH 5.000000 -3.545000 1.057000\n";
    let here = roots(mol, CisSpin::Singlet, 5);
    let there = roots(moved, CisSpin::Singlet, 5);
    for (a, b) in here.iter().zip(&there) {
        assert!(
            (a.energy_ev - b.energy_ev).abs() < 1.0e-9,
            "excitation energy moved with the molecule: {:.10} vs {:.10} eV",
            a.energy_ev,
            b.energy_ev
        );
        assert!(
            (a.oscillator_strength - b.oscillator_strength).abs() < 1.0e-9,
            "oscillator strength moved with the molecule: {:.10} vs {:.10}",
            a.oscillator_strength,
            b.oscillator_strength
        );
    }
}

/// How much a truncated CI window moves the roots it keeps.
///
/// `ZindoOptions::active_occupied` / `active_virtual` change every excitation
/// energy, so a study that used the full space for its small molecules and a
/// window for its large ones would be comparing two different methods and
/// reporting the difference as chemistry. The window is fixed across the whole
/// set for that reason; this measures what fixing it costs.
fn windowed_roots(xyz: &str, occupied: usize, virtual_: usize) -> Vec<f64> {
    let mol = Molecule::from_xyz_str(xyz, 0.0).expect("geometry");
    let options = ZindoOptions {
        reference: Reference::Rhf,
        n_states: 3,
        active_occupied: Some(occupied),
        active_virtual: Some(virtual_),
        ..ZindoOptions::default()
    };
    zindo_s_cis_spin(&mol, params(), &options, CisSpin::Singlet)
        .expect("CIS")
        .states
        .iter()
        .map(|s| s.energy_ev)
        .collect()
}

#[test]
fn a_truncated_ci_window_only_ever_raises_a_root() {
    // Variational, and therefore exact rather than a trend: CIS in a window is
    // CIS in the full space restricted to a subspace, so every root it returns
    // is a Rayleigh quotient over fewer configurations and can only lie at or
    // above the corresponding full-space root. A window that *lowered* a root
    // would mean the truncation was not a subspace of the full problem.
    //
    // The size of the rise is the number the study reports: it is the
    // methodological error a fixed window imposes on the whole set, and it is
    // not small.
    let mut worst = 0.0_f64;
    let mut worst_label = String::new();
    for chromophore in chromophores() {
        let full = roots(&chromophore.xyz, CisSpin::Singlet, 3);
        let windowed = windowed_roots(&chromophore.xyz, 3, 3);
        for (index, (f, w)) in full.iter().zip(&windowed).enumerate() {
            let shift = w - f.energy_ev;
            assert!(
                shift > -1.0e-9,
                "{}: root {} fell from {:.6} to {:.6} eV when the window was \
                 narrowed, but a subspace cannot lower a variational root",
                chromophore.name,
                index + 1,
                f.energy_ev,
                w
            );
            if shift > worst {
                worst = shift;
                worst_label = format!("{} root {}", chromophore.name, index + 1);
            }
        }
    }
    println!("widest 3x3-window shift: {worst:.4} eV at {worst_label}");
}

#[test]
fn the_results_table_is_written_and_is_what_the_report_is_rendered_from() {
    let rows = survey();
    assert!(
        rows.len() >= 80,
        "the survey collapsed to {} rows; the set is nine molecules times two \
         spin sectors times up to {SURVEY_ROOTS} roots",
        rows.len()
    );

    let mut out = String::new();
    out.push_str(
        "# ZINDO/S-CIS over a planar C/N/O/H chromophore set. Generated by\n\
         # tests/photoexcitation.rs; regenerate with XNDO_STUDY_WRITE=1.\n\
         # polarisation is exact for a planar molecule: reflection in the\n\
         # molecular plane separates the in-plane dipole components from the\n\
         # normal one, so no transition can be polarised partly along both.\n",
    );
    out.push_str("molecule\tspin\troot\tenergy_ev\toscillator\tpolarisation\td_ct_angstrom\tdiff_dipole_debye\twhy\n");
    for row in &rows {
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{:.4}\t{:.4}\t{}",
            row.molecule,
            row.spin,
            row.root,
            row.energy_ev,
            row.oscillator,
            row.polarisation,
            row.d_ct_angstrom,
            row.difference_dipole_debye,
            row.why,
        );
    }

    let path = std::path::Path::new("studies/photoexcitation/results.tsv");
    if std::env::var("XNDO_STUDY_WRITE").is_ok() {
        std::fs::create_dir_all(path.parent().unwrap()).expect("studies directory");
        std::fs::write(path, &out).expect("write results");
        return;
    }
    // Otherwise the committed file must be what this run produces. The study is
    // a regression suite as much as a report: if an optimisation moves one of
    // these numbers, that is the thing to notice.
    let committed = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!(
            "{}: {e}; regenerate with XNDO_STUDY_WRITE=1",
            path.display()
        )
    });
    assert_eq!(
        committed.replace("\r\n", "\n"),
        out,
        "studies/photoexcitation/results.tsv is stale; regenerate with \
         XNDO_STUDY_WRITE=1 and review the diff"
    );
}
