// SPDX-License-Identifier: GPL-3.0-or-later

//! Proof that a written Molden file describes the wavefunction the engine
//! computed.
//!
//! The load-bearing tests here read the file back with a parser written from
//! the Molden format alone -- it shares no code with the writer beyond the
//! Gaussian integrals themselves -- and rebuild the basis from the text. A test
//! that asked the writer what it had written would pass on a file no other
//! program could read.

use xndo_rs::gto::{shell_overlap, Shell};
use xndo_rs::math::Vec3;
use xndo_rs::molden::{molden_string, MoldenCoefficients};
use xndo_rs::sto::PrimitiveGaussian;
use xndo_rs::{constants::ANGSTROM_TO_BOHR, run_method, Method, Molecule, NddoOptions};

const WATER: &str =
    "3\nwater\nO 0.0000 0.0000 0.0000\nH 0.9584 0.0000 0.0000\nH -0.2400 0.9278 0.0000\n";
const FORMALDEHYDE: &str = "4\nformaldehyde\nC 0.000000 0.000000 0.000000\n\
     O 0.000000 0.000000 1.203000\nH 0.943000 0.000000 -0.545000\n\
     H -0.943000 0.000000 -0.545000\n";
/// Hydrogen sulfide: sulfur carries d functions under MNDO/d, so this is the
/// case that exercises the `[5D]` block and the d permutation.
const H2S: &str = "3\nhydrogen sulfide\nS 0.000000 0.000000 0.000000\n\
     H 1.334800 0.000000 0.000000\nH -0.048700 1.334000 0.000000\n";

/// A Molden file, parsed back from its text.
struct Parsed {
    /// Atomic numbers and positions in Bohr.
    atoms: Vec<(u8, Vec3)>,
    shells: Vec<Shell>,
    /// `[orbital][ao]`.
    coefficients: Vec<Vec<f64>>,
    occupations: Vec<f64>,
    spins: Vec<String>,
    spherical_d: bool,
}

/// Read a Molden file using only the format's own rules.
fn parse(text: &str) -> Parsed {
    let mut atoms: Vec<(u8, Vec3)> = Vec::new();
    let mut shells = Vec::new();
    let mut coefficients: Vec<Vec<f64>> = Vec::new();
    let mut occupations = Vec::new();
    let mut spins = Vec::new();
    let mut spherical_d = false;

    let mut section = String::new();
    let mut angstrom = false;
    let mut current_atom: Option<usize> = None;
    let mut pending: Option<(usize, usize)> = None; // (l, remaining primitives)
    let mut primitives: Vec<PrimitiveGaussian> = Vec::new();

    for raw in text.lines() {
        let line = raw.trim();
        if line.starts_with('[') {
            let tag = line.to_ascii_lowercase();
            if tag.starts_with("[5d") {
                spherical_d = true;
            }
            section = tag
                .split(']')
                .next()
                .unwrap_or("")
                .trim_start_matches('[')
                .to_string();
            angstrom = tag.contains("angs");
            continue;
        }
        if line.is_empty() {
            // A blank line closes an atom's [GTO] block.
            if section == "gto" {
                current_atom = None;
            }
            continue;
        }
        match section.as_str() {
            "atoms" => {
                let f: Vec<&str> = line.split_whitespace().collect();
                let z: u8 = f[2].parse().unwrap();
                let scale = if angstrom { ANGSTROM_TO_BOHR } else { 1.0 };
                atoms.push((
                    z,
                    Vec3::new(
                        f[3].parse::<f64>().unwrap() * scale,
                        f[4].parse::<f64>().unwrap() * scale,
                        f[5].parse::<f64>().unwrap() * scale,
                    ),
                ));
            }
            "gto" => {
                let f: Vec<&str> = line.split_whitespace().collect();
                if let Some((l, remaining)) = pending {
                    let exponent: f64 = f[0].replace('D', "E").parse().unwrap();
                    let coefficient: f64 = f[1].replace('D', "E").parse().unwrap();
                    // The file gives the contraction coefficient; the reader
                    // supplies the Cartesian normalisation. Doing that here,
                    // from the format's definition, is what makes this an
                    // independent check of the writer's `primitive_norm`.
                    let norm = (2.0 / std::f64::consts::PI * exponent).powf(0.75)
                        * (4.0 * exponent).sqrt().powi(l as i32)
                        / if l == 2 { 3.0_f64.sqrt() } else { 1.0 };
                    primitives.push(PrimitiveGaussian {
                        exponent,
                        coefficient: coefficient * norm,
                    });
                    if remaining == 1 {
                        shells.push(Shell {
                            centre: atoms[current_atom.unwrap()].1,
                            l,
                            primitives: std::mem::take(&mut primitives),
                        });
                        pending = None;
                    } else {
                        pending = Some((l, remaining - 1));
                    }
                } else if f.len() == 2 && f[1] == "0" {
                    current_atom = Some(f[0].parse::<usize>().unwrap() - 1);
                } else {
                    let l = match f[0] {
                        "s" => 0,
                        "p" => 1,
                        "d" => 2,
                        other => panic!("unexpected shell label {other:?}"),
                    };
                    pending = Some((l, f[1].parse().unwrap()));
                }
            }
            "mo" => {
                if let Some(rest) = line.strip_prefix("Sym=") {
                    let _ = rest;
                    coefficients.push(Vec::new());
                } else if let Some(rest) = line.strip_prefix("Occup=") {
                    occupations.push(rest.trim().parse().unwrap());
                } else if let Some(rest) = line.strip_prefix("Spin=") {
                    spins.push(rest.trim().to_string());
                } else if line.starts_with("Ene=") {
                    // energies are not needed by these checks
                } else {
                    let f: Vec<&str> = line.split_whitespace().collect();
                    coefficients
                        .last_mut()
                        .unwrap()
                        .push(f[1].replace('D', "E").parse().unwrap());
                }
            }
            _ => {}
        }
    }
    Parsed {
        atoms,
        shells,
        coefficients,
        occupations,
        spins,
        spherical_d,
    }
}

/// The overlap matrix a Molden reader would build from these shells.
///
/// Written from the **format's** definition rather than reusing
/// `gto::spherical_overlap_matrix`, because the two orders differ and that
/// difference is exactly what needs checking. `[5D]` declares the components in
/// the order `D 0, D+1, D-1, D+2, D-2`, which as real solid harmonics are
/// `z2, xz, yz, x2-y2, xy`; `gto` uses the engine's internal order, which is
/// `x2-y2, xz, z2, yz, xy`. A test that took the order from the writer's own
/// module would agree with the writer about a convention no reader shares.
///
/// The Cartesian integrals themselves are shared -- independence here is about
/// the format and the ordering, not about re-deriving the Gaussian product
/// theorem.
fn molden_order_overlap(shells: &[Shell]) -> Vec<Vec<f64>> {
    // Coefficients on normalised Cartesians (xx, yy, zz, xy, xz, yz).
    const SQRT3_2: f64 = 0.866_025_403_784_438_6;
    let rows = |l: usize| -> Vec<Vec<f64>> {
        match l {
            0 => vec![vec![1.0]],
            1 => vec![
                vec![1.0, 0.0, 0.0],
                vec![0.0, 1.0, 0.0],
                vec![0.0, 0.0, 1.0],
            ],
            2 => vec![
                vec![-0.5, -0.5, 1.0, 0.0, 0.0, 0.0],        // D 0  = z2
                vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0],          // D+1  = xz
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0],          // D-1  = yz
                vec![SQRT3_2, -SQRT3_2, 0.0, 0.0, 0.0, 0.0], // D+2 = x2-y2
                vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0],          // D-2  = xy
            ],
            _ => unreachable!(),
        }
    };
    let sizes: Vec<usize> = shells.iter().map(|s| 2 * s.l + 1).collect();
    let mut offsets = Vec::with_capacity(shells.len());
    let mut total = 0;
    for n in &sizes {
        offsets.push(total);
        total += n;
    }
    let mut out = vec![vec![0.0; total]; total];
    for (ia, a) in shells.iter().enumerate() {
        for (ib, b) in shells.iter().enumerate() {
            let cart = shell_overlap(a, b);
            for (i, ti) in rows(a.l).iter().enumerate() {
                for (j, tj) in rows(b.l).iter().enumerate() {
                    let mut value = 0.0;
                    for (ca, wa) in ti.iter().enumerate() {
                        for (cb, wb) in tj.iter().enumerate() {
                            value += wa * wb * cart[ca][cb];
                        }
                    }
                    out[offsets[ia] + i][offsets[ib] + j] = value;
                }
            }
        }
    }
    out
}

fn write(xyz: &str, method: Method, mode: MoldenCoefficients) -> (String, Molecule) {
    let molecule = Molecule::from_xyz_str(xyz, 0.0).unwrap();
    let text = molden_string(&molecule, method, &NddoOptions::default(), mode)
        .unwrap_or_else(|e| panic!("{method}: {e}"));
    (text, molecule)
}

const CASES: [(&str, &str, Method); 5] = [
    ("water", WATER, Method::Mndo),
    ("water", WATER, Method::Mindo3),
    ("formaldehyde", FORMALDEHYDE, Method::ZindoS),
    ("water", WATER, Method::Cndo2),
    ("h2s", H2S, Method::MndoD),
];

#[test]
fn the_file_describes_the_geometry_and_basis_it_was_given() {
    for (name, xyz, method) in CASES {
        let (text, molecule) = write(xyz, method, MoldenCoefficients::default());
        let parsed = parse(&text);
        assert!(parsed.spherical_d, "{name}/{method}: [5D] is missing");
        assert_eq!(parsed.atoms.len(), molecule.atoms.len(), "{name}/{method}");
        for (index, (z, position)) in parsed.atoms.iter().enumerate() {
            let atom = &molecule.atoms[index];
            assert_eq!(*z, atom.z, "{name}/{method}: atom {index}");
            let d = *position - atom.position;
            assert!(
                d.norm() < 1.0e-9,
                "{name}/{method}: atom {index} moved by {} Bohr",
                d.norm()
            );
        }
        let nao: usize = parsed.shells.iter().map(|s| 2 * s.l + 1).sum();
        assert_eq!(
            nao,
            parsed.coefficients[0].len(),
            "{name}/{method}: the shells give {nao} functions but an MO has {}",
            parsed.coefficients[0].len()
        );
    }
}

#[test]
fn the_written_orbitals_are_orthonormal_over_the_basis_in_the_file() {
    // The claim the default mode exists to make. `C^T S C = I`, with both `C`
    // and `S` taken from the file -- `S` rebuilt from the parsed exponents and
    // contraction coefficients, not from anything the writer kept.
    for (name, xyz, method) in CASES {
        let (text, _) = write(xyz, method, MoldenCoefficients::default());
        let parsed = parse(&text);
        let s = molden_order_overlap(&parsed.shells);
        let n = s.len();
        let orbitals: Vec<&Vec<f64>> = parsed
            .coefficients
            .iter()
            .zip(&parsed.spins)
            .filter(|(_, spin)| spin.as_str() == "Alpha")
            .map(|(c, _)| c)
            .collect();
        let mut worst = 0.0_f64;
        for (i, ci) in orbitals.iter().enumerate() {
            for (j, cj) in orbitals.iter().enumerate() {
                let mut value = 0.0;
                for a in 0..n {
                    for b in 0..n {
                        value += ci[a] * s[a][b] * cj[b];
                    }
                }
                let want = if i == j { 1.0 } else { 0.0 };
                worst = worst.max((value - want).abs());
            }
        }
        assert!(
            worst < 1.0e-8,
            "{name}/{method}: C^T S C differs from the identity by {worst:.3e}"
        );
    }
}

#[test]
fn the_raw_mode_is_not_orthonormal_and_the_two_modes_really_differ() {
    // Without this, the previous test could pass on a file whose basis happened
    // to be near-orthogonal anyway, and the back-transformation would be doing
    // nothing. The raw mode must fail the same check by a wide margin.
    let (text, _) = write(WATER, Method::Mndo, MoldenCoefficients::RawZdo);
    let parsed = parse(&text);
    let s = molden_order_overlap(&parsed.shells);
    let n = s.len();
    let mut worst = 0.0_f64;
    for (i, ci) in parsed.coefficients.iter().enumerate() {
        for (j, cj) in parsed.coefficients.iter().enumerate() {
            let mut value = 0.0;
            for a in 0..n {
                for b in 0..n {
                    value += ci[a] * s[a][b] * cj[b];
                }
            }
            worst = worst.max((value - if i == j { 1.0 } else { 0.0 }).abs());
        }
    }
    assert!(
        worst > 1.0e-3,
        "the raw coefficients are orthonormal over the Gaussian basis to {worst:.3e}, \
         so the back-transformation is not doing anything and the default mode's \
         test proves nothing"
    );
    assert!(
        text.contains("WARNING: raw ZDO coefficients"),
        "a raw-mode file must say so in its title"
    );
}

#[test]
fn the_lowdin_charges_recomputed_from_the_file_are_the_engine_s_charges() {
    // The strongest assertion here, and the only one that checks AO `k` belongs
    // to atom `A`. Orthonormality is invariant under any permutation of the
    // basis functions, so a file that scrambled the AO order would pass the
    // test above and fail this one.
    //
    // For a Lowdin-orthogonalised basis the atomic population is the diagonal
    // of `S^(1/2) P S^(1/2)` summed over an atom's functions, and by
    // construction that equals the engine's own NDO population.
    for (name, xyz, method) in CASES {
        let (text, molecule) = write(xyz, method, MoldenCoefficients::default());
        let parsed = parse(&text);
        let s = molden_order_overlap(&parsed.shells);
        let n = s.len();

        // P = sum_i occ_i c_i c_i^T over every orbital in the file.
        let mut p = vec![vec![0.0; n]; n];
        for (c, occ) in parsed.coefficients.iter().zip(&parsed.occupations) {
            if *occ == 0.0 {
                continue;
            }
            for a in 0..n {
                for b in 0..n {
                    p[a][b] += occ * c[a] * c[b];
                }
            }
        }
        // S^(1/2) P S^(1/2), whose diagonal is the Lowdin population.
        let root = matrix_sqrt(&s);
        let mut populations = vec![0.0; n];
        for a in 0..n {
            let mut value = 0.0;
            for i in 0..n {
                for j in 0..n {
                    value += root[a][i] * p[i][j] * root[j][a];
                }
            }
            populations[a] = value;
        }

        // Which atom each function belongs to, read off the parsed shells.
        let mut owner = Vec::new();
        for shell in &parsed.shells {
            let atom = parsed
                .atoms
                .iter()
                .position(|(_, centre)| (*centre - shell.centre).norm() < 1.0e-9)
                .expect("every shell sits on an atom in the file");
            for _ in 0..(2 * shell.l + 1) {
                owner.push(atom);
            }
        }

        let engine = charges(xyz, method);
        let mut recovered = vec![0.0; parsed.atoms.len()];
        for (ao, atom) in owner.iter().enumerate() {
            recovered[*atom] += populations[ao];
        }
        for (index, atom) in molecule.atoms.iter().enumerate() {
            let core = core_charge(atom.z, method);
            let charge = core - recovered[index];
            assert!(
                (charge - engine[index]).abs() < 1.0e-6,
                "{name}/{method}: atom {index} has charge {charge:.8} from the file \
                 and {:.8} from the engine",
                engine[index]
            );
        }
    }
}

fn charges(xyz: &str, method: Method) -> Vec<f64> {
    let molecule = Molecule::from_xyz_str(xyz, 0.0).unwrap();
    match run_method(&molecule, method, &NddoOptions::default()).unwrap() {
        xndo_rs::CalculationResult::Nddo(r) => r.charges,
        xndo_rs::CalculationResult::Mindo3(r) => r.charges,
        xndo_rs::CalculationResult::CndoIndo(r) => r.charges,
        xndo_rs::CalculationResult::ZindoS(r) => r.charges,
    }
}

fn core_charge(z: u8, method: Method) -> f64 {
    match method {
        Method::Mndo | Method::MndoD => {
            xndo_rs::NddoParameters::cached_for_method(method)
                .unwrap()
                .element(z)
                .unwrap()
                .core_charge
        }
        Method::Mindo3 => xndo_rs::mindo3_element(z).unwrap().core_charge,
        Method::Cndo2 | Method::Indo => xndo_rs::cndo_indo_element(z).unwrap().core_charge,
        Method::ZindoS => {
            xndo_rs::ZindoParameters::cached()
                .unwrap()
                .element(z)
                .unwrap()
                .core_charge
        }
        other => panic!("no core charge for {other}"),
    }
}

/// `S^(1/2)` by eigendecomposition, written here rather than reused so that the
/// check does not lean on the same routine the writer inverted.
fn matrix_sqrt(s: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = s.len();
    let mut m = xndo_rs::Matrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            m[(i, j)] = s[i][j];
        }
    }
    let (values, vectors) = xndo_rs::linalg::symmetric_eigen(&m).unwrap();
    let mut out = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut acc = 0.0;
            for (k, &lambda) in values.iter().enumerate() {
                acc += vectors[(i, k)] * lambda.sqrt() * vectors[(j, k)];
            }
            out[i][j] = acc;
        }
    }
    out
}

#[test]
fn the_sto_expansion_stays_close_to_the_slater_orbital_it_replaces() {
    // Not an equivalence claim -- STO-6G is a fit, and the residual is real.
    // This is a regression guard on the *mapping*: getting `n` or `zeta` wrong,
    // or reading the wrong row of Stewart's table, moves an overlap by O(0.1),
    // while the fit's own error stays three orders below that. A loose ceiling
    // that catches the first without pretending the second is zero.
    let molecule = Molecule::from_xyz_str(WATER, 0.0).unwrap();
    let text = molden_string(
        &molecule,
        Method::Mndo,
        &NddoOptions::default(),
        MoldenCoefficients::default(),
    )
    .unwrap();
    let parsed = parse(&text);
    let gaussian = molden_order_overlap(&parsed.shells);

    // The Slater side, built from the engine's own diatomic overlap routine --
    // a completely different algorithm (analytic A/B integrals in the diatomic
    // frame, then rotated) from the Gaussian product theorem above.
    let params = xndo_rs::NddoParameters::cached_for_method(Method::Mndo).unwrap();
    let mut worst = 0.0_f64;
    let offsets: Vec<usize> = {
        let mut acc = 0;
        molecule
            .atoms
            .iter()
            .map(|atom| {
                let start = acc;
                acc += params.element(atom.z).unwrap().n_orb;
                start
            })
            .collect()
    };
    for (ia, a) in molecule.atoms.iter().enumerate() {
        for (ib, b) in molecule.atoms.iter().enumerate().skip(ia + 1) {
            let ea = params.element(a.z).unwrap();
            let eb = params.element(b.z).unwrap();
            let d = b.position - a.position;
            let block = xndo_rs::overlap::diatom_overlap_spd::<f64>(ea, eb, [d.x, d.y, d.z]);
            for i in 0..ea.n_orb {
                for j in 0..eb.n_orb {
                    let g = gaussian[offsets[ia] + i][offsets[ib] + j];
                    worst = worst.max((g - block[i][j]).abs());
                }
            }
        }
    }
    println!("max |S_STOnG - S_Slater| = {worst:.3e}");
    assert!(
        worst < 5.0e-3,
        "the Gaussian expansion differs from the Slater overlap by {worst:.3e}, \
         which is far beyond the STO-6G fit's own residual: the shell mapping is wrong"
    );
}
