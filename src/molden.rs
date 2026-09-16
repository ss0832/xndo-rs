// SPDX-License-Identifier: GPL-3.0-or-later

//! Molden wavefunction output.
//!
//! # The problem this file exists to solve
//!
//! Writing a Molden file from a semiempirical calculation is not a formatting
//! exercise, and treating it as one produces a file that looks right and is
//! wrong.
//!
//! Every engine in this crate assumes its AO basis is **orthonormal**. That is
//! what zero differential overlap means; it is why the working equations are
//! `F C = C eps` with no `S` in them, and why there is no Pulay term in any
//! gradient here. The MO coefficients an SCF returns are therefore coefficients
//! on an orthonormal set.
//!
//! A Molden file describes real Gaussians, which are **not** orthonormal, and
//! every program that reads one forms the density as `P = C n C^T` over that
//! non-orthogonal set. Hand it the raw coefficients and it will integrate the
//! density to something that is not the electron count, and draw orbitals that
//! are not normalised. MOPAC's own orbital output makes exactly this
//! identification, which is not an argument that it is right.
//!
//! # What is done instead
//!
//! The Slater basis is expanded into Gaussians ([`crate::sto`]), the overlap of
//! *those* Gaussians is computed ([`crate::gto`]), and the coefficients are
//! back-transformed
//!
//! ```text
//!     C_written = S^(-1/2) C_engine
//! ```
//!
//! which is the unique transformation that treats the engine's orthonormal set
//! as the **L??wdin**-orthogonalised version of the real one. Two things then
//! hold exactly, and both are tested:
//!
//! * `C^T S C = I` for the written coefficients, so the density integrates to
//!   the electron count;
//! * the L??wdin populations recomputed from the file reproduce the engine's own
//!   atomic charges. This is the stronger of the two, because it is the one
//!   that fails if AO `k` is attributed to the wrong atom.
//!
//! [`MoldenCoefficients::RawZdo`] writes the untransformed coefficients for
//! comparison with programs that make the identification above. It is not the
//! default, and a file written that way says so in its title.

use std::fmt::Write as _;
use std::path::Path;

use crate::cndo_indo::element as cndo_indo_element;
use crate::constants::{ModelConstants, BOHR_TO_ANGSTROM};
use crate::error::{Result, XndoError};
use crate::gto::{spherical_overlap_matrix, Shell};
use crate::linalg::{symmetric_eigen, Matrix};
use crate::math::Vec3;
use crate::method::Method;
use crate::mindo3::element as mindo3_element;
use crate::orbitals::OrbitalEnergies;
use crate::params::NddoParameters;
use crate::scf::NddoOptions;
use crate::sto::slater_to_gauss;
use crate::system::{z_to_symbol, Molecule};
use crate::xndo::{run_method, CalculationResult};
use crate::zindo::ZindoParameters;

/// Element symbol, or the bare atomic number when there is no name for it --
/// system.rs carries pseudo-atom placeholders above Z=97 that have no symbol.
fn symbol(z: u8) -> String {
    z_to_symbol(z).map_or_else(|| format!("Z{z}"), str::to_string)
}

/// Number of primitives per Slater orbital.
///
/// Six is the largest Stewart tabulated and the only size for which 6s and 6p
/// exist at all -- MNDO reaches bismuth, so those are needed.
pub const DEFAULT_NPRIM: usize = 6;

/// Which coefficients to write.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum MoldenCoefficients {
    /// Back-transform with `S^(-1/2)`, so the written orbitals are orthonormal
    /// over the Gaussians in the file. The default, and the only setting for
    /// which the density in the file is the density the engine computed.
    #[default]
    LowdinBackTransformed,
    /// Write the engine's coefficients unchanged, identifying the orthonormal
    /// ZDO basis with the Gaussian one. For comparison with programs that make
    /// the same identification. The orbitals in such a file are **not**
    /// orthonormal over its own basis.
    RawZdo,
}

/// One atom's shells: principal quantum number and Slater exponent per `l`.
struct AtomShells {
    position: Vec3,
    /// `(l, n, zeta)`, in the order the engine's AO block uses.
    shells: Vec<(usize, u8, f64)>,
}

/// Molden's `[5D]` ordering: `D 0, D+1, D-1, D+2, D-2`, given as offsets into
/// the engine's d block.
///
/// The engine stores d AOs as `d(x2-y2), d(xz), d(z2), d(yz), d(xy)`
/// (`rotations.rs`), so the map is `z2, xz, yz, x2-y2, xy` -- a **pure
/// permutation, with no sign changes**. That was derived by evaluating
/// `rotations.rs`'s own rotation matrices rather than assumed: MOPAC's `d(z2)`
/// carries the standard sign, proportional to `(3z^2 - r^2)/2`.
///
/// Worth stating because the obvious place to copy this from has a sign in it.
/// gfn2-rs's equivalent map carries `-1` on `D 0`, because *its* in-tree `dz2`
/// is stored with the solid harmonic's sign reversed. Copying that here would
/// invert every d(z2) coefficient in the file, and no test of orthonormality
/// would notice, because a sign flip on a basis function preserves `C^T S C`.
const MOLDEN_D_ORDER: [usize; 5] = [2, 1, 3, 0, 4];

/// Shells for one atom under one method.
fn atom_shells(z: u8, position: Vec3, method: Method, nddo: &NddoParameters) -> Result<AtomShells> {
    let unsupported = |what: &str| {
        Err(XndoError::InvalidInput(format!(
            "Molden output is not available for {} under {method}: {what}",
            symbol(z),
        )))
    };
    let shells = match method {
        Method::Mndo | Method::MndoD => {
            let e = nddo.element(z)?;
            let mut v = vec![(0usize, e.n_s, e.zeta_s)];
            if e.n_orb >= 4 {
                v.push((1, e.n_p, e.zeta_p));
            }
            if e.n_orb >= 9 {
                v.push((2, e.n_d, e.zeta_d));
            }
            v
        }
        Method::Mindo3 => {
            // MINDO/3's exponents were fitted against MOPAC7's Bohr radius, and
            // the geometry written to the file is in the crate's internal
            // (CODATA) Bohr, so the exponents have to be restated in it. Same
            // correction as the overlap integral in `mindo3.rs`; leaving it out
            // would write orbitals 1.93e-5 too diffuse.
            let e = mindo3_element(z)?;
            let k = ModelConstants::HISTORICAL.length_scale();
            let mut v = vec![(0usize, e.n_s, e.zeta_s * k)];
            if e.n_orb >= 4 {
                v.push((1, e.n_p, e.zeta_p * k));
            }
            v
        }
        Method::Cndo2 | Method::Indo => {
            // CNDO/2 and INDO use one principal quantum number for every shell
            // on an atom -- `valence_shell` -- which is the implementation's own
            // definition rather than an approximation to something else. Noted
            // in the file's title so a reader is not left to infer it.
            let e = cndo_indo_element(z)?;
            let k = ModelConstants::MOLDS.length_scale();
            let n = e.valence_shell;
            let mut v = vec![(0usize, n, e.zeta_s * k)];
            if e.n_orb >= 4 {
                v.push((1, n, e.zeta_p * k));
            }
            if e.n_orb >= 9 {
                v.push((2, n, e.zeta_d * k));
            }
            v
        }
        Method::ZindoS => {
            let params = ZindoParameters::cached()?;
            let e = params.element(z)?;
            if e.n_orb > 4 {
                // The engine rejects the d branch, so there are no d
                // coefficients to write. Emitting the s and p shells and
                // silently dropping d would produce a file describing a
                // different molecule.
                return unsupported("the ZINDO/S d branch is not implemented");
            }
            // s and p share one exponent in ZINDO/S. Also noted in the title.
            let mut v = vec![(0usize, e.n_s, e.zeta_sp)];
            if e.n_orb >= 4 {
                v.push((1, e.n_p, e.zeta_sp));
            }
            v
        }
        other => return Err(XndoError::InvalidInput(other.execution_error())),
    };
    for &(l, n, zeta) in &shells {
        if n == 0 || (n as usize) <= l || zeta <= 0.0 {
            return unsupported(&format!(
                "its l={l} shell has n={n} and zeta={zeta}, which is not a real orbital"
            ));
        }
    }
    Ok(AtomShells { position, shells })
}

/// Build the Gaussian shells, and the permutation from the engine's AO order to
/// the file's.
fn build_basis(atoms: &[AtomShells]) -> Result<(Vec<Shell>, Vec<usize>)> {
    let mut shells = Vec::new();
    let mut permutation = Vec::new();
    let mut engine_ao = 0usize;
    let mut file_ao = 0usize;
    for atom in atoms {
        let block_start = engine_ao;
        for &(l, n, zeta) in &atom.shells {
            shells.push(Shell {
                centre: atom.position,
                l,
                primitives: slater_to_gauss(DEFAULT_NPRIM, n, l, zeta, true)?,
            });
            match l {
                // s and p are in the same order on both sides: s, then px, py, pz.
                0 | 1 => {
                    let count = 2 * l + 1;
                    for k in 0..count {
                        permutation.push(block_start + engine_offset(atom, l) + k);
                        let _ = file_ao;
                        file_ao += 1;
                    }
                }
                2 => {
                    for &offset in &MOLDEN_D_ORDER {
                        permutation.push(block_start + engine_offset(atom, 2) + offset);
                        file_ao += 1;
                    }
                }
                _ => unreachable!("only s, p and d shells are built"),
            }
        }
        engine_ao += atom
            .shells
            .iter()
            .map(|&(l, _, _)| 2 * l + 1)
            .sum::<usize>();
    }
    Ok((shells, permutation))
}

/// First AO index of shell `l` within an atom's engine block: s at 0, p at 1,
/// d at 4.
fn engine_offset(atom: &AtomShells, l: usize) -> usize {
    atom.shells
        .iter()
        .take_while(|&&(other, _, _)| other < l)
        .map(|&(other, _, _)| 2 * other + 1)
        .sum()
}

/// `S^(-1/2)` for a symmetric positive-definite `S`.
fn inverse_sqrt(s: &Matrix) -> Result<Matrix> {
    let (values, vectors) = symmetric_eigen(s)?;
    let n = s.rows;
    let mut scaled = Matrix::zeros(n, n);
    for (k, &lambda) in values.iter().enumerate() {
        if lambda <= 1.0e-10 {
            return Err(XndoError::InvalidInput(format!(
                "the Gaussian overlap matrix is singular (eigenvalue {lambda:.3e}); the \
                 STO-nG expansion of this basis is linearly dependent"
            )));
        }
        let factor = lambda.powf(-0.5);
        for i in 0..n {
            scaled[(i, k)] = vectors[(i, k)] * factor;
        }
    }
    let mut out = Matrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let mut acc = 0.0;
            for k in 0..n {
                acc += scaled[(i, k)] * vectors[(j, k)];
            }
            out[(i, j)] = acc;
        }
    }
    Ok(out)
}

fn sci(value: f64) -> String {
    format!("{value:.10e}")
}

/// Render a Molden file for `molecule` under `method`.
pub fn molden_string(
    molecule: &Molecule,
    method: Method,
    options: &NddoOptions,
    coefficients: MoldenCoefficients,
) -> Result<String> {
    let result = run_method(molecule, method, options)?;
    let orbitals = result.orbitals();
    let nddo = match method {
        Method::Mndo | Method::MndoD => NddoParameters::cached_for_method(method)?,
        // Only the NDDO branch reads this; the others have their own tables.
        _ => NddoParameters::cached_for_method(Method::Mndo)?,
    };

    let atoms = molecule
        .atoms
        .iter()
        .map(|atom| atom_shells(atom.z, atom.position, method, nddo))
        .collect::<Result<Vec<_>>>()?;
    let (shells, permutation) = build_basis(&atoms)?;
    let overlap = spherical_overlap_matrix(&shells);
    let nao = overlap.rows;
    if nao != orbitals.n_orbitals() {
        return Err(XndoError::InvalidInput(format!(
            "the Gaussian basis has {nao} functions but the engine produced \
             {} orbitals; the shell map is wrong",
            orbitals.n_orbitals()
        )));
    }

    let (alpha, beta) = mo_blocks(&result);
    let transform = match coefficients {
        MoldenCoefficients::LowdinBackTransformed => Some(inverse_sqrt(&overlap)?),
        MoldenCoefficients::RawZdo => None,
    };

    let mut out = String::new();
    out.push_str("[Molden Format]\n[Title]\n");
    let _ = writeln!(out, " xndo-rs {} / {method}", env!("CARGO_PKG_VERSION"));
    match coefficients {
        MoldenCoefficients::LowdinBackTransformed => out.push_str(
            " MO coefficients are back-transformed with S^(-1/2), so they are\n \
             orthonormal over the Gaussians in this file.\n",
        ),
        MoldenCoefficients::RawZdo => out.push_str(
            " WARNING: raw ZDO coefficients. They are orthonormal over the\n \
             engine's assumed basis, NOT over the Gaussians in this file, so the\n \
             density here does not integrate to the electron count.\n",
        ),
    }
    let _ = writeln!(
        out,
        " Slater basis expanded as STO-{DEFAULT_NPRIM}G (Stewart, JCP 52, 431 (1970))."
    );
    match method {
        Method::Cndo2 | Method::Indo => out.push_str(
            " NOTE: CNDO/2 and INDO use one principal quantum number per atom for\n \
             every shell, which is this model's own definition.\n",
        ),
        Method::ZindoS => {
            out.push_str(" NOTE: ZINDO/S uses one Slater exponent for the s and p shells.\n")
        }
        _ => {}
    }

    out.push_str("[Atoms] Angs\n");
    for (index, atom) in molecule.atoms.iter().enumerate() {
        let p = atom.position * BOHR_TO_ANGSTROM;
        let _ = writeln!(
            out,
            " {:<3} {:>5} {:>5} {:>18.10} {:>18.10} {:>18.10}",
            symbol(atom.z),
            index + 1,
            atom.z,
            p.x,
            p.y,
            p.z
        );
    }

    out.push_str("[GTO]\n");
    let mut shell_index = 0usize;
    for (index, atom) in atoms.iter().enumerate() {
        let _ = writeln!(out, "{:>4} 0", index + 1);
        for _ in &atom.shells {
            let shell = &shells[shell_index];
            let label = ["s", "p", "d"][shell.l];
            let _ = writeln!(out, " {label} {:>4} 1.00", shell.primitives.len());
            for primitive in &shell.primitives {
                // Molden wants the contraction coefficients, and
                // `slater_to_gauss` returned them with the Cartesian
                // normalisation already folded in. Dividing it back out is what
                // makes a reader's own normalisation reproduce these functions
                // rather than square the factor.
                let raw = primitive.coefficient / primitive_norm(shell.l, primitive.exponent);
                let _ = writeln!(out, "{:>22} {:>22}", sci(primitive.exponent), sci(raw));
            }
            shell_index += 1;
        }
        out.push('\n');
    }

    out.push_str("[5D]\n");
    out.push_str("[MO]\n");
    write_mo_block(
        &mut out,
        &alpha,
        &orbitals.alpha_ev,
        orbitals.n_alpha,
        &permutation,
        transform.as_ref(),
        "Alpha",
        orbitals.unrestricted_reference(),
    );
    if let (Some(beta), Some(energies)) = (beta, orbitals.beta_ev.as_ref()) {
        write_mo_block(
            &mut out,
            &beta,
            energies,
            orbitals.n_beta,
            &permutation,
            transform.as_ref(),
            "Beta",
            true,
        );
    }
    Ok(out)
}

/// The Cartesian normalisation `slater_to_gauss` folded into each coefficient,
/// so the `[GTO]` block can report the contraction coefficient itself.
fn primitive_norm(l: usize, exponent: f64) -> f64 {
    const DFACTORIAL: [f64; 3] = [1.0, 1.0, 3.0];
    (2.0 / std::f64::consts::PI * exponent).powf(0.75) * (4.0 * exponent).sqrt().powi(l as i32)
        / DFACTORIAL[l].sqrt()
}

#[allow(clippy::too_many_arguments)]
fn write_mo_block(
    out: &mut String,
    coefficients: &Matrix,
    energies: &[f64],
    n_occupied: usize,
    permutation: &[usize],
    transform: Option<&Matrix>,
    spin: &str,
    unrestricted: bool,
) {
    let nao = coefficients.rows;
    let occupancy = if unrestricted { 1.0 } else { 2.0 };
    for orbital in 0..nao {
        let _ = writeln!(out, " Sym= {}a", orbital + 1);
        let _ = writeln!(out, " Ene= {:>18.10}", energies[orbital]);
        let _ = writeln!(out, " Spin= {spin}");
        let _ = writeln!(
            out,
            " Occup= {:>10.6}",
            if orbital < n_occupied { occupancy } else { 0.0 }
        );
        for (file_ao, &engine_ao) in permutation.iter().enumerate() {
            let value = match transform {
                Some(t) => (0..nao)
                    .map(|k| t[(engine_ao, k)] * coefficients[(k, orbital)])
                    .sum(),
                None => coefficients[(engine_ao, orbital)],
            };
            let _ = writeln!(out, "{:>4} {:>22}", file_ao + 1, sci(value));
        }
    }
}

/// Alpha and (for an unrestricted reference) beta MO coefficient matrices.
fn mo_blocks(result: &CalculationResult) -> (Matrix, Option<Matrix>) {
    match result {
        CalculationResult::Nddo(r) => (r.mo_coeff.clone(), r.mo_coeff_beta.clone()),
        CalculationResult::Mindo3(r) => (r.mo_coeff.clone(), r.mo_coeff_beta.clone()),
        CalculationResult::CndoIndo(r) => (r.mo_coeff.clone(), r.mo_coeff_beta.clone()),
        CalculationResult::ZindoS(r) => (r.mo_coeff.clone(), r.mo_coeff_beta.clone()),
    }
}

/// Write a Molden file.
pub fn write_molden(
    path: &Path,
    molecule: &Molecule,
    method: Method,
    options: &NddoOptions,
    coefficients: MoldenCoefficients,
) -> Result<()> {
    let text = molden_string(molecule, method, options, coefficients)?;
    std::fs::write(path, text)
        .map_err(|e| XndoError::InvalidInput(format!("could not write {}: {e}", path.display())))
}

/// The Gaussian overlap matrix a Molden file for this system would be written
/// over, in the file's AO order.
///
/// Public because the tests reconstruct the overlap from the written text and
/// compare, and because a caller doing its own analysis of the file needs the
/// same metric.
pub fn gaussian_overlap(
    molecule: &Molecule,
    method: Method,
) -> Result<(Matrix, Vec<usize>, OrbitalEnergies)> {
    let nddo = match method {
        Method::Mndo | Method::MndoD => NddoParameters::cached_for_method(method)?,
        _ => NddoParameters::cached_for_method(Method::Mndo)?,
    };
    let atoms = molecule
        .atoms
        .iter()
        .map(|atom| atom_shells(atom.z, atom.position, method, nddo))
        .collect::<Result<Vec<_>>>()?;
    let (shells, permutation) = build_basis(&atoms)?;
    let overlap = spherical_overlap_matrix(&shells);
    // Reordered into the file's AO order, which is the order a reader sees.
    let n = overlap.rows;
    let mut reordered = Matrix::zeros(n, n);
    for (i, &ei) in permutation.iter().enumerate() {
        for (j, &ej) in permutation.iter().enumerate() {
            reordered[(i, j)] = overlap[(ei, ej)];
        }
    }
    Ok((
        reordered,
        permutation,
        OrbitalEnergies::restricted(vec![], 0),
    ))
}
