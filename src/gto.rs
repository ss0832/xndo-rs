// SPDX-License-Identifier: GPL-3.0-or-later

//! Overlap integrals over contracted Cartesian Gaussians.
//!
//! This exists for one reason, and it is the reason Molden output is a physics
//! problem rather than a formatting one.
//!
//! Every engine in this crate assumes its AO basis is **orthonormal** -- that is
//! what zero differential overlap means, and it is why the working equations are
//! `F C = C eps` with no `S` in them. So the coefficients an SCF returns are
//! coefficients on an orthonormal set. A Molden file, by contrast, describes
//! real Gaussians, which are **not** orthonormal, and every program that reads
//! one evaluates the density as `P = C n C^T` over that non-orthogonal set. Hand
//! such a program the raw coefficients and it computes a density that does not
//! integrate to the number of electrons and orbitals that are not normalised.
//!
//! The fix needs the overlap matrix of the Gaussians actually written to the
//! file, which is what this module computes. `molden.rs` then back-transforms
//! with `S^(-1/2)`, which is the unique transformation that treats the
//! orthonormal set as the Löwdin-orthogonalised version of the real one -- the
//! choice that keeps each output orbital as close as possible to the input.
//!
//! s, p and d only: no engine here has f functions, and an integral routine for
//! shells the crate cannot produce would be untested code.

use crate::linalg::Matrix;
use crate::math::Vec3;
use crate::sto::PrimitiveGaussian;

/// A contracted Cartesian Gaussian shell: one centre, one angular momentum, and
/// the primitives the Slater orbital was expanded into.
#[derive(Clone, Debug)]
pub struct Shell {
    /// Centre, in Bohr.
    pub centre: Vec3,
    /// 0 = s, 1 = p, 2 = d.
    pub l: usize,
    /// Exponents and contraction coefficients, already carrying the Cartesian
    /// normalisation `slater_to_gauss` applies when asked for it.
    pub primitives: Vec<PrimitiveGaussian>,
}

impl Shell {
    /// Number of Cartesian components: 1, 3 or 6.
    pub fn n_cartesian(&self) -> usize {
        (self.l + 1) * (self.l + 2) / 2
    }
}

/// Cartesian exponent triples for a shell, in the order this module uses.
///
/// s: `1`. p: `x, y, z`. d: `xx, yy, zz, xy, xz, yz` -- which is the order
/// Molden's `[6D]` uses, chosen so that the mapping to `[5D]` in `molden.rs` is
/// written against a fixed, stated convention rather than an implicit one.
pub fn cartesian_powers(l: usize) -> &'static [[usize; 3]] {
    match l {
        0 => &[[0, 0, 0]],
        1 => &[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        2 => &[
            [2, 0, 0],
            [0, 2, 0],
            [0, 0, 2],
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
        ],
        _ => &[],
    }
}

/// `(2k-1)!!` for small `k`, with `(-1)!! = 1`.
fn double_factorial_odd(k: i32) -> f64 {
    let mut out = 1.0;
    let mut i = 2 * k - 1;
    while i > 1 {
        out *= i as f64;
        i -= 2;
    }
    out
}

/// Per-component normalisation, relative to the one already in the contraction.
///
/// A Cartesian Gaussian with powers `(i, j, k)` normalises with
/// `1 / sqrt((2i-1)!! (2j-1)!! (2k-1)!!)`, which is **not** the same for every
/// component of a shell: for d it is `1/sqrt(3)` for `xx` and `1` for `xy`.
/// `slater_to_gauss` folds in the `xx`-type factor `1/sqrt((2l-1)!!)` once for
/// the whole shell, which is the convention basis-set files and Molden use --
/// one set of contraction coefficients per shell, with the reader expected to
/// supply the rest. This supplies the rest:
///
/// ```text
///   sqrt( (2l-1)!! / ((2i-1)!! (2j-1)!! (2k-1)!!) )
/// ```
///
/// which is 1 for `xx`, `yy`, `zz` and `sqrt(3)` for `xy`, `xz`, `yz`. Leaving
/// it out makes the off-diagonal d components come back with a self-overlap of
/// exactly 1/3, which is what happened.
fn cartesian_normalisation(l: usize, powers: [usize; 3]) -> f64 {
    let denominator: f64 = powers
        .iter()
        .map(|&p| double_factorial_odd(p as i32))
        .product();
    (double_factorial_odd(l as i32) / denominator).sqrt()
}

/// One Cartesian direction of the overlap of two primitive Gaussians.
///
/// The Gaussian product theorem says two Gaussians multiply to a third centred
/// between them, so the integral factorises over x, y and z and each factor is a
/// polynomial times `sqrt(pi/p)`. This evaluates one factor directly from the
/// binomial expansion about the product centre:
///
/// ```text
///   S_i = sum_{t even} C(t; i, j, PA, PB) * (2t-1)!! / (2p)^t * sqrt(pi/p)
/// ```
///
/// Odd `t` vanishes because the remaining integrand is odd about the product
/// centre. Written out rather than recursed: for l <= 2 the sums have at most
/// three terms, and the closed form is easier to check against the algebra than
/// an Obara-Saika recursion would be.
fn overlap_1d(i: usize, j: usize, pa: f64, pb: f64, p: f64) -> f64 {
    let mut total = 0.0;
    for t in 0..=(i + j) {
        if t % 2 != 0 {
            continue;
        }
        // Binomial prefactor: the coefficient of (x-P)^t in (x-A)^i (x-B)^j.
        let mut f = 0.0;
        for k in 0..=i.min(t) {
            let m = t - k;
            if m > j {
                continue;
            }
            f +=
                binomial(i, k) * binomial(j, m) * pa.powi((i - k) as i32) * pb.powi((j - m) as i32);
        }
        total += f * double_factorial_odd(t as i32 / 2) / (2.0 * p).powi(t as i32 / 2);
    }
    total * (std::f64::consts::PI / p).sqrt()
}

fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let mut out = 1.0;
    for step in 0..k {
        out = out * (n - step) as f64 / (step + 1) as f64;
    }
    out
}

/// Overlap of two contracted Cartesian Gaussian shells: an
/// `n_cartesian(a) x n_cartesian(b)` block.
pub fn shell_overlap(a: &Shell, b: &Shell) -> Vec<Vec<f64>> {
    let pa_powers = cartesian_powers(a.l);
    let pb_powers = cartesian_powers(b.l);
    let d = b.centre - a.centre;
    let r2 = d.x * d.x + d.y * d.y + d.z * d.z;
    let mut block = vec![vec![0.0; pb_powers.len()]; pa_powers.len()];
    for ga in &a.primitives {
        for gb in &b.primitives {
            let p = ga.exponent + gb.exponent;
            let mu = ga.exponent * gb.exponent / p;
            let pre = ga.coefficient * gb.coefficient * (-mu * r2).exp();
            // Product centre, measured from each of the two centres.
            let pv = [
                (ga.exponent * a.centre.x + gb.exponent * b.centre.x) / p,
                (ga.exponent * a.centre.y + gb.exponent * b.centre.y) / p,
                (ga.exponent * a.centre.z + gb.exponent * b.centre.z) / p,
            ];
            let av = [a.centre.x, a.centre.y, a.centre.z];
            let bv = [b.centre.x, b.centre.y, b.centre.z];
            for (ia, powers_a) in pa_powers.iter().enumerate() {
                for (ib, powers_b) in pb_powers.iter().enumerate() {
                    let mut value = pre
                        * cartesian_normalisation(a.l, *powers_a)
                        * cartesian_normalisation(b.l, *powers_b);
                    for axis in 0..3 {
                        value *= overlap_1d(
                            powers_a[axis],
                            powers_b[axis],
                            pv[axis] - av[axis],
                            pv[axis] - bv[axis],
                            p,
                        );
                    }
                    block[ia][ib] += value;
                }
            }
        }
    }
    block
}

/// Real solid harmonics as combinations of this module's **normalised**
/// Cartesian components, in the *engine's* d order: `x2-y2, xz, z2, yz, xy`.
///
/// Constants, not a fit. For a shell whose Cartesian components are normalised,
/// the on-centre overlaps are exactly `<xx|xx> = 1`, `<xx|yy> = <xx|zz> = 1/3`
/// and zero between an `xx`-type and an `xy`-type, whatever the contraction is.
/// From those:
///
/// * `z2 = (2z^2 - x^2 - y^2)/2` has norm `1/4 + 1/4 + 1 + 2(1/4)(1/3)
///   - 2(1/2)(1/3) - 2(1/2)(1/3) = 1`, so it is already normalised;
/// * `x2-y2` has norm `1 + 1 - 2/3 = 4/3`, hence the `sqrt(3)/2`.
///
/// Deriving the factor from each shell's own computed overlap instead would
/// also *renormalise away the STO-nG fit's residual error* -- 1.2e-10 on a
/// 6s contraction -- which is a real property of the basis being written, not
/// noise to be corrected. It also made s and p disagree with the Cartesian
/// routine, which they must not: for them this is the identity.
fn spherical_rows(l: usize) -> &'static [&'static [f64]] {
    const SQRT3_2: f64 = 0.866_025_403_784_438_6; // sqrt(3)/2
    match l {
        0 => &[&[1.0]],
        1 => &[&[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0], &[0.0, 0.0, 1.0]],
        //       xx       yy      zz   xy   xz   yz
        2 => &[
            &[SQRT3_2, -SQRT3_2, 0.0, 0.0, 0.0, 0.0], // x2 - y2
            &[0.0, 0.0, 0.0, 0.0, 1.0, 0.0],          // xz
            &[-0.5, -0.5, 1.0, 0.0, 0.0, 0.0],        // z2
            &[0.0, 0.0, 0.0, 0.0, 0.0, 1.0],          // yz
            &[0.0, 0.0, 0.0, 1.0, 0.0, 0.0],          // xy
        ],
        _ => &[],
    }
}

/// Overlap matrix over the **real solid harmonic** (spherical) functions, which
/// is the basis every engine here actually uses: five d functions, not six.
///
/// The Cartesian block is computed first and contracted on both sides with
/// [`spherical_transform`].
pub fn spherical_overlap_matrix(shells: &[Shell]) -> Matrix {
    let transforms: Vec<&'static [&'static [f64]]> =
        shells.iter().map(|s| spherical_rows(s.l)).collect();
    let sizes: Vec<usize> = shells.iter().map(|s| 2 * s.l + 1).collect();
    let offsets: Vec<usize> = sizes
        .iter()
        .scan(0, |acc, n| {
            let start = *acc;
            *acc += n;
            Some(start)
        })
        .collect();
    let total: usize = sizes.iter().sum();
    let mut out = Matrix::zeros(total, total);
    for (ia, a) in shells.iter().enumerate() {
        for (ib, b) in shells.iter().enumerate().skip(ia) {
            let cart = shell_overlap(a, b);
            for (i, ti) in transforms[ia].iter().enumerate() {
                for (j, tj) in transforms[ib].iter().enumerate() {
                    let mut value = 0.0;
                    for (ca, wa) in ti.iter().enumerate() {
                        if *wa == 0.0 {
                            continue;
                        }
                        for (cb, wb) in tj.iter().enumerate() {
                            value += wa * wb * cart[ca][cb];
                        }
                    }
                    out[(offsets[ia] + i, offsets[ib] + j)] = value;
                    out[(offsets[ib] + j, offsets[ia] + i)] = value;
                }
            }
        }
    }
    out
}

/// Overlap matrix of a list of shells, in the Cartesian AO order the shells
/// imply.
pub fn overlap_matrix(shells: &[Shell]) -> Matrix {
    let offsets: Vec<usize> = shells
        .iter()
        .scan(0, |acc, s| {
            let start = *acc;
            *acc += s.n_cartesian();
            Some(start)
        })
        .collect();
    let n =
        offsets.last().copied().unwrap_or(0) + shells.last().map(|s| s.n_cartesian()).unwrap_or(0);
    let mut out = Matrix::zeros(n, n);
    for (ia, a) in shells.iter().enumerate() {
        for (ib, b) in shells.iter().enumerate().skip(ia) {
            let block = shell_overlap(a, b);
            for (i, row) in block.iter().enumerate() {
                for (j, value) in row.iter().enumerate() {
                    out[(offsets[ia] + i, offsets[ib] + j)] = *value;
                    out[(offsets[ib] + j, offsets[ia] + i)] = *value;
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sto::slater_to_gauss;

    fn shell(l: usize, n: u8, zeta: f64, centre: Vec3) -> Shell {
        Shell {
            centre,
            l,
            primitives: slater_to_gauss(6, n, l, zeta, true).unwrap(),
        }
    }

    #[test]
    fn a_normalised_shell_overlaps_itself_to_one() {
        // The contraction `slater_to_gauss` returns is normalised, so every
        // diagonal element of the overlap must be 1. Not to within a tolerance
        // that hides a missing factor: a dropped sqrt(3) on a d function shows
        // up here as 0.333 or 3.0, and a wrong double factorial as a small
        // integer ratio. The bound is tight enough to catch those and loose
        // enough for the STO-6G fit's own residual.
        for (l, n, zeta) in [
            (0usize, 1u8, 1.3),
            (0, 2, 1.8),
            (1, 2, 1.7),
            (1, 3, 1.4),
            (2, 3, 2.2),
        ] {
            let s = shell(l, n, zeta, Vec3::zero());
            let block = shell_overlap(&s, &s);
            for (i, row) in block.iter().enumerate() {
                assert!(
                    (row[i] - 1.0).abs() < 1.0e-9,
                    "l={l} n={n}: component {i} has self-overlap {}, not 1",
                    row[i]
                );
            }
        }
    }

    #[test]
    fn the_cartesian_components_of_a_shell_are_orthogonal_at_one_centre() {
        // On the same centre, two Cartesian components of the same shell are
        // orthogonal whenever their exponent difference is odd in any axis --
        // x vs y, xy vs xz, and so on. `xx` and `yy` are *not* orthogonal, and
        // this test would be wrong to claim they are: that non-orthogonality is
        // exactly why the six Cartesian d functions carry a spurious s
        // component, and why Molden has a `[5D]` form at all.
        let s = shell(2, 3, 2.2, Vec3::zero());
        let block = shell_overlap(&s, &s);
        let powers = cartesian_powers(2);
        for (i, pa) in powers.iter().enumerate() {
            for (j, pb) in powers.iter().enumerate() {
                let odd = (0..3).any(|axis| (pa[axis] + pb[axis]) % 2 != 0);
                if odd {
                    assert!(
                        block[i][j].abs() < 1.0e-12,
                        "{pa:?} and {pb:?} differ by an odd power and must not overlap, got {}",
                        block[i][j]
                    );
                }
            }
        }
        // And the statement the comment above makes, asserted rather than
        // assumed: xx and yy do overlap.
        assert!(block[0][1].abs() > 0.1);
    }

    #[test]
    fn two_s_functions_reproduce_the_closed_form_gaussian_overlap() {
        // An independent check of the machinery against algebra that can be
        // done by hand: for two normalised s primitives the overlap is
        //     (2 sqrt(a b) / (a + b))^{3/2} exp(-a b R^2 / (a + b)).
        // Different route entirely -- no binomial expansion, no product centre.
        let (a, b, r) = (0.7_f64, 1.9_f64, 1.3_f64);
        let norm = |e: f64| (2.0 * e / std::f64::consts::PI).powf(0.75);
        let left = Shell {
            centre: Vec3::zero(),
            l: 0,
            primitives: vec![PrimitiveGaussian {
                exponent: a,
                coefficient: norm(a),
            }],
        };
        let right = Shell {
            centre: Vec3::new(0.0, 0.0, r),
            l: 0,
            primitives: vec![PrimitiveGaussian {
                exponent: b,
                coefficient: norm(b),
            }],
        };
        let got = shell_overlap(&left, &right)[0][0];
        let want = (2.0 * (a * b).sqrt() / (a + b)).powf(1.5) * (-a * b * r * r / (a + b)).exp();
        assert!(
            (got - want).abs() < 1.0e-12,
            "got {got}, closed form {want}"
        );
    }

    #[test]
    fn the_five_spherical_d_functions_are_orthonormal_on_one_centre() {
        // The check that the *shapes* in `raw_spherical_rows` are right, which
        // normalising the rows cannot tell you: five functions each normalised
        // to one are still wrong if they are not mutually orthogonal. `z2` is
        // the row that would fail -- it is orthogonal to `x2-y2` only because
        // its `xx` and `yy` coefficients are both exactly -1/2, and the six
        // Cartesian d functions are famously *not* orthogonal among themselves.
        // The two halves get different bounds, and the difference is the point.
        // Orthogonality between two spherical d functions is a symmetry
        // statement and holds to machine precision whatever the radial part is.
        // Normalisation is not: it inherits the STO-6G least-squares fit's own
        // residual, about 2e-10 here, because `spherical_rows` is analytic and
        // deliberately does not renormalise that away. A single tolerance loose
        // enough for the diagonal would stop checking the off-diagonal at all.
        let s = shell(2, 3, 2.2, Vec3::zero());
        let m = spherical_overlap_matrix(std::slice::from_ref(&s));
        assert_eq!(m.rows, 5);
        for i in 0..5 {
            assert!(
                (m[(i, i)] - 1.0).abs() < 1.0e-8,
                "spherical d function {i} has norm {}",
                m[(i, i)]
            );
            for j in 0..5 {
                if i == j {
                    continue;
                }
                assert!(
                    m[(i, j)].abs() < 1.0e-12,
                    "spherical d functions {i} and {j} overlap by {}",
                    m[(i, j)]
                );
            }
        }
    }

    #[test]
    fn spherical_and_cartesian_agree_wherever_there_is_no_d_shell() {
        // s and p have nothing to transform, so the two routines must return
        // the same matrix element for element. A difference would mean the
        // spherical path had picked up a stray normalisation on the way past.
        let shells = vec![
            shell(0, 2, 1.8, Vec3::zero()),
            shell(1, 2, 1.7, Vec3::zero()),
            shell(0, 1, 1.3, Vec3::new(0.0, 0.0, 2.0)),
        ];
        let a = overlap_matrix(&shells);
        let b = spherical_overlap_matrix(&shells);
        assert_eq!(a.rows, b.rows);
        let mut worst = 0.0_f64;
        let mut where_ = (0, 0);
        for i in 0..a.rows {
            for j in 0..a.cols {
                let d = (a[(i, j)] - b[(i, j)]).abs();
                if d > worst {
                    worst = d;
                    where_ = (i, j);
                }
            }
        }
        assert!(
            worst < 1.0e-14,
            "worst difference {worst:.3e} at {where_:?}: {} vs {}",
            a[where_],
            b[where_]
        );
    }

    #[test]
    fn the_overlap_matrix_is_symmetric_and_unit_diagonal() {
        let shells = vec![
            shell(0, 2, 1.8, Vec3::zero()),
            shell(1, 2, 1.7, Vec3::zero()),
            shell(0, 1, 1.3, Vec3::new(0.0, 0.0, 2.0)),
        ];
        let s = overlap_matrix(&shells);
        assert_eq!(s.rows, 5);
        for i in 0..s.rows {
            assert!((s[(i, i)] - 1.0).abs() < 1.0e-9);
            for j in 0..s.cols {
                assert!((s[(i, j)] - s[(j, i)]).abs() < 1.0e-14);
            }
        }
    }
}
