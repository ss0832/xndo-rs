// SPDX-License-Identifier: GPL-3.0-or-later

//! SCF convergence acceleration, shared by every engine in the crate.
//!
//! All four engines work in an orthonormal AO basis, so the Pulay error vector
//! is simply the commutator `[F, P]` with no `S` anywhere, and one
//! implementation serves the NDDO driver and the three ZDO drivers alike.
//!
//! The policy the drivers use is A-DIIS first, CDIIS once the commutator norm
//! falls below a switch threshold: A-DIIS (Hu & Yang 2010) is robust far from
//! convergence where CDIIS can extrapolate onto a saddle, and CDIIS converges
//! faster once the error is small.
//!
//! This module was extracted from `scf.rs`, which had the only implementation;
//! the ZDO engines had no accelerator at all, which is what let ZINDO/S settle
//! on a symmetry-broken solution for BH and fail to converge on benzene (see
//! `tests/data/ORACLE_NOTES.md` items 18 and 19).

use crate::error::Result;
use crate::linalg::Matrix;

/// SCF charge-convergence accelerator.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ScfAccelerator {
    /// No extrapolation (plain iteration).
    None,
    /// Pulay CDIIS on the `[F,P]` commutator throughout.
    Cdiis,
    /// **A-DIIS** (Hu & Yang, *J. Chem. Phys.* **132**, 054109 (2010)) while far from
    /// convergence, switching to CDIIS once the commutator error drops below a threshold
    /// the robust hybrid recommended for hard cases (radicals, small gaps, poor guesses).
    AdiisCdiis,
}

/// Superposition of atomic densities for an s/p valence basis.
///
/// `atoms` is `(first AO index, AO count, core charge)` per atom. Each atom's
/// valence electrons are placed spherically: up to two in `s`, the remainder
/// split equally over the three `p` functions, so the guess carries no
/// orientation of its own. That matters more than it sounds -- the AO basis is
/// fixed to the global axes while the molecule is not, so a guess that favours
/// particular p functions starts the iteration from a physically different
/// place depending on how the molecule happens to be oriented.
///
/// The trace is the valence electron count, so the result is a total density
/// and can be used as the RHF starting `P` directly.
pub(crate) fn sad_density_sp(nao: usize, atoms: &[(usize, usize, f64)]) -> Matrix {
    let mut p = Matrix::zeros(nao, nao);
    for &(offset, n_orb, core_charge) in atoms {
        if n_orb == 0 {
            continue;
        }
        let zc = core_charge.max(0.0);
        if n_orb == 1 {
            p[(offset, offset)] = zc.min(2.0);
            continue;
        }
        let ns = zc.min(2.0);
        p[(offset, offset)] = ns;
        let np = (zc - ns).max(0.0) / 3.0;
        for q in 1..4.min(n_orb) {
            p[(offset + q, offset + q)] = np;
        }
        // Any d functions share whatever is left, again spherically.
        if n_orb > 4 {
            let placed = ns + np * 3.0;
            let nd = (zc - placed).max(0.0) / (n_orb - 4) as f64;
            for q in 4..n_orb {
                p[(offset + q, offset + q)] = nd;
            }
        }
    }
    p
}

/// Valence density with every AO of an atom equally occupied.
///
/// The companion to [`sad_density_sp`], and the other end of the same axis. SAD
/// is the free *atom*: carbon gets `s^2 p^{2/3}`. This is the fully hybridised
/// limit: carbon gets 1.0 in each of s, px, py, pz. A bonded carbon is between
/// the two (MOPAC7 puts 1.34 in methyl's carbon 2s), and neither end is a
/// reliable sole starting point for an open shell -- see
/// [`open_shell_starts`] for why that matters and what is done about it.
pub(crate) fn uniform_valence_density(nao: usize, atoms: &[(usize, usize, f64)]) -> Matrix {
    let mut p = Matrix::zeros(nao, nao);
    for &(offset, n_orb, core_charge) in atoms {
        if n_orb == 0 {
            continue;
        }
        let share = (core_charge.max(0.0) / n_orb as f64).min(2.0);
        for q in 0..n_orb {
            p[(offset + q, offset + q)] = share;
        }
    }
    p
}

/// The starting densities an unrestricted SCF should try, best first.
///
/// A closed shell gets exactly one: `P/2` in each channel from the SAD density,
/// so nothing about a restricted or spin-symmetric run changes.
///
/// An open shell gets two, and the caller is expected to run them all and keep
/// the variationally **lowest** converged solution. That is not belt and
/// braces; it is load-bearing, because of where an open-shell SCF makes its
/// most consequential decision:
///
/// > which orbital the unpaired electron occupies is decided by the aufbau
/// > ordering of the *first* Fock matrix -- the one built from the guess, the
/// > worst density of the entire run -- and once that electron is placed,
/// > self-consistency defends the choice.
///
/// The frontier ordering is exactly what a guess gets wrong. SAD is the free
/// atom, so carbon arrives with 2.0 electrons in its 2s where a bonded carbon
/// has about 1.34. With a one-centre `gsp` of 11.47 eV that surplus pushes the
/// carbon `p_z` up by several eV, and in planar methyl it lifts `p_z` above the
/// C-H antibonding `a1*`: the odd electron goes into `a1*` at iteration one and
/// stays. The result is a real, aufbau-satisfying, fully converged SCF solution
/// 3.07 eV above the right one, with zero spin population on `p_z` where MOPAC7
/// has exactly 1.0. No accelerator and no damping factor escapes it -- they all
/// converge faithfully to the same wrong fixed point, which is the point: this
/// is not a convergence failure and cannot be fixed by converging harder.
///
/// The second start comes from the opposite end of the same axis
/// ([`uniform_valence_density`]), so it does not carry the free-atom s/p
/// imbalance and orders the frontier differently. Running both and keeping the
/// lower is monotone -- it can never return a worse answer than the first start
/// alone -- and costs one extra SCF, and only for open shells.
pub(crate) fn open_shell_starts(
    sad: &Matrix,
    uniform: &Matrix,
    n_alpha: usize,
    n_beta: usize,
) -> Vec<(Matrix, Matrix)> {
    let first = split_by_spin(sad, n_alpha, n_beta);
    if n_alpha == n_beta {
        return vec![first];
    }
    vec![first, split_by_spin(uniform, n_alpha, n_beta)]
}

/// Run an SCF from each start and keep the variationally lowest solution.
///
/// The companion to [`open_shell_starts`]. Ties go to the earlier start, with a
/// margin far wider than the noise between two runs that reached the same fixed
/// point and far narrower than the gap between two different ones, so the
/// answer does not depend on which of two identical solutions came first.
///
/// A start that fails to converge is not fatal while another succeeds; the
/// first error is kept and returned only if every start fails.
pub(crate) fn lowest_solution<R>(
    starts: Vec<(Matrix, Matrix)>,
    mut solve: impl FnMut((Matrix, Matrix)) -> Result<R>,
    energy: impl Fn(&R) -> f64,
) -> Result<R> {
    debug_assert!(!starts.is_empty());
    let mut best: Option<R> = None;
    let mut first_error = None;
    for start in starts {
        match solve(start) {
            Ok(candidate) => {
                let better = match &best {
                    Some(current) => energy(&candidate) < energy(current) - 1.0e-9,
                    None => true,
                };
                if better {
                    best = Some(candidate);
                }
            }
            Err(error) => {
                if first_error.is_none() {
                    first_error = Some(error);
                }
            }
        }
    }
    match (best, first_error) {
        (Some(result), _) => Ok(result),
        (None, Some(error)) => Err(error),
        (None, None) => unreachable!("lowest_solution was given no starts"),
    }
}

/// Split a total density into spin channels, putting the excess spin where
/// there is room for it.
///
/// The closed-shell part is `P/2` in each channel; the `n_alpha - n_beta`
/// excess electrons go into alpha, spread over the AOs in proportion to
///
/// ```text
///     cap_mu = min(P_mu, 2 - P_mu)
/// ```
///
/// which is exactly how much spin polarisation AO `mu` can carry: an empty AO
/// has nothing to polarise, a doubly occupied one has no room, and the most an
/// AO in between can hold is whichever of those two is smaller. Weighting by it
/// automatically keeps `0 <= P_beta` and `P_alpha <= 1` per AO without a
/// separate clamp.
///
/// **Why not in proportion to `P` itself.** That is what this did before, and
/// it is backwards: it puts the *largest* share of the unpaired spin on the
/// *most* occupied AO. For a carbon in a planar radical that is the doubly
/// occupied 2s, so the guess arrives carrying a hole in the sigma system, and
/// the iteration converges to an aufbau solution with the odd electron in a
/// C-H antibonding orbital instead of the out-of-plane p. Planar methyl came
/// out 3.07 eV above the MOPAC7 MINDO/3 answer that way, with a p_z spin
/// population of exactly zero where the oracle has 1.0. Nothing downstream can
/// recover from it -- every accelerator and every damping factor reached the
/// same wrong fixed point, because it is a real fixed point, just not the
/// lowest one. Distorting the molecule out of plane, which lets the p mix and
/// destroys the trap, made the same code agree to 6.7e-4 eV.
///
/// The rule has no free parameters and is not fitted to that case: applied to a
/// free atom it reproduces Hund's rule exactly. Carbon gets `s(1,1) p(2/3,0)`,
/// nitrogen `s(1,1) p(1,0)`, oxygen `s(1,1) p(1,1/3)` -- the spherically
/// averaged maximum-multiplicity configurations, which is what "put the excess
/// spin where there are holes" means for an isolated atom.
///
/// A closed-shell system is untouched: `n_alpha == n_beta` gives `P/2` in each
/// channel, the same as before.
pub(crate) fn split_by_spin(total: &Matrix, n_alpha: usize, n_beta: usize) -> (Matrix, Matrix) {
    let mut pa = total.clone();
    for value in pa.as_mut_slice() {
        *value *= 0.5;
    }
    let mut pb = pa.clone();
    let excess = n_alpha as f64 - n_beta as f64;
    if excess.abs() < 1.0e-12 {
        return (pa, pb);
    }
    let n = total.rows;
    let cap: Vec<f64> = (0..n)
        .map(|i| total[(i, i)].min(2.0 - total[(i, i)]).max(0.0))
        .collect();
    let capacity: f64 = cap.iter().sum();
    if capacity <= 0.0 {
        // No AO can carry spin (every one empty or exactly doubly occupied).
        // Nothing better is available, so leave the channels equal and let the
        // first diagonalisation separate them.
        return (pa, pb);
    }
    // Clamped at 1: if the guess cannot hold all the requested spin, carry what
    // it can rather than pushing an AO past single occupancy.
    let fill = (excess.abs() / capacity).min(1.0);
    let sign = excess.signum();
    for i in 0..n {
        let s = sign * fill * cap[i];
        pa[(i, i)] += 0.5 * s;
        pb[(i, i)] -= 0.5 * s;
    }
    (pa, pb)
}

/// One SCF convergence accelerator, driving the A-DIIS then CDIIS policy.
///
/// Every engine follows the same recipe, so it lives here once: feed the plain
/// Fock built from the current density, get back the Fock to diagonalise.
///
/// A-DIIS (Hu & Yang 2010) runs while the commutator norm is above `switch`,
/// because far from convergence CDIIS can extrapolate onto a saddle; CDIIS takes
/// over below it, because it converges faster once the error is small.
pub(crate) struct Accelerator {
    mode: ScfAccelerator,
    switch: f64,
    keep_densities: bool,
    history: AccelHistory,
}

impl Accelerator {
    /// Accelerator for a single (closed-shell) Fock/density pair.
    pub(crate) fn new(nao: usize, mode: ScfAccelerator, switch: f64, budget_mb: usize) -> Self {
        Self::build(nao, 1, mode, switch, budget_mb)
    }

    /// Accelerator for a spin-unrestricted pair of channels.
    ///
    /// The two channels go into **one** history with **one** set of
    /// coefficients, because they are not independent: `F_alpha` is built from
    /// the total density and so depends on `P_beta` as well. Extrapolating each
    /// channel against its own history would drive them to points that are not
    /// mutually consistent. Stacking `[F_a; F_b]`, `[P_a; P_b]` and
    /// `[[F_a,P_a]; [F_b,P_b]]` and running the ordinary machinery on the stack
    /// is the textbook UHF-DIIS error vector.
    ///
    /// A-DIIS works on the stack unchanged, and needs no new energy functional
    /// to do so: the Frobenius product of the stacks is
    /// `Tr[P_a F_a] + Tr[P_b F_b]`, which is exactly the pairing the UHF energy
    /// uses, just as the single-channel product is the pairing RHF uses. (An
    /// earlier comment in `scf.rs` claimed otherwise and left UHF unaccelerated
    /// by default.)
    pub(crate) fn new_uhf(nao: usize, mode: ScfAccelerator, switch: f64, budget_mb: usize) -> Self {
        Self::build(nao, 2, mode, switch, budget_mb)
    }

    /// `channels` is 1 for RHF, 2 for UHF; a slot holds that many Fock matrices,
    /// that many error matrices, and (for A-DIIS) that many densities.
    ///
    /// When the memory budget cannot hold even a depth-2 A-DIIS history, this
    /// degrades to CDIIS rather than blowing the budget: the weaker accelerator
    /// reaches the same fixed point.
    fn build(
        nao: usize,
        channels: usize,
        mode: ScfAccelerator,
        switch: f64,
        budget_mb: usize,
    ) -> Self {
        let mut keep_densities = mode == ScfAccelerator::AdiisCdiis;
        if keep_densities && !diis_depth_fits(nao, 3 * channels, budget_mb) {
            keep_densities = false;
        }
        let mode = if mode == ScfAccelerator::AdiisCdiis && !keep_densities {
            ScfAccelerator::Cdiis
        } else {
            mode
        };
        let depth = diis_depth(
            nao,
            if keep_densities { 3 } else { 2 } * channels,
            budget_mb,
        );
        Self {
            mode,
            switch,
            keep_densities,
            history: AccelHistory::new(depth, keep_densities),
        }
    }

    pub(crate) fn is_active(&self) -> bool {
        self.mode != ScfAccelerator::None
    }

    /// Effective history depth after the memory budget is applied; reported by
    /// the `XNDO_TIMING` trace.
    pub(crate) fn depth(&self) -> usize {
        self.history.depth
    }

    /// Consume the plain Fock and the density it was built from; return the Fock
    /// to diagonalise.
    pub(crate) fn step(&mut self, fock: Matrix, density: &Matrix) -> Matrix {
        if self.mode == ScfAccelerator::None {
            return fock;
        }
        let error = commutator(&fock, density);
        self.extrapolate(fock, error, density.clone())
    }

    /// The spin-unrestricted step: both channels, one extrapolation.
    pub(crate) fn step_uhf(
        &mut self,
        fa: Matrix,
        fb: Matrix,
        pa: &Matrix,
        pb: &Matrix,
    ) -> (Matrix, Matrix) {
        if self.mode == ScfAccelerator::None {
            return (fa, fb);
        }
        let nao = fa.rows;
        let error = stack_channels(&commutator(&fa, pa), &commutator(&fb, pb));
        let fock = stack_channels(&fa, &fb);
        let density = stack_channels(pa, pb);
        let out = self.extrapolate(fock, error, density);
        split_channels(&out, nao)
    }

    /// Shared tail of `step` / `step_uhf`: push the iterate and extrapolate.
    fn extrapolate(&mut self, fock: Matrix, error: Matrix, density: Matrix) -> Matrix {
        let error_norm = self
            .history
            .push(fock, error, self.keep_densities.then_some(density))
            .sqrt();
        let extrapolated = if self.mode == ScfAccelerator::AdiisCdiis && error_norm > self.switch {
            self.history.adiis()
        } else {
            self.history.cdiis()
        };
        // No usable extrapolation (a singular Gram matrix, or a single-entry
        // history): fall back to the newest plain Fock.
        extrapolated.unwrap_or_else(|| self.history.focks[self.history.len() - 1].clone())
    }
}

/// Stack two `nao x nao` spin channels into one `2nao x nao` matrix.
fn stack_channels(a: &Matrix, b: &Matrix) -> Matrix {
    let nao = a.rows;
    let mut out = Matrix::zeros(2 * nao, a.cols);
    out.as_mut_slice()[..nao * a.cols].copy_from_slice(a.as_slice());
    out.as_mut_slice()[nao * a.cols..].copy_from_slice(b.as_slice());
    out
}

/// Inverse of [`stack_channels`].
fn split_channels(stacked: &Matrix, nao: usize) -> (Matrix, Matrix) {
    let cols = stacked.cols;
    let mut a = Matrix::zeros(nao, cols);
    let mut b = Matrix::zeros(nao, cols);
    a.as_mut_slice()
        .copy_from_slice(&stacked.as_slice()[..nao * cols]);
    b.as_mut_slice()
        .copy_from_slice(&stacked.as_slice()[nao * cols..]);
    (a, b)
}

/// Hard cap on the DIIS history, independent of the memory budget.
pub(crate) const MAX_DIIS_DEPTH: usize = 8;

/// Effective DIIS depth under the configured memory budget.
///
/// A depth-`k` history holds `k` Fock matrices, `k` error matrices and (for
/// A-DIIS) `k` densities  `3 k nao2` doubles. At `nao = 4800` the full depth
/// is 4.4 GiB, which is an out-of-memory abort on most machines; trimming the
/// depth costs a few extra iterations instead. Never drops below 2, the
/// minimum for any extrapolation.
pub(crate) fn diis_depth(nao: usize, matrices_per_slot: usize, budget_mb: usize) -> usize {
    if budget_mb == 0 {
        return MAX_DIIS_DEPTH;
    }
    let per_slot = slot_bytes(nao, matrices_per_slot);
    if per_slot == 0 {
        return MAX_DIIS_DEPTH;
    }
    let budget = budget_mb.saturating_mul(1024 * 1024);
    (budget / per_slot).clamp(2, MAX_DIIS_DEPTH)
}

pub(crate) fn slot_bytes(nao: usize, matrices_per_slot: usize) -> usize {
    matrices_per_slot
        .saturating_mul(nao)
        .saturating_mul(nao)
        .saturating_mul(std::mem::size_of::<f64>())
}

/// Whether the minimum useful history (depth 2) fits the budget.
pub(crate) fn diis_depth_fits(nao: usize, matrices_per_slot: usize, budget_mb: usize) -> bool {
    budget_mb == 0 || 2 * slot_bytes(nao, matrices_per_slot) <= budget_mb * 1024 * 1024
}

/// DIIS error `[F,P] = FP  PF`.
///
/// `F` and `P` are both symmetric, so `PF = (FP)T` and the second `O(nao3)`
/// matrix product is redundant: form `FP` once and antisymmetrize it in place.
pub(crate) fn commutator(f: &Matrix, p: &Matrix) -> Matrix {
    debug_assert_eq!(f.rows, p.rows);
    let mut e = f.matmul(p);
    let n = e.rows;
    let data = e.as_mut_slice();
    for i in 0..n {
        data[i * n + i] = 0.0;
        for j in 0..i {
            let upper = data[j * n + i];
            let lower = data[i * n + j];
            data[i * n + j] = lower - upper;
            data[j * n + i] = upper - lower;
        }
    }
    e
}

/// Rolling history for the SCF convergence accelerators.
///
/// The Gram matrices the extrapolators need, `b[i][j] = <E_i,E_j>` for CDIIS
/// and `m[i][j] = <D_i,F_j>` for A-DIIS, are maintained **incrementally**:
/// only the new row and column are evaluated each iteration (`2k - 1` dot
/// products instead of `2k^2`), and the A-DIIS difference matrices `D_i - D_n`
/// and `F_j - F_n` are never materialized. The previous formulation allocated
/// `2k` full `nao x nao` temporaries per iteration, about 416 MiB of allocator churn
/// per iteration at `nao = 1800`, which dominated the SCF wall time.
pub(crate) struct AccelHistory {
    pub(crate) depth: usize,
    pub(crate) keep_densities: bool,
    pub(crate) focks: Vec<Matrix>,
    pub(crate) errors: Vec<Matrix>,
    pub(crate) densities: Vec<Matrix>,
    /// `E_i, E_j`, same order as `errors`.
    pub(crate) b: Vec<Vec<f64>>,
    /// `D_i, F_j`, same order as `densities`/`focks`.
    pub(crate) m: Vec<Vec<f64>>,
}

impl AccelHistory {
    pub(crate) fn new(depth: usize, keep_densities: bool) -> Self {
        Self {
            depth,
            keep_densities,
            focks: Vec::with_capacity(depth),
            errors: Vec::with_capacity(depth),
            densities: Vec::with_capacity(depth),
            b: Vec::with_capacity(depth),
            m: Vec::with_capacity(depth),
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.focks.len()
    }

    /// Append one iterate, updating the Gram matrices with only the new
    /// row/column, then drop the oldest entry if the depth is exceeded.
    /// `density` is required exactly when the history keeps densities (A-DIIS).
    /// Returns `E_new2`, which the Gram update computes anyway.
    pub(crate) fn push(&mut self, fock: Matrix, error: Matrix, density: Option<Matrix>) -> f64 {
        let mut b_row: Vec<f64> = self
            .errors
            .iter()
            .map(|e| e.frobenius_dot(&error))
            .collect();
        let error_norm_sq = error.frobenius_dot(&error);
        b_row.push(error_norm_sq);
        for (i, row) in self.b.iter_mut().enumerate() {
            row.push(b_row[i]);
        }
        self.b.push(b_row);

        if let Some(density) = density.filter(|_| self.keep_densities) {
            // New row `D_new, F_j` and new column `D_i, F_new`.
            let mut m_row: Vec<f64> = self
                .focks
                .iter()
                .map(|f| density.frobenius_dot(f))
                .collect();
            m_row.push(density.frobenius_dot(&fock));
            for (row, d) in self.m.iter_mut().zip(&self.densities) {
                row.push(d.frobenius_dot(&fock));
            }
            self.m.push(m_row);
            self.densities.push(density);
        }
        self.focks.push(fock);
        self.errors.push(error);

        while self.focks.len() > self.depth {
            self.focks.remove(0);
            self.errors.remove(0);
            self.b.remove(0);
            for row in &mut self.b {
                row.remove(0);
            }
            if self.keep_densities {
                self.densities.remove(0);
                self.m.remove(0);
                for row in &mut self.m {
                    row.remove(0);
                }
            }
        }
        error_norm_sq
    }

    /// Pulay CDIIS: extrapolated Fock from the stored `E_i,E_j` Gram matrix.
    pub(crate) fn cdiis(&self) -> Option<Matrix> {
        let coeffs = diis_coeffs_from_gram(&self.b)?;
        Some(combine(&self.focks, &coeffs))
    }

    /// A-DIIS (Hu & Yang 2010) from the stored `D_i,F_j` Gram matrix.
    pub(crate) fn adiis(&self) -> Option<Matrix> {
        let n = self.len();
        if n < 2 || !self.keep_densities {
            return None;
        }
        let last = n - 1;
        let mnn = self.m[last][last];
        // d_i = D_i  D_n, F_n and s_ij = D_i  D_n, F_j  F_n, expanded
        // from the stored inner products (no difference matrices formed).
        let d: Vec<f64> = (0..n).map(|i| self.m[i][last] - mnn).collect();
        let s: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| self.m[i][j] - self.m[i][last] - self.m[last][j] + mnn)
                    .collect()
            })
            .collect();
        let c = solve_adiis_simplex(&d, &s);
        Some(combine(&self.focks, &c))
    }
}

/// Linear combination ` c_i F_i` of the history, chunked over rayon: the
/// history is up to 8 matrices of `nao2` doubles, so this is a pure
/// memory-bandwidth pass that is worth splitting for large bases.
pub(crate) fn combine(fs: &[Matrix], coeffs: &[f64]) -> Matrix {
    use rayon::prelude::*;
    let (r, c) = (fs[0].rows, fs[0].cols);
    let mut out = Matrix::zeros(r, c);
    let n = r * c;
    if n < (1 << 18) {
        for (i, f) in fs.iter().enumerate() {
            let ci = coeffs[i];
            for (o, v) in out.as_mut_slice().iter_mut().zip(f.as_slice()) {
                *o += ci * v;
            }
        }
        return out;
    }
    const CHUNK: usize = 1 << 16;
    out.as_mut_slice()
        .par_chunks_mut(CHUNK)
        .enumerate()
        .for_each(|(k, dst)| {
            let lo = k * CHUNK;
            for (i, f) in fs.iter().enumerate() {
                let ci = coeffs[i];
                let src = &f.as_slice()[lo..lo + dst.len()];
                for (o, v) in dst.iter_mut().zip(src) {
                    *o += ci * v;
                }
            }
        });
    out
}

/// Solve the Pulay DIIS coefficient system from a precomputed `E_i,E_j` Gram matrix.
pub(crate) fn diis_coeffs_from_gram(gram: &[Vec<f64>]) -> Option<Vec<f64>> {
    let n = gram.len();
    if n < 2 {
        return None;
    }
    let dim = n + 1;
    let mut b = Matrix::zeros(dim, dim);
    for i in 0..n {
        for j in 0..n {
            b[(i, j)] = gram[i][j];
        }
        b[(i, n)] = -1.0;
        b[(n, i)] = -1.0;
    }
    let mut rhs = vec![0.0; dim];
    rhs[n] = -1.0;
    // The DIIS matrix (a small bordered saddle-point system) becomes singular near
    // convergence when the error vectors turn linearly dependent. A pivot-guarded
    // Gaussian elimination returns `None` there so the caller falls back to the plain
    // Fock  faer's LU instead returns a degenerate solution that derails DIIS. The heavy
    // O(n^3) eigendecomposition still uses faer; only this tiny solve is bespoke.
    solve_bordered_small(&b, &rhs)
}

/// Gaussian elimination with partial pivoting for the small DIIS system; returns `None`
/// if the matrix is (near-)singular.
pub(crate) fn solve_bordered_small(a: &Matrix, b: &[f64]) -> Option<Vec<f64>> {
    let n = a.rows;
    let mut m = a.clone();
    let mut rhs = b.to_vec();
    for col in 0..n {
        let mut pivot = col;
        let mut best = m[(col, col)].abs();
        for row in (col + 1)..n {
            let v = m[(row, col)].abs();
            if v > best {
                best = v;
                pivot = row;
            }
        }
        if best < 1.0e-12 {
            return None;
        }
        if pivot != col {
            for j in 0..n {
                let t = m[(col, j)];
                m[(col, j)] = m[(pivot, j)];
                m[(pivot, j)] = t;
            }
            rhs.swap(col, pivot);
        }
        for row in (col + 1)..n {
            let factor = m[(row, col)] / m[(col, col)];
            if factor == 0.0 {
                continue;
            }
            for j in col..n {
                let v = m[(col, j)];
                m[(row, j)] -= factor * v;
            }
            rhs[row] -= factor * rhs[col];
        }
    }
    let mut x = vec![0.0; n];
    for col in (0..n).rev() {
        let mut sum = rhs[col];
        for j in (col + 1)..n {
            sum -= m[(col, j)] * x[j];
        }
        x[col] = sum / m[(col, col)];
    }
    Some(x)
}

/// Projected-gradient minimization of the A-DIIS quadratic
/// `f(c) = 2  c_i D_iD_n|F_n +  c_i c_j D_iD_n|F_jF_n`
/// (Hu & Yang 2010) on the probability simplex `{c  0, c = 1}`. The
/// nonnegative weights prevent the runaway extrapolation plain DIIS can produce
/// far from convergence.
pub(crate) fn solve_adiis_simplex(d: &[f64], s: &[Vec<f64>]) -> Vec<f64> {
    let n = d.len();
    // Lipschitz estimate for the step size from (S + ST).
    let mut l: f64 = 1.0e-12;
    for (i, s_i) in s.iter().enumerate().take(n) {
        let mut row = 0.0;
        for (j, &s_ij) in s_i.iter().enumerate().take(n) {
            row += (s_ij + s[j][i]).abs();
        }
        l = l.max(row);
    }
    let lr = 1.0 / l;
    // Start from the latest point (all weight on the newest Fock/density).
    let mut c = vec![0.0; n];
    c[n - 1] = 1.0;
    for _ in 0..400 {
        // grad_k = 2 d_k + _j (s_kj + s_jk) c_j
        let mut g = vec![0.0; n];
        for k in 0..n {
            let mut acc = 2.0 * d[k];
            for j in 0..n {
                acc += (s[k][j] + s[j][k]) * c[j];
            }
            g[k] = acc;
        }
        let trial: Vec<f64> = (0..n).map(|i| c[i] - lr * g[i]).collect();
        let proj = simplex_project(&trial);
        let mut delta = 0.0;
        for i in 0..n {
            delta += (proj[i] - c[i]).abs();
        }
        c = proj;
        if delta < 1.0e-12 {
            break;
        }
    }
    c
}

/// Euclidean projection of `v` onto the probability simplex `{c  0, c = 1}`.
pub(crate) fn simplex_project(v: &[f64]) -> Vec<f64> {
    let mut u = v.to_vec();
    u.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let mut css = 0.0;
    let mut rho = 0;
    let mut theta = 0.0;
    for (j, &uj) in u.iter().enumerate() {
        css += uj;
        let t = (css - 1.0) / (j as f64 + 1.0);
        if uj - t > 0.0 {
            rho = j + 1;
            theta = t;
        }
    }
    let _ = rho;
    v.iter().map(|&vi| (vi - theta).max(0.0)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diis_depth_respects_the_memory_budget() {
        // Unlimited budget keeps the full depth.
        assert_eq!(diis_depth(1800, 3, 0), MAX_DIIS_DEPTH);
        // 512 MiB / (3 * 1800^2 * 8 B = 74.2 MiB) = 6 slots.
        assert_eq!(diis_depth(1800, 3, 512), 6);
        // A huge basis is floored at 2 rather than 0.
        assert_eq!(diis_depth(8000, 3, 512), 2);
        // ...and at that size a depth-2 A-DIIS history does not fit, so the
        // caller degrades to CDIIS instead of blowing the budget.
        assert!(!diis_depth_fits(8000, 3, 512));
        assert!(diis_depth_fits(1800, 3, 512));
        // A small system is unaffected.
        assert_eq!(diis_depth(60, 3, 512), MAX_DIIS_DEPTH);
    }

    #[test]
    fn a_uhf_slot_costs_twice_a_closed_shell_slot() {
        // Both channels live in one history, so a slot is twice the size and
        // the budget buys half the depth. Getting this wrong would silently
        // double the accelerator's memory beyond the configured budget.
        let rhf = Accelerator::new(1800, ScfAccelerator::AdiisCdiis, 0.1, 512);
        let uhf = Accelerator::new_uhf(1800, ScfAccelerator::AdiisCdiis, 0.1, 512);
        assert_eq!(rhf.depth(), 6);
        assert_eq!(uhf.depth(), 3);
    }

    #[test]
    fn stacking_spin_channels_round_trips() {
        let nao = 5;
        let mut a = Matrix::zeros(nao, nao);
        let mut b = Matrix::zeros(nao, nao);
        for i in 0..nao {
            for j in 0..nao {
                a[(i, j)] = (i * nao + j) as f64;
                b[(i, j)] = -((i * nao + j) as f64) - 0.5;
            }
        }
        let stacked = stack_channels(&a, &b);
        assert_eq!(stacked.rows, 2 * nao);
        // The stacked Frobenius product is the UHF pairing Tr[P_a F_a] + Tr[P_b F_b],
        // which is what makes A-DIIS correct on the stack with no new functional.
        let paired = stack_channels(&b, &a);
        assert!(
            (stacked.frobenius_dot(&paired) - (a.frobenius_dot(&b) + b.frobenius_dot(&a))).abs()
                < 1.0e-12
        );
        let (a2, b2) = split_channels(&stacked, nao);
        assert!(a2.rms_difference(&a) < 1.0e-15);
        assert!(b2.rms_difference(&b) < 1.0e-15);
    }

    #[test]
    fn an_inactive_accelerator_returns_the_fock_untouched() {
        let nao = 3;
        let mut f = Matrix::zeros(nao, nao);
        f[(0, 0)] = 1.0;
        let p = Matrix::zeros(nao, nao);
        let mut accel = Accelerator::new(nao, ScfAccelerator::None, 0.1, 512);
        assert!(!accel.is_active());
        let out = accel.step(f.clone(), &p);
        assert!(out.rms_difference(&f) < 1.0e-15);
        let mut uhf = Accelerator::new_uhf(nao, ScfAccelerator::None, 0.1, 512);
        let (fa, fb) = uhf.step_uhf(f.clone(), f.clone(), &p, &p);
        assert!(fa.rms_difference(&f) < 1.0e-15);
        assert!(fb.rms_difference(&f) < 1.0e-15);
    }

    #[test]
    fn sad_fills_s_before_p_and_keeps_the_electron_count() {
        // Carbon: 4 valence electrons over 4 AOs. The old guess put 1.0 in each;
        // SAD puts 2 in s and 2/3 in each p.
        let p = sad_density_sp(4, &[(0, 4, 4.0)]);
        assert!((p[(0, 0)] - 2.0).abs() < 1.0e-12);
        for q in 1..4 {
            assert!((p[(q, q)] - 2.0 / 3.0).abs() < 1.0e-12);
        }
        let trace: f64 = (0..4).map(|i| p[(i, i)]).sum();
        assert!((trace - 4.0).abs() < 1.0e-12);
    }

    fn trace(m: &Matrix) -> f64 {
        (0..m.rows).map(|i| m[(i, i)]).sum()
    }

    #[test]
    fn splitting_by_spin_gives_each_channel_its_electron_count() {
        // A doublet: 5 alpha, 4 beta. Each channel's trace is its own count, and
        // the two differ, which is what breaks spin symmetry in the guess.
        let total = sad_density_sp(4, &[(0, 4, 4.0)]);
        let mut scaled = total.clone();
        for v in scaled.as_mut_slice() {
            *v *= 9.0 / 4.0;
        }
        let (pa, pb) = split_by_spin(&scaled, 5, 4);
        assert!((trace(&pa) - 5.0).abs() < 1.0e-12);
        assert!((trace(&pb) - 4.0).abs() < 1.0e-12);
        assert!(pa.rms_difference(&pb) > 1.0e-3);
    }

    #[test]
    fn the_spin_split_reproduces_hunds_rule_for_a_free_atom() {
        // Not a fitted rule: weighting the excess spin by each AO's capacity
        // `min(P, 2-P)` gives the spherically averaged maximum-multiplicity
        // configuration of every free atom. If this ever stops holding, the
        // weight is no longer the physical capacity.
        //
        //   (valence electrons, 2S, expected alpha/beta on s, on each p)
        for (zc, unpaired, s_ab, p_ab) in [
            (4.0, 2.0, (1.0, 1.0), (2.0 / 3.0, 0.0)), // C: s2 px1 py1
            (5.0, 3.0, (1.0, 1.0), (1.0, 0.0)),       // N: s2 px1 py1 pz1
            (6.0, 2.0, (1.0, 1.0), (1.0, 1.0 / 3.0)), // O: s2 px2 py1 pz1
            (7.0, 1.0, (1.0, 1.0), (1.0, 2.0 / 3.0)), // F: s2 px2 py2 pz1
        ] {
            let total = sad_density_sp(4, &[(0, 4, zc)]);
            let n_beta = ((zc - unpaired) / 2.0) as usize;
            let n_alpha = n_beta + unpaired as usize;
            let (pa, pb) = split_by_spin(&total, n_alpha, n_beta);
            assert!(
                (pa[(0, 0)] - s_ab.0).abs() < 1.0e-12 && (pb[(0, 0)] - s_ab.1).abs() < 1.0e-12,
                "Z_v={zc}: s channel is ({}, {}), wanted {s_ab:?}",
                pa[(0, 0)],
                pb[(0, 0)]
            );
            for q in 1..4 {
                assert!(
                    (pa[(q, q)] - p_ab.0).abs() < 1.0e-12 && (pb[(q, q)] - p_ab.1).abs() < 1.0e-12,
                    "Z_v={zc}: p channel is ({}, {}), wanted {p_ab:?}",
                    pa[(q, q)],
                    pb[(q, q)]
                );
            }
            assert!((trace(&pa) - n_alpha as f64).abs() < 1.0e-12, "Z_v={zc}");
            assert!((trace(&pb) - n_beta as f64).abs() < 1.0e-12, "Z_v={zc}");
        }
    }

    #[test]
    fn no_excess_spin_reaches_a_doubly_occupied_orbital() {
        // The defect this rule replaced: weighting by `P` put the largest share
        // of the unpaired spin on the *most* occupied AO. Carbon's 2s is full in
        // the guess, so it must carry none.
        let total = sad_density_sp(4, &[(0, 4, 4.0)]);
        let (pa, pb) = split_by_spin(&total, 3, 1);
        assert!(
            (pa[(0, 0)] - pb[(0, 0)]).abs() < 1.0e-12,
            "the doubly occupied s picked up spin: {} vs {}",
            pa[(0, 0)],
            pb[(0, 0)]
        );
        assert!(pa[(1, 1)] - pb[(1, 1)] > 0.1, "the p shell carries none");
    }

    #[test]
    fn the_spin_split_never_leaves_a_channel_outside_its_bounds() {
        // Per AO a spin channel holds between 0 and 1 electron. Capacity
        // weighting is supposed to guarantee that without a separate clamp, so
        // push it: more unpaired electrons than the guess can hold.
        let total = sad_density_sp(8, &[(0, 4, 4.0), (4, 4, 1.0)]);
        for (n_alpha, n_beta) in [(5, 0), (4, 1), (3, 2), (1, 4)] {
            let (pa, pb) = split_by_spin(&total, n_alpha, n_beta);
            for i in 0..8 {
                for (label, m) in [("alpha", &pa), ("beta", &pb)] {
                    assert!(
                        m[(i, i)] >= -1.0e-12 && m[(i, i)] <= 1.0 + 1.0e-12,
                        "{label}[{i}] = {} outside [0, 1] for ({n_alpha}, {n_beta})",
                        m[(i, i)]
                    );
                }
            }
            // The total is preserved whatever the split does.
            assert!((trace(&pa) + trace(&pb) - trace(&total)).abs() < 1.0e-12);
        }
    }

    #[test]
    fn a_closed_shell_split_is_exactly_half_and_half() {
        let total = sad_density_sp(8, &[(0, 4, 6.0), (4, 4, 4.0)]);
        let (pa, pb) = split_by_spin(&total, 5, 5);
        assert!(pa.rms_difference(&pb) < 1.0e-15);
        for i in 0..8 {
            assert!((pa[(i, i)] - 0.5 * total[(i, i)]).abs() < 1.0e-15);
        }
    }
}
