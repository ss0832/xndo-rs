// SPDX-License-Identifier: GPL-3.0-or-later

//! Orbital energies, and the frontier quantities read off them.
//!
//! Each engine returns its own result type, and before this module each of them
//! reported the orbital information differently: the ZDO engines carried the
//! whole spectrum but no frontier pair, while the NDDO engine carried a frontier
//! pair but discarded the beta spectrum entirely. [`OrbitalEnergies`] is the one
//! shape all four are converted to, so `homo_ev()` means the same thing for
//! every method.
//!
//! **What "the HOMO" means for an unrestricted reference.** Alpha and beta
//! eigenvalues come from two *different* Fock operators, so they are not levels
//! of one orbital diagram and pairing them into a single diagram is a fiction.
//! What is well defined is the set of occupied spin orbitals, so [`Self::homo_ev`]
//! returns the highest occupied over both channels and [`Self::lumo_ev`] the
//! lowest unoccupied over both. For a doublet radical the beta LUMO is the
//! partner of the singly occupied orbital and normally lies *below* the alpha
//! LUMO, which is why taking the alpha channel alone (what this crate did before
//! v0.3.0) reports a LUMO that is not the lowest unoccupied orbital. The
//! per-channel values are public too, so a caller who wants one spin's diagram
//! never has to go through the combined pair.

/// The converged reference determinant's orbital energies, in eV.
///
/// `alpha_ev` is ascending, as it comes off the eigensolver. `beta_ev` is
/// `None` for a restricted reference — not a copy of `alpha_ev` — so that
/// "restricted" stays distinguishable from "unrestricted that happens to be
/// spin-symmetric".
#[derive(Clone, Debug, PartialEq)]
pub struct OrbitalEnergies {
    /// Alpha (or, for a restricted reference, the doubly occupied) spectrum.
    pub alpha_ev: Vec<f64>,
    /// Beta spectrum; `None` for a restricted reference.
    pub beta_ev: Option<Vec<f64>>,
    /// Number of occupied alpha spin orbitals.
    pub n_alpha: usize,
    /// Number of occupied beta spin orbitals.
    pub n_beta: usize,
}

impl OrbitalEnergies {
    /// A restricted reference: one spectrum, `n_occ` doubly occupied orbitals.
    pub fn restricted(alpha_ev: Vec<f64>, n_occ: usize) -> Self {
        Self {
            alpha_ev,
            beta_ev: None,
            n_alpha: n_occ,
            n_beta: n_occ,
        }
    }

    /// An unrestricted reference: two spectra with their own occupation counts.
    pub fn unrestricted(
        alpha_ev: Vec<f64>,
        n_alpha: usize,
        beta_ev: Vec<f64>,
        n_beta: usize,
    ) -> Self {
        Self {
            alpha_ev,
            beta_ev: Some(beta_ev),
            n_alpha,
            n_beta,
        }
    }

    /// Build from the optional beta spectrum the engines carry: `Some` selects
    /// the unrestricted form, `None` the restricted one.
    pub fn new(
        alpha_ev: Vec<f64>,
        n_alpha: usize,
        beta_ev: Option<Vec<f64>>,
        n_beta: usize,
    ) -> Self {
        match beta_ev {
            Some(beta) => Self::unrestricted(alpha_ev, n_alpha, beta, n_beta),
            None => Self {
                alpha_ev,
                beta_ev: None,
                n_alpha,
                n_beta,
            },
        }
    }

    /// Whether this came from an unrestricted reference.
    pub fn unrestricted_reference(&self) -> bool {
        self.beta_ev.is_some()
    }

    /// Number of orbitals in each spectrum (the size of the AO basis).
    pub fn n_orbitals(&self) -> usize {
        self.alpha_ev.len()
    }

    /// The beta spectrum as the caller should read it: the alpha spectrum for a
    /// restricted reference, where the two are equal by construction.
    pub fn beta_or_alpha_ev(&self) -> &[f64] {
        self.beta_ev.as_deref().unwrap_or(&self.alpha_ev)
    }

    /// Highest occupied alpha orbital energy; `None` if no alpha electrons.
    pub fn homo_alpha_ev(&self) -> Option<f64> {
        self.n_alpha
            .checked_sub(1)
            .and_then(|i| self.alpha_ev.get(i))
            .copied()
    }

    /// Lowest unoccupied alpha orbital energy; `None` if the alpha shell is full.
    pub fn lumo_alpha_ev(&self) -> Option<f64> {
        self.alpha_ev.get(self.n_alpha).copied()
    }

    /// Highest occupied beta orbital energy; `None` if no beta electrons.
    pub fn homo_beta_ev(&self) -> Option<f64> {
        self.n_beta
            .checked_sub(1)
            .and_then(|i| self.beta_or_alpha_ev().get(i))
            .copied()
    }

    /// Lowest unoccupied beta orbital energy; `None` if the beta shell is full.
    pub fn lumo_beta_ev(&self) -> Option<f64> {
        self.beta_or_alpha_ev().get(self.n_beta).copied()
    }

    /// Highest occupied spin orbital over both channels (the HOMO).
    ///
    /// For a restricted reference this is just the highest doubly occupied
    /// orbital. See the module docs for why the maximum, and not the alpha
    /// channel alone, is the right reading for an unrestricted one.
    pub fn homo_ev(&self) -> Option<f64> {
        max_opt(self.homo_alpha_ev(), self.homo_beta_ev())
    }

    /// Lowest unoccupied spin orbital over both channels (the LUMO).
    pub fn lumo_ev(&self) -> Option<f64> {
        min_opt(self.lumo_alpha_ev(), self.lumo_beta_ev())
    }

    /// HOMO-LUMO gap in eV; `None` unless both frontier orbitals exist.
    ///
    /// Never negative for a converged aufbau solution, but it is *not* clamped:
    /// a negative gap means the SCF settled on a non-aufbau occupation, and
    /// hiding that would hide a real result.
    pub fn gap_ev(&self) -> Option<f64> {
        Some(self.lumo_ev()? - self.homo_ev()?)
    }

    /// Per-channel gap, alpha then beta. Each is `None` when that channel is
    /// missing a frontier orbital.
    pub fn gap_alpha_ev(&self) -> Option<f64> {
        Some(self.lumo_alpha_ev()? - self.homo_alpha_ev()?)
    }

    /// See [`Self::gap_alpha_ev`].
    pub fn gap_beta_ev(&self) -> Option<f64> {
        Some(self.lumo_beta_ev()? - self.homo_beta_ev()?)
    }

    /// Orbital occupation numbers for the alpha channel (1.0 or 0.0 when
    /// unrestricted, 2.0 or 0.0 when restricted), in spectrum order.
    ///
    /// Molden's `Occup=` lines want exactly this.
    pub fn occupations_alpha(&self) -> Vec<f64> {
        let full = if self.unrestricted_reference() {
            1.0
        } else {
            2.0
        };
        (0..self.alpha_ev.len())
            .map(|i| if i < self.n_alpha { full } else { 0.0 })
            .collect()
    }

    /// Beta occupations; `None` for a restricted reference, whose occupations
    /// are already carried by [`Self::occupations_alpha`].
    pub fn occupations_beta(&self) -> Option<Vec<f64>> {
        let beta = self.beta_ev.as_ref()?;
        Some(
            (0..beta.len())
                .map(|i| if i < self.n_beta { 1.0 } else { 0.0 })
                .collect(),
        )
    }
}

/// `f64::max` propagates through `None` the wrong way for this use: a missing
/// channel must not win, and two `None`s must stay `None`.
fn max_opt(a: Option<f64>, b: Option<f64>) -> Option<f64> {
    match (a, b) {
        (Some(x), Some(y)) => Some(x.max(y)),
        (x, None) => x,
        (None, y) => y,
    }
}

fn min_opt(a: Option<f64>, b: Option<f64>) -> Option<f64> {
    match (a, b) {
        (Some(x), Some(y)) => Some(x.min(y)),
        (x, None) => x,
        (None, y) => y,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_restricted_reference_reads_its_frontier_pair_off_one_spectrum() {
        let o = OrbitalEnergies::restricted(vec![-30.0, -20.0, -10.0, 2.0, 5.0], 3);
        assert_eq!(o.homo_ev(), Some(-10.0));
        assert_eq!(o.lumo_ev(), Some(2.0));
        assert_eq!(o.gap_ev(), Some(12.0));
        assert!(!o.unrestricted_reference());
        // Both channels report the same thing, because there is only one.
        assert_eq!(o.homo_alpha_ev(), o.homo_beta_ev());
        assert_eq!(o.occupations_alpha(), vec![2.0, 2.0, 2.0, 0.0, 0.0]);
        assert_eq!(o.occupations_beta(), None);
    }

    #[test]
    fn a_doublet_takes_its_lumo_from_the_beta_channel() {
        // The signature of an open-shell doublet: the singly occupied orbital is
        // occupied in alpha and empty in beta, so the beta LUMO sits below the
        // alpha LUMO. Reading only the alpha channel (what the NDDO engine did
        // before v0.3.0) reports 4.0 here, which is not the lowest unoccupied
        // orbital.
        let alpha = vec![-25.0, -12.0, -3.0, 4.0];
        let beta = vec![-24.0, -11.0, -2.5, 4.5];
        let o = OrbitalEnergies::unrestricted(alpha, 3, beta, 2);
        assert_eq!(o.homo_alpha_ev(), Some(-3.0));
        assert_eq!(o.lumo_alpha_ev(), Some(4.0));
        assert_eq!(o.homo_beta_ev(), Some(-11.0));
        assert_eq!(o.lumo_beta_ev(), Some(-2.5));
        assert_eq!(o.homo_ev(), Some(-3.0));
        assert_eq!(o.lumo_ev(), Some(-2.5));
        assert_eq!(o.gap_ev(), Some(0.5));
        assert_eq!(o.occupations_alpha(), vec![1.0, 1.0, 1.0, 0.0]);
        assert_eq!(o.occupations_beta(), Some(vec![1.0, 1.0, 0.0, 0.0]));
    }

    #[test]
    fn a_full_or_empty_shell_reports_none_rather_than_an_invented_number() {
        // Every orbital occupied: there is no LUMO. Indexing past the end (or
        // subtracting one from zero) is the bug this guards.
        let full = OrbitalEnergies::restricted(vec![-9.0, -4.0], 2);
        assert_eq!(full.homo_ev(), Some(-4.0));
        assert_eq!(full.lumo_ev(), None);
        assert_eq!(full.gap_ev(), None);

        let empty = OrbitalEnergies::restricted(vec![-9.0, -4.0], 0);
        assert_eq!(empty.homo_ev(), None);
        assert_eq!(empty.lumo_ev(), Some(-9.0));
        assert_eq!(empty.gap_ev(), None);

        // A bare proton: no orbitals at all, so no frontier pair either.
        let none = OrbitalEnergies::restricted(Vec::new(), 0);
        assert_eq!(none.homo_ev(), None);
        assert_eq!(none.lumo_ev(), None);
    }

    #[test]
    fn a_high_spin_reference_with_an_empty_beta_channel_still_has_a_homo() {
        // A triplet with every beta electron removed, e.g. a bare cation shell:
        // the beta channel contributes no HOMO, and the combined answer must be
        // the alpha one, not `None` and not the beta LUMO.
        let o = OrbitalEnergies::unrestricted(vec![-20.0, -8.0], 2, vec![-19.0, -7.0], 0);
        assert_eq!(o.homo_beta_ev(), None);
        assert_eq!(o.homo_ev(), Some(-8.0));
        assert_eq!(o.lumo_ev(), Some(-19.0));
        assert_eq!(o.gap_alpha_ev(), None);
        assert_eq!(o.gap_beta_ev(), None);
    }

    #[test]
    fn the_gap_is_reported_negative_rather_than_clamped() {
        // A non-aufbau fixed point. Clamping to zero would turn a result worth
        // investigating into a plausible-looking number.
        let o = OrbitalEnergies::restricted(vec![-5.0, -1.0, -3.0], 2);
        assert_eq!(o.homo_ev(), Some(-1.0));
        assert_eq!(o.lumo_ev(), Some(-3.0));
        assert_eq!(o.gap_ev(), Some(-2.0));
    }
}
