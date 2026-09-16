// SPDX-License-Identifier: GPL-3.0-or-later

//! Physical constants and unit conversions used by the NDDO engines.
//!
//! OpenMOPAC v23.2.5 uses the 2018 CODATA constants in `conref_C.F90`.
//! xndo-rs freezes the same values so MNDO and MNDO/d calculations use the
//! exact upstream model convention.
//!
//!
//! Internal unit policy: energies in **eV**,
//! distances in **Bohr** throughout the computational core. The Rust and
//! Python-native API boundaries convert to Hartree/Bohr (atomic units); the
//! ASE calculator converts to eV/A.
//!
//! PROVENANCE: derived from MOPAC (Molecular Orbital PACkage) v23.2.5,
//! Copyright 2021 Virginia Polytechnic Institute and State University,
//! licensed under the Apache License, Version 2.0.
//! UPSTREAM: src/conref_C.F90.
//! MODIFIED for xndo-rs v0.3.0 on 2026-09-14:
//! transcribed the two constant sets into `ModelConstants` and added the
//! MolDS and INDO/S sets, which upstream keeps elsewhere.
//! Retained notices: NOTICE; per-file record: THIRD_PARTY_NOTICES.md.

/// One Hartree in eV (2018 CODATA, MOPAC `fpc(4)` = `fpcref(1,4)`).
pub const HARTREE_TO_EV: f64 = 27.211386245988;

/// Bohr radius in Angstrom (2018 CODATA, MOPAC `fpc(3)` = `fpcref(1,3)`).
pub const BOHR_TO_ANGSTROM: f64 = 0.529177210903;

/// One eV in kcal/mol (2018 CODATA derived, MOPAC `fpc(9)` = `fpcref(1,9)`).
pub const EV_TO_KCAL: f64 = 23.060_547_830_619_03;

/// `a0 * ev`  one (elementary charge)2/A in eV (MOPAC `fpc(2)` = `fpcref(1,2)`).
/// MOPAC's two-electron integral kernels use this product directly.
pub const A0_TIMES_EV: f64 = 14.399645478456;

pub const EV_TO_HARTREE: f64 = 1.0 / HARTREE_TO_EV;
pub const ANGSTROM_TO_BOHR: f64 = 1.0 / BOHR_TO_ANGSTROM;
pub const KCAL_TO_EV: f64 = 1.0 / EV_TO_KCAL;
pub const HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM: f64 = HARTREE_TO_EV / BOHR_TO_ANGSTROM;

/// Dipole conversion: one electronBohr in Debye.
/// MOPAC (`dipole.F90`) computes point-charge dipoles as `4.803 * charge * A`;
/// the atomic-unit equivalent used here is `ea0  Debye`.
pub const AU_DIPOLE_TO_DEBYE: f64 = 2.541746473;

/// The set of fundamental physical constants a parameterisation was fitted with.
///
/// MOPAC carries two sets (`conref_C.F90`): the 2018 CODATA values, and a
/// historical set from before the 2019 SI revision. `readmo.F90:411` selects the
/// historical set for the `OLDFPC` **and `MNDOD`** keywords -- so MOPAC evaluates
/// MNDO/d with a Bohr radius of 0.529167 A and a Hartree of 27.21 eV, while every
/// other method uses CODATA. That is a real part of the MNDO/d model, not a
/// rounding artefact: Thiel and Voityuk fitted the parameters against those
/// constants, and the relative offsets (5.1e-5 in the Hartree, 1.9e-5 in the Bohr
/// radius) accumulate over pairs. It is why MOPAC's own `MNDO` and `MNDOD` differ
/// by 2.25e-3 kcal/mol on H2 -- a molecule whose parameters are identical in the
/// two tables. See `tests/data/ORACLE_NOTES.md` item 13.
///
/// # How a non-CODATA set is applied
///
/// The crate's internal length unit is the **CODATA** Bohr everywhere, because the
/// geometry is converted once, before any method is chosen. A parameterisation
/// fitted with a different Bohr radius is therefore folded into the parameters
/// rather than into the coordinates. Writing `k = a0_CODATA / a0_model`:
///
/// * every Slater exponent (an inverse length) is scaled by `k`;
/// * every energy prefactor uses `ev_model / k` ([`Self::effective_hartree_ev`]).
///
/// Both substitutions are exact. For the two-centre integrals,
/// `ev/sqrt(r_model^2 + rho_model^2)` in the model's own Bohr equals
/// `(ev/k)/sqrt(r_CODATA^2 + rho_CODATA^2)`, because `rho = 0.5 ev / G` scales the
/// same way `r` does. For the one-centre Slater-Condon integrals, which are
/// linear in the exponent, `rsc(zeta) * ev_model == rsc(k zeta) * ev_model / k`.
/// Distances in Angstrom -- the core-core exponents `alpha` are per Angstrom --
/// are physical and unaffected.
///
/// The CODATA set has `k == 1` exactly, so MNDO is bit-identical to before.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModelConstants {
    /// One Hartree in eV, as the parameterisation used it.
    pub hartree_to_ev: f64,
    /// Bohr radius in Angstrom, as the parameterisation used it.
    pub bohr_to_angstrom: f64,
    /// One eV in kcal/mol, as the parameterisation used it.
    pub ev_to_kcal: f64,
    /// One electron-Angstrom in Debye, as the parameterisation used it.
    ///
    /// MOPAC computes this rather than tabulating it: `dipole.F90:153` forms
    /// `factdm = fpc(8) * fpc(1) * 1e-10` from the speed of light and the
    /// elementary charge, so it inherits whichever constant set is active. The
    /// historical speed of light is 5.5e-5 low, which is the whole of the
    /// remaining MNDO/d dipole disagreement -- 3.8e-4 D on NaH's 6 D dipole.
    pub debye_per_e_angstrom: f64,
}

impl ModelConstants {
    /// 2018 CODATA: MOPAC `fpcref(1,*)`. Every method but MNDO/d.
    pub const CODATA_2018: Self = Self {
        hartree_to_ev: HARTREE_TO_EV,
        bohr_to_angstrom: BOHR_TO_ANGSTROM,
        ev_to_kcal: EV_TO_KCAL,
        // 2.99792458e10 cm/s * 1.602176634e-19 C * 1e-10.
        debye_per_e_angstrom: 4.803_204_712_570_264,
    };

    /// Pre-2019 constants: MOPAC `fpcref(2,*)`, selected by `OLDFPC` and `MNDOD`.
    pub const HISTORICAL: Self = Self {
        hartree_to_ev: 27.21,
        bohr_to_angstrom: 0.529167,
        ev_to_kcal: 23.061,
        // 2.99776e10 * 1.60217733e-19 * 1e-10.
        debye_per_e_angstrom: 4.802_943_112_780_8,
    };

    /// The constants MOPAC's INDO module carries privately
    /// (`src/INDO/reimers_C.F90:93`: `au2ev, au2ang, au2cm, debye`), which are
    /// neither of the two sets `conref_C.F90` offers. The Debye conversion there
    /// is the historical one rounded to six figures; the Hartree and the Bohr
    /// radius are near-CODATA and move the energies by about 3e-6 eV, which is
    /// below what MOPAC prints.
    pub const INDO_S: Self = Self {
        hartree_to_ev: 27.2114,
        bohr_to_angstrom: 0.529177,
        ev_to_kcal: EV_TO_KCAL,
        debye_per_e_angstrom: 4.80294,
    };

    /// The constants MolDS 0.3.1 uses (`src/base/Parameters.cpp:45-57`).
    ///
    /// MolDS works in atomic units internally and converts on output with
    /// `eV2AU = 0.03674903`, which is a Hartree of 27.211602592 eV -- **7.95e-6
    /// above** CODATA. Its `angstrom2AU = 1/0.5291772` and
    /// `debye2AU = 0.393430191` are within 2e-7 of CODATA and contribute
    /// nothing measurable, but they are carried here so the set is the whole
    /// of what MolDS uses rather than the part that happened to matter.
    ///
    /// The Hartree offset is small but it multiplies the core repulsion, which
    /// is hundreds of eV: on ethane it is 5.75e-3 eV, which was the largest
    /// disagreement in the CNDO/2 suite.
    pub const MOLDS: Self = Self {
        hartree_to_ev: 27.211_602_591_961_4,
        bohr_to_angstrom: 0.5291772,
        ev_to_kcal: EV_TO_KCAL,
        // 1 / 0.393430191, per Angstrom.
        debye_per_e_angstrom: 4.803_205_770_146_6,
    };

    /// `a0_CODATA / a0_model`: the factor that takes an inverse length from the
    /// model's Bohr to the crate's internal Bohr. Exactly 1 for CODATA.
    pub fn length_scale(&self) -> f64 {
        BOHR_TO_ANGSTROM / self.bohr_to_angstrom
    }

    /// The energy prefactor to use in place of [`HARTREE_TO_EV`] once lengths are
    /// expressed in the internal (CODATA) Bohr. Exactly [`HARTREE_TO_EV`] for CODATA.
    pub fn effective_hartree_ev(&self) -> f64 {
        self.hartree_to_ev / self.length_scale()
    }

    /// Dipole prefactor for an internal (CODATA) Bohr dipole, in place of
    /// [`AU_DIPOLE_TO_DEBYE`]. Exactly that constant for CODATA, to 2e-8.
    pub fn debye_per_e_bohr(&self) -> f64 {
        self.debye_per_e_angstrom * BOHR_TO_ANGSTROM
    }
}

impl Default for ModelConstants {
    fn default() -> Self {
        Self::CODATA_2018
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pin the MOPAC 2018-CODATA model constants (extraction checkpoint).
    #[test]
    fn openmopac_codata_2018_constants() {
        assert_eq!(HARTREE_TO_EV, 27.211386245988);
        assert_eq!(BOHR_TO_ANGSTROM, 0.529177210903);
        assert_eq!(EV_TO_KCAL, 23.060_547_830_619_03);
        // fpc(2) is the product a0*ev, tabulated separately by MOPAC.
        assert!((A0_TIMES_EV - BOHR_TO_ANGSTROM * HARTREE_TO_EV).abs() < 1e-9);
    }

    /// Pin MOPAC's historical set, `fpcref(2,*)`, which `MNDOD` selects.
    #[test]
    fn openmopac_historical_constants() {
        let h = ModelConstants::HISTORICAL;
        assert_eq!(h.hartree_to_ev, 27.21);
        assert_eq!(h.bohr_to_angstrom, 0.529167);
        assert_eq!(h.ev_to_kcal, 23.061);
    }

    /// Pin MolDS's set, which is neither of MOPAC's.
    #[test]
    fn molds_constants() {
        let m = ModelConstants::MOLDS;
        // MolDS stores the inverse; this is the Hartree it implies.
        assert!((m.hartree_to_ev - 1.0 / 0.03674903).abs() < 1e-9);
        assert_eq!(m.bohr_to_angstrom, 0.5291772);
        assert!((m.debye_per_e_angstrom * 0.5291772 - 1.0 / 0.393430191).abs() < 1e-7);
        let rel = (m.hartree_to_ev - HARTREE_TO_EV) / HARTREE_TO_EV;
        assert!(
            (7.9e-6..8.0e-6).contains(&rel),
            "MolDS Hartree offset {rel:e}"
        );
    }

    /// MOPAC forms the Debye conversion from the speed of light and the
    /// elementary charge rather than tabulating it, so it must be reproduced
    /// that way or it will drift from whichever set is active.
    #[test]
    fn the_debye_conversion_is_c_times_e() {
        let codata = 2.997_924_58e10 * 1.602_176_634 * 1.0e-10;
        assert!((ModelConstants::CODATA_2018.debye_per_e_angstrom - codata).abs() < 1e-12);
        let historical = 2.997_76e10 * 1.602_177_33 * 1.0e-10;
        assert!((ModelConstants::HISTORICAL.debye_per_e_angstrom - historical).abs() < 1e-12);
        // The historical speed of light is what makes this 5.5e-5 low.
        let rel = (ModelConstants::HISTORICAL.debye_per_e_angstrom - codata) / codata;
        assert!((-5.6e-5..-5.4e-5).contains(&rel), "Debye offset {rel:e}");
    }

    /// The existing `AU_DIPOLE_TO_DEBYE` and the CODATA set must agree, or the
    /// MNDO dipoles would move when they are routed through the new path.
    #[test]
    fn the_codata_dipole_prefactor_matches_the_standalone_constant() {
        let via_set = ModelConstants::CODATA_2018.debye_per_e_bohr();
        assert!(
            (via_set - AU_DIPOLE_TO_DEBYE).abs() < 1e-7,
            "{via_set} vs {AU_DIPOLE_TO_DEBYE}"
        );
    }

    /// The CODATA set must be the exact identity, so every method but MNDO/d is
    /// bit-identical to a build that never knew about model constants.
    #[test]
    fn the_codata_set_is_the_identity_transformation() {
        let c = ModelConstants::CODATA_2018;
        assert_eq!(c.length_scale(), 1.0);
        assert_eq!(c.effective_hartree_ev(), HARTREE_TO_EV);
    }

    /// The historical Hartree is 5.09e-5 low and the historical Bohr radius
    /// 1.93e-5 low. The prefactor that actually reaches the integrals carries
    /// both, so it is 7.02e-5 low -- the two offsets add rather than cancel --
    /// and that, accumulated over pairs, is the whole of ORACLE_NOTES item 13.
    #[test]
    fn the_historical_set_shifts_the_energy_prefactor_by_seven_parts_in_a_hundred_thousand() {
        let h = ModelConstants::HISTORICAL;
        let rel = (h.effective_hartree_ev() - HARTREE_TO_EV) / HARTREE_TO_EV;
        assert!(
            (-7.1e-5..-6.9e-5).contains(&rel),
            "effective Hartree offset {rel:e}"
        );
        let d_ev = (h.hartree_to_ev - HARTREE_TO_EV) / HARTREE_TO_EV;
        let d_a0 = (h.bohr_to_angstrom - BOHR_TO_ANGSTROM) / BOHR_TO_ANGSTROM;
        assert!((rel - (d_ev + d_a0)).abs() < 1.0e-8);
    }
}
