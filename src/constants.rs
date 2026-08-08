// SPDX-License-Identifier: GPL-3.0-or-later

//! Physical constants and unit conversions used by the NDDO engines.
//!
//! OpenMOPAC v23.2.5 uses the 2018 CODATA constants in `conref_C.F90`.
//! xndo-rs freezes the same values so MNDO and MNDO/d calculations use the
//! exact upstream model convention.
//!
//! PROVENANCE: openmopac/mopac v23.2.5 `src/conref_C.F90` (Apache-2.0).
//! See THIRD_PARTY_NOTICES.md.
//!
//! Internal unit policy: energies in **eV**,
//! distances in **Bohr** throughout the computational core. The Rust and
//! Python-native API boundaries convert to Hartree/Bohr (atomic units); the
//! ASE calculator converts to eV/A.

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
}
