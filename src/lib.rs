// SPDX-License-Identifier: GPL-3.0-or-later

#![forbid(unsafe_code)]

//! # xndo-rs
//!
//! `xndo-rs` is a Rust implementation framework for legacy NDO-family
//! semiempirical molecular-orbital methods.  Numerical engines included in
//! this release are CNDO/2, ground-state INDO, MNDO, MNDO/d, and MINDO/3
//! (RHF/UHF with analytic derivatives), and the s/p branch of Zerner INDO/S
//! (ZINDO/S; RHF/UHF derivatives, spin-adapted RHF-CIS, and spin-orbital UCIS
//! properties and state-specific derivatives).
//!
//! Additional historical method names are registered so that their eventual
//! implementations have a stable API: CNDO/1, INDO/1, INDO/2, ZINDO/1,
//! ZINDO/2, MINDO/1, MINDO/2, SINDO1, and MSINDO.  A registered name
//! is not silently treated as another Hamiltonian; if its complete,
//! provenance-audited parameter set is not shipped, execution returns an
//! explicit error. CNDO/3, INDO/3, and ZINDO/3 are retained only as requested
//! compatibility labels and are explicitly marked non-canonical.
//!
//! Internal units of the inherited NDDO engine are eV and Bohr.  Public
//! Python input coordinates are Angstrom. See `docs/methods.md` and
//! `docs/oracles.md` for scope and validation policy.

pub mod basis;
pub mod cndo_indo;
pub mod constants;
pub mod data_tables;
pub mod dual;
pub mod dual2;
pub mod error;
pub mod fock;
pub mod frame;
pub mod gradient;
pub mod hamiltonian;
pub mod hessian;
pub mod integrals;
pub mod integrals_d;
pub mod linalg;
pub mod math;
pub mod method;
pub mod mindo3;
pub mod onecenter;
pub mod optimizer;
pub mod overlap;
pub mod overlap_numeric;
pub mod params;
pub mod repulsion;
pub mod rotations;
pub mod scf;
pub mod system;
pub mod xndo;
mod zdo_gradient;
pub mod zindo;

#[cfg(feature = "python")]
pub mod python;

pub use cndo_indo::{
    element as cndo_indo_element, run_cndo_indo, CndoIndoElement, CndoIndoOptions, CndoIndoResult,
};
pub use error::{Result, XndoError};
pub use gradient::{analytic_gradient, closed_form_gradient, numerical_gradient, GradientResult};
pub use hessian::{
    analytic_hessian, numerical_hessian, vibrational_analysis, vibrational_analysis_from_hessian,
    VibrationalModes,
};
pub use linalg::Matrix;
pub use math::{Mat3, Vec3};
pub use method::{methods_for_api, Method, MethodStatus};
pub use mindo3::{
    element as mindo3_element, pair as mindo3_pair, run_mindo3, Mindo3Element, Mindo3Options,
    Mindo3Result,
};
pub use optimizer::{optimize, OptOptions, OptResult, OptStep};
pub use params::{NddoElement, NddoParameters, PairParams};
pub use scf::{
    run_nddo_with_parameters, NddoCalculator, NddoOptions, NddoResult, Reference, ScfAccelerator,
};
pub use system::{symbol_to_z, z_to_symbol, Atom, Molecule};
pub use xndo::{
    run_gradient, run_hessian, run_method, run_nddo, CalculationResult, MethodGradientResult,
    MethodHessianResult,
};
pub use zindo::{
    run_zindo_s, uv_vis_spectrum, zindo_s_cis, zindo_s_cis_gradients, zindo_s_cis_hessians,
    zindo_s_cis_spin, zindo_s_ucis, zindo_s_ucis_gradients, zindo_s_ucis_hessians, CiContribution,
    CisSpin, ExcitedState, ExcitedStateGradient, ExcitedStateHessian, ZindoCisGradientResult,
    ZindoCisHessianResult, ZindoElement, ZindoOptions, ZindoParameters, ZindoResult, ZindoSpectrum,
    DEBYE_PER_E_BOHR, EV_TO_WAVENUMBER_CM1,
};
