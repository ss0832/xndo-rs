// SPDX-License-Identifier: GPL-3.0-or-later

//! High-level method dispatch for xndo-rs.

use crate::cndo_indo::{
    analytic_ground_gradient as cndo_indo_gradient, analytic_ground_hessian as cndo_indo_hessian,
    run_cndo_indo, CndoIndoOptions, CndoIndoResult,
};
use crate::error::{Result, XndoError};
use crate::gradient::closed_form_gradient;
use crate::hessian::analytic_hessian as nddo_hessian;
use crate::linalg::Matrix;
use crate::math::Vec3;
use crate::method::Method;
use crate::mindo3::{
    analytic_ground_gradient as mindo3_gradient, analytic_ground_hessian as mindo3_hessian,
    run_mindo3, Mindo3Options, Mindo3Result,
};
use crate::orbitals::OrbitalEnergies;
use crate::params::NddoParameters;
use crate::scf::{run_nddo_with_parameters, NddoOptions, NddoResult};
use crate::system::Molecule;
use crate::zindo::{
    analytic_ground_gradient as zindo_gradient, analytic_ground_hessian as zindo_hessian,
    run_zindo_s, ZindoOptions, ZindoParameters, ZindoResult,
};

/// Result type for the native fixed-geometry ground/reference-state engines.
#[derive(Clone, Debug)]
pub enum CalculationResult {
    CndoIndo(CndoIndoResult),
    Nddo(NddoResult),
    Mindo3(Mindo3Result),
    ZindoS(ZindoResult),
}

impl CalculationResult {
    /// The converged reference determinant's orbital energies, whichever engine
    /// produced them.
    ///
    /// The four result types spell the same quantity four ways
    /// (`mo_energies` vs `mo_energies_ev`, beta present or absent); this is the
    /// one accessor that does not care which engine ran.
    pub fn orbitals(&self) -> OrbitalEnergies {
        match self {
            CalculationResult::CndoIndo(r) => r.orbitals(),
            CalculationResult::Nddo(r) => r.orbitals(),
            CalculationResult::Mindo3(r) => r.orbitals(),
            CalculationResult::ZindoS(r) => r.orbitals(),
        }
    }

    /// Ground-state total energy in eV, whichever engine produced it.
    pub fn total_ev(&self) -> f64 {
        match self {
            CalculationResult::CndoIndo(r) => r.total_ev,
            CalculationResult::Nddo(r) => r.total_ev,
            CalculationResult::Mindo3(r) => r.total_ev,
            CalculationResult::ZindoS(r) => r.total_ev,
        }
    }

    /// Heat of formation in kcal/mol, for the methods parameterised to produce
    /// one. `None` for CNDO/2, INDO and ZINDO/S, which have no atomic heat
    /// terms, so their total energy is the comparable quantity.
    pub fn heat_of_formation_kcal(&self) -> Option<f64> {
        match self {
            CalculationResult::Nddo(r) => Some(r.heat_of_formation_kcal),
            CalculationResult::Mindo3(r) => Some(r.heat_of_formation_kcal),
            CalculationResult::CndoIndo(_) | CalculationResult::ZindoS(_) => None,
        }
    }
}

/// Orbital energies and the frontier pair for any method with a native engine.
///
/// A single-point run under the hood; use it when the orbitals are what you
/// want and the rest of the result is not.
pub fn run_orbitals(
    molecule: &Molecule,
    method: Method,
    options: &NddoOptions,
) -> Result<OrbitalEnergies> {
    Ok(run_method(molecule, method, options)?.orbitals())
}

/// Unified analytic ground-state nuclear gradient.
#[derive(Clone, Debug)]
pub struct MethodGradientResult {
    pub method: Method,
    pub energy_ev: f64,
    pub gradient: Vec<Vec3>,
    pub forces: Vec<Vec3>,
    pub iterations: usize,
    pub unrestricted: bool,
    pub heat_of_formation_kcal: Option<f64>,
}

/// Unified analytic ground-state Cartesian Hessian.
#[derive(Clone, Debug)]
pub struct MethodHessianResult {
    pub method: Method,
    pub hessian: Matrix,
}

fn cndo_options(options: &NddoOptions) -> CndoIndoOptions {
    CndoIndoOptions {
        charge: options.charge,
        multiplicity: options.multiplicity,
        reference: options.reference,
        max_scf: options.max_scf,
        e_tol_ev: options.e_tol,
        p_tol: options.p_tol,
        damping: if options.damping == 0.0 {
            0.20
        } else {
            options.damping
        },
        accelerator: options.accelerator,
        adiis_switch: options.adiis_switch,
        scf_memory_mb: options.scf_memory_mb,
    }
}

fn mindo_options(options: &NddoOptions) -> Mindo3Options {
    Mindo3Options {
        charge: options.charge,
        multiplicity: options.multiplicity,
        reference: options.reference,
        max_scf: options.max_scf,
        e_tol_ev: options.e_tol,
        p_tol: options.p_tol,
        damping: if options.damping == 0.0 {
            0.20
        } else {
            options.damping
        },
        accelerator: options.accelerator,
        adiis_switch: options.adiis_switch,
        scf_memory_mb: options.scf_memory_mb,
    }
}

fn zindo_options(options: &NddoOptions) -> ZindoOptions {
    ZindoOptions {
        charge: options.charge,
        multiplicity: options.multiplicity,
        reference: options.reference,
        max_scf: options.max_scf,
        e_tol_ev: options.e_tol,
        p_tol: options.p_tol,
        damping: if options.damping == 0.0 {
            0.20
        } else {
            options.damping
        },
        accelerator: options.accelerator,
        adiis_switch: options.adiis_switch,
        scf_memory_mb: options.scf_memory_mb,
        ..ZindoOptions::default()
    }
}

/// Run the analytic ground-state gradient for any native method that supports it.
pub fn run_gradient(
    molecule: &Molecule,
    method: Method,
    options: &NddoOptions,
) -> Result<MethodGradientResult> {
    let (energy_ev, gradient, iterations, unrestricted, heat_of_formation_kcal) = match method {
        Method::Cndo2 | Method::Indo => {
            let (result, gradient) = cndo_indo_gradient(molecule, method, &cndo_options(options))?;
            (
                result.total_ev,
                gradient,
                result.iterations,
                result.unrestricted,
                None,
            )
        }
        Method::Mndo | Method::MndoD => {
            let params = NddoParameters::cached_for_method(method)?;
            let result = closed_form_gradient(molecule, params, options)?;
            (
                result.energy_ev,
                result.gradient,
                result.scf.iterations,
                result.scf.unrestricted,
                Some(result.scf.heat_of_formation_kcal),
            )
        }
        Method::Mindo3 => {
            let (result, gradient) = mindo3_gradient(molecule, &mindo_options(options))?;
            (
                result.total_ev,
                gradient,
                result.iterations,
                result.unrestricted,
                Some(result.heat_of_formation_kcal),
            )
        }
        Method::ZindoS => {
            let params = ZindoParameters::cached()?;
            let (result, gradient) = zindo_gradient(molecule, params, &zindo_options(options))?;
            (
                result.total_ev,
                gradient,
                result.iterations,
                result.unrestricted,
                None,
            )
        }
        other => return Err(XndoError::InvalidInput(other.execution_error())),
    };
    let forces = gradient.iter().map(|g| *g * -1.0).collect();
    Ok(MethodGradientResult {
        method,
        energy_ev,
        gradient,
        forces,
        iterations,
        unrestricted,
        heat_of_formation_kcal,
    })
}

/// Run the analytic ground-state Hessian for any native derivative method.
pub fn run_hessian(
    molecule: &Molecule,
    method: Method,
    options: &NddoOptions,
) -> Result<MethodHessianResult> {
    let hessian = match method {
        Method::Cndo2 | Method::Indo => {
            cndo_indo_hessian(molecule, method, &cndo_options(options))?.2
        }
        Method::Mndo | Method::MndoD => {
            let params = NddoParameters::cached_for_method(method)?;
            nddo_hessian(molecule, params, options, 1.0e-3)?
        }
        Method::Mindo3 => mindo3_hessian(molecule, &mindo_options(options))?.2,
        Method::ZindoS => {
            zindo_hessian(
                molecule,
                ZindoParameters::cached()?,
                &zindo_options(options),
            )?
            .2
        }
        other => return Err(XndoError::InvalidInput(other.execution_error())),
    };
    Ok(MethodHessianResult { method, hessian })
}

/// Run MNDO or MNDO/d through the mature NDDO engine.
pub fn run_nddo(molecule: &Molecule, method: Method, options: &NddoOptions) -> Result<NddoResult> {
    if !matches!(method, Method::Mndo | Method::MndoD) {
        return Err(XndoError::InvalidInput(format!(
            "{} is not an MNDO/NDDO method; {}",
            method,
            method.implementation_note()
        )));
    }
    let params = NddoParameters::cached_for_method(method)?;
    run_nddo_with_parameters(molecule, params, options)
}

/// Unified single-point dispatcher for methods with a native numerical engine.
///
/// Spectroscopic excited states for ZINDO/S are obtained separately through
/// [`crate::zindo_s_cis`].  Registered legacy methods without a complete,
/// redistributable parameterization return an explicit error.
pub fn run_method(
    molecule: &Molecule,
    method: Method,
    nddo_options: &NddoOptions,
) -> Result<CalculationResult> {
    match method {
        Method::Cndo2 | Method::Indo => {
            let opts = cndo_options(nddo_options);
            Ok(CalculationResult::CndoIndo(run_cndo_indo(
                molecule, method, &opts,
            )?))
        }
        Method::Mndo | Method::MndoD => Ok(CalculationResult::Nddo(run_nddo(
            molecule,
            method,
            nddo_options,
        )?)),
        Method::Mindo3 => {
            let opts = mindo_options(nddo_options);
            Ok(CalculationResult::Mindo3(run_mindo3(molecule, &opts)?))
        }
        Method::ZindoS => {
            let params = ZindoParameters::cached()?;
            let opts = zindo_options(nddo_options);
            Ok(CalculationResult::ZindoS(run_zindo_s(
                molecule, params, &opts,
            )?))
        }
        other => Err(XndoError::InvalidInput(other.execution_error())),
    }
}
