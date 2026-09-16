// SPDX-License-Identifier: GPL-3.0-or-later

//! Python bindings for the `xndo-rs-python` distribution (`xndo_rs._native`).
//!
//! Coordinates are accepted in A.  Native energies/gradients are exposed in
//! atomic units as well as the convenient eV / eV A1 representations used by
//! ASE. Ground-state analytic gradients, geometry optimization, and Hessians
//! are available for every native method.

use crate::cndo_indo::{run_cndo_indo, CndoIndoOptions};
use crate::constants::{ANGSTROM_TO_BOHR, BOHR_TO_ANGSTROM, EV_TO_HARTREE};
use crate::data_tables;
use crate::method::Method;
use crate::mindo3::{run_mindo3, Mindo3Options};
use crate::molden::{molden_string, MoldenCoefficients};
use crate::optimizer::{optimize as opt_geom, OptOptions};
use crate::orbitals::OrbitalEnergies;
use crate::params::NddoParameters;
use crate::scf::{run_nddo_with_parameters, NddoOptions, Reference};
use crate::system::{Atom, Molecule};
use crate::xndo::{run_gradient, run_hessian, run_method};
use crate::zindo::{
    run_zindo_s, zindo_s_cis_gradients, zindo_s_cis_hessians, zindo_s_cis_spin, zindo_s_ucis,
    zindo_s_ucis_gradients, zindo_s_ucis_hessians, CisSpin, ZindoOptions, ZindoParameters,
    ZindoSpectrum, ZINDO_AU2EV,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

fn to_py_err(e: crate::error::XndoError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn parse_reference(reference: &str) -> PyResult<Reference> {
    match reference.to_ascii_lowercase().as_str() {
        "auto" | "" => Ok(Reference::Auto),
        "rhf" | "r" | "restricted" => Ok(Reference::Rhf),
        "uhf" | "u" | "unrestricted" => Ok(Reference::Uhf),
        other => Err(PyValueError::new_err(format!(
            "reference must be 'auto', 'rhf' or 'uhf' (got {other:?})"
        ))),
    }
}

fn parse_method(method: &str) -> PyResult<Method> {
    Method::parse(method).ok_or_else(|| {
        PyValueError::new_err(format!(
            "unknown xndo-rs method {method:?}; call available_methods() for registered names"
        ))
    })
}

fn build_molecule(
    numbers: &[u8],
    positions: &[Vec<f64>],
    charge: f64,
    mult: usize,
) -> PyResult<Molecule> {
    if numbers.len() != positions.len() {
        return Err(PyValueError::new_err(
            "numbers and positions length mismatch",
        ));
    }
    let mut atoms = Vec::with_capacity(numbers.len());
    for (z, p) in numbers.iter().zip(positions) {
        if p.len() != 3 {
            return Err(PyValueError::new_err(
                "each position must have 3 components",
            ));
        }
        atoms.push(Atom {
            z: *z,
            position: crate::math::Vec3::new(p[0], p[1], p[2]) * ANGSTROM_TO_BOHR,
        });
    }
    Ok(Molecule {
        atoms,
        charge,
        multiplicity: mult.max(1),
    })
}

fn nddo_options(charge: f64, multiplicity: usize, reference: Reference) -> NddoOptions {
    NddoOptions {
        charge,
        multiplicity,
        reference,
        ..NddoOptions::default()
    }
}

fn nddo_params(method: Method) -> PyResult<&'static NddoParameters> {
    match method {
        Method::Mndo | Method::MndoD => {
            NddoParameters::cached_for_method(method).map_err(to_py_err)
        }
        other => Err(PyValueError::new_err(format!(
            "{} is not handled by the MNDO/NDDO derivative engine: {}",
            other,
            other.implementation_note()
        ))),
    }
}

fn matrix_rows(m: &crate::linalg::Matrix) -> Vec<Vec<f64>> {
    (0..m.rows)
        .map(|i| (0..m.cols).map(|j| m[(i, j)]).collect())
        .collect()
}

/// Write the orbital block into a result dict.
///
/// Every engine's `single_point` dict carries the same keys, so a caller does
/// not have to know which one ran. `homo_ev` / `lumo_ev` are the frontier pair
/// over *both* spin channels; the per-channel values are there too, because for
/// an unrestricted reference the two channels come from different Fock
/// operators and a caller may legitimately want one diagram rather than the
/// combined pair. A frontier orbital that does not exist (a full or empty
/// shell) is `None`, never a substituted number.
fn set_orbital_items(d: &Bound<'_, PyDict>, o: &OrbitalEnergies) -> PyResult<()> {
    d.set_item("mo_energies_ev", o.alpha_ev.clone())?;
    d.set_item("mo_energies_beta_ev", o.beta_ev.clone())?;
    d.set_item("n_occ", o.n_alpha)?;
    d.set_item("n_alpha", o.n_alpha)?;
    d.set_item("n_beta", o.n_beta)?;
    d.set_item("homo_ev", o.homo_ev())?;
    d.set_item("lumo_ev", o.lumo_ev())?;
    d.set_item("homo_lumo_gap_ev", o.gap_ev())?;
    d.set_item("homo_alpha_ev", o.homo_alpha_ev())?;
    d.set_item("lumo_alpha_ev", o.lumo_alpha_ev())?;
    d.set_item("homo_beta_ev", o.homo_beta_ev())?;
    d.set_item("lumo_beta_ev", o.lumo_beta_ev())?;
    d.set_item("occupations", o.occupations_alpha())?;
    d.set_item("occupations_beta", o.occupations_beta())?;
    Ok(())
}

fn parse_cis_spins(state_type: &str) -> PyResult<Vec<CisSpin>> {
    match state_type.trim().to_ascii_lowercase().as_str() {
        "singlet" => Ok(vec![CisSpin::Singlet]),
        "triplet" => Ok(vec![CisSpin::Triplet]),
        "both" => Ok(vec![CisSpin::Singlet, CisSpin::Triplet]),
        _ => Err(PyValueError::new_err(
            "state_type must be exactly 'singlet', 'triplet', or 'both'",
        )),
    }
}

#[pyfunction]
#[pyo3(signature = (numbers, positions, charge=0.0, multiplicity=1, reference="auto", method="mndo"))]
fn single_point(
    py: Python<'_>,
    numbers: Vec<u8>,
    positions: Vec<Vec<f64>>,
    charge: f64,
    multiplicity: usize,
    reference: &str,
    method: &str,
) -> PyResult<PyObject> {
    let meth = parse_method(method)?;
    let mol = build_molecule(&numbers, &positions, charge, multiplicity)?;
    let d = PyDict::new(py);
    match meth {
        Method::Cndo2 | Method::Indo => {
            let opts = CndoIndoOptions {
                charge,
                multiplicity,
                reference: parse_reference(reference)?,
                ..CndoIndoOptions::default()
            };
            let r = run_cndo_indo(&mol, meth, &opts).map_err(to_py_err)?;
            d.set_item("method", meth.as_str())?;
            d.set_item("energy_hartree", r.total_ev * EV_TO_HARTREE)?;
            d.set_item("energy_ev", r.total_ev)?;
            d.set_item("electronic_ev", r.electronic_ev)?;
            d.set_item("core_ev", r.core_ev)?;
            d.set_item("charges", r.charges.clone())?;
            d.set_item(
                "dipole_debye",
                [r.dipole_debye[0], r.dipole_debye[1], r.dipole_debye[2]],
            )?;
            set_orbital_items(&d, &r.orbitals())?;
            d.set_item("unrestricted", r.unrestricted)?;
            if let Some(sd) = &r.spin_density {
                d.set_item("spin_density", matrix_rows(sd))?;
            } else {
                d.set_item("spin_density", py.None())?;
            }
            d.set_item("iterations", r.iterations)?;
            d.set_item("converged", r.converged)?;
        }
        Method::Mndo | Method::MndoD => {
            let params = nddo_params(meth)?;
            let r = run_nddo_with_parameters(
                &mol,
                params,
                &nddo_options(charge, multiplicity, parse_reference(reference)?),
            )
            .map_err(to_py_err)?;
            d.set_item("method", meth.as_str())?;
            d.set_item("energy_hartree", r.total_ev * EV_TO_HARTREE)?;
            d.set_item("energy_ev", r.total_ev)?;
            d.set_item("heat_of_formation_kcal", r.heat_of_formation_kcal)?;
            d.set_item("electronic_ev", r.electronic_ev)?;
            d.set_item("core_ev", r.core_ev)?;
            d.set_item("charges", r.charges.clone())?;
            d.set_item(
                "dipole_debye",
                [r.dipole_debye.x, r.dipole_debye.y, r.dipole_debye.z],
            )?;
            set_orbital_items(&d, &r.orbitals())?;
            d.set_item("iterations", r.iterations)?;
            d.set_item("converged", r.converged)?;
            d.set_item("unrestricted", r.unrestricted)?;
            if let Some(sd) = &r.spin_density {
                d.set_item("spin_density", matrix_rows(sd))?;
            } else {
                d.set_item("spin_density", py.None())?;
            }
        }
        Method::Mindo3 => {
            let opts = Mindo3Options {
                charge,
                multiplicity,
                reference: parse_reference(reference)?,
                ..Mindo3Options::default()
            };
            let r = run_mindo3(&mol, &opts).map_err(to_py_err)?;
            d.set_item("method", "MINDO/3")?;
            d.set_item("energy_hartree", r.total_ev * EV_TO_HARTREE)?;
            d.set_item("energy_ev", r.total_ev)?;
            d.set_item("heat_of_formation_kcal", r.heat_of_formation_kcal)?;
            d.set_item("electronic_ev", r.electronic_ev)?;
            d.set_item("core_ev", r.core_ev)?;
            d.set_item("charges", r.charges.clone())?;
            set_orbital_items(&d, &r.orbitals())?;
            d.set_item("unrestricted", r.unrestricted)?;
            if let Some(sd) = &r.spin_density {
                d.set_item("spin_density", matrix_rows(sd))?;
            } else {
                d.set_item("spin_density", py.None())?;
            }
            d.set_item("iterations", r.iterations)?;
            d.set_item("converged", r.converged)?;
        }
        Method::ZindoS => {
            let params = ZindoParameters::cached().map_err(to_py_err)?;
            let opts = ZindoOptions {
                charge,
                multiplicity,
                reference: parse_reference(reference)?,
                ..ZindoOptions::default()
            };
            let r = run_zindo_s(&mol, params, &opts).map_err(to_py_err)?;
            d.set_item("method", "ZINDO/S")?;
            d.set_item("energy_hartree", r.total_ev / ZINDO_AU2EV)?;
            d.set_item("energy_ev", r.total_ev)?;
            d.set_item("electronic_ev", r.electronic_ev)?;
            d.set_item("core_ev", r.core_ev)?;
            d.set_item("charges", r.charges.clone())?;
            d.set_item(
                "dipole_debye",
                [r.dipole_debye[0], r.dipole_debye[1], r.dipole_debye[2]],
            )?;
            set_orbital_items(&d, &r.orbitals())?;
            d.set_item("unrestricted", r.unrestricted)?;
            if let Some(spin_density) = &r.spin_density {
                d.set_item("spin_density", matrix_rows(spin_density))?;
            } else {
                d.set_item("spin_density", py.None())?;
            }
            d.set_item("iterations", r.iterations)?;
            d.set_item("converged", r.converged)?;
        }
        other => return Err(PyValueError::new_err(other.execution_error())),
    }
    Ok(d.into())
}

/// Orbital energies and the HOMO-LUMO pair for any method with a native engine.
///
/// Returns the orbital block of `single_point` on its own: the alpha and beta
/// spectra in eV, the occupation counts, the frontier pair over both spin
/// channels and per channel, and the gap.
///
/// `homo_ev` and `lumo_ev` are read over *both* spin channels, so for an
/// open-shell doublet the LUMO is usually the beta partner of the singly
/// occupied orbital rather than the lowest unoccupied alpha orbital. Use
/// `lumo_alpha_ev` if the alpha diagram alone is what you want. A frontier
/// orbital that does not exist is `None`.
#[pyfunction]
#[pyo3(signature = (numbers, positions, charge=0.0, multiplicity=1, reference="auto", method="mndo"))]
fn orbital_energies(
    py: Python<'_>,
    numbers: Vec<u8>,
    positions: Vec<Vec<f64>>,
    charge: f64,
    multiplicity: usize,
    reference: &str,
    method: &str,
) -> PyResult<PyObject> {
    let meth = parse_method(method)?;
    let mol = build_molecule(&numbers, &positions, charge, multiplicity)?;
    let result = run_method(
        &mol,
        meth,
        &nddo_options(charge, multiplicity, parse_reference(reference)?),
    )
    .map_err(to_py_err)?;
    let d = PyDict::new(py);
    d.set_item("method", meth.as_str())?;
    d.set_item("energy_ev", result.total_ev())?;
    set_orbital_items(&d, &result.orbitals())?;
    Ok(d.into())
}

/// A Molden wavefunction file, as a string.
///
/// The MO coefficients are back-transformed with `S^(-1/2)` so that they are
/// orthonormal over the Gaussians in the file. That matters: every engine here
/// assumes an orthonormal AO basis, while a Molden file describes real
/// Gaussians, which are not orthonormal, and a reader that took the raw
/// coefficients would compute a density that does not integrate to the electron
/// count. Pass `coefficients="raw"` to get the untransformed coefficients for
/// comparison with programs that make that identification; such a file says so
/// in its title.
///
/// The Slater basis is expanded as STO-6G (Stewart, *J. Chem. Phys.* **52**,
/// 431 (1970)), and d shells are written in Molden's `[5D]` order.
#[pyfunction]
#[pyo3(signature = (numbers, positions, charge=0.0, multiplicity=1, reference="auto",
                    method="mndo", coefficients="lowdin"))]
#[allow(clippy::too_many_arguments)]
fn molden(
    numbers: Vec<u8>,
    positions: Vec<Vec<f64>>,
    charge: f64,
    multiplicity: usize,
    reference: &str,
    method: &str,
    coefficients: &str,
) -> PyResult<String> {
    let meth = parse_method(method)?;
    let mol = build_molecule(&numbers, &positions, charge, multiplicity)?;
    let mode = match coefficients.trim().to_ascii_lowercase().as_str() {
        "lowdin" | "" => MoldenCoefficients::LowdinBackTransformed,
        "raw" | "raw_zdo" => MoldenCoefficients::RawZdo,
        other => {
            return Err(PyValueError::new_err(format!(
                "coefficients must be 'lowdin' or 'raw' (got {other:?})"
            )))
        }
    };
    molden_string(
        &mol,
        meth,
        &nddo_options(charge, multiplicity, parse_reference(reference)?),
        mode,
    )
    .map_err(to_py_err)
}

/// The third-party licence and attribution documents this build embeds.
///
/// One dict per document, with `path`, `role` and the verbatim `text`. These
/// are the notices Apache-2.0 4(c) and GPL-3.0 5(a) require to be carried, so
/// they travel with the installed wheel rather than only with the source tree.
#[pyfunction]
fn third_party_licenses(py: Python<'_>) -> PyResult<PyObject> {
    let out = PyList::empty(py);
    for doc in crate::licenses::third_party_licenses() {
        let d = PyDict::new(py);
        d.set_item("path", doc.path)?;
        d.set_item("role", doc.role)?;
        d.set_item("text", doc.text)?;
        out.append(d)?;
    }
    Ok(out.into())
}

#[pyfunction]
#[pyo3(signature = (numbers, positions, charge=0.0, multiplicity=1, reference="auto", method="mndo"))]
fn gradient(
    py: Python<'_>,
    numbers: Vec<u8>,
    positions: Vec<Vec<f64>>,
    charge: f64,
    multiplicity: usize,
    reference: &str,
    method: &str,
) -> PyResult<PyObject> {
    let meth = parse_method(method)?;
    if !meth.supports_gradient() {
        return Err(PyValueError::new_err(meth.execution_error()));
    }
    let mol = build_molecule(&numbers, &positions, charge, multiplicity)?;
    let g = run_gradient(
        &mol,
        meth,
        &nddo_options(charge, multiplicity, parse_reference(reference)?),
    )
    .map_err(to_py_err)?;
    let grad_au: Vec<[f64; 3]> = g
        .gradient
        .iter()
        .map(|v| {
            [
                v.x * EV_TO_HARTREE,
                v.y * EV_TO_HARTREE,
                v.z * EV_TO_HARTREE,
            ]
        })
        .collect();
    let grad_ev_ang: Vec<[f64; 3]> = g
        .gradient
        .iter()
        .map(|v| {
            [
                v.x * ANGSTROM_TO_BOHR,
                v.y * ANGSTROM_TO_BOHR,
                v.z * ANGSTROM_TO_BOHR,
            ]
        })
        .collect();
    let d = PyDict::new(py);
    d.set_item("method", meth.as_str())?;
    d.set_item("energy_hartree", g.energy_ev * EV_TO_HARTREE)?;
    d.set_item("energy_ev", g.energy_ev)?;
    if let Some(heat) = g.heat_of_formation_kcal {
        d.set_item("heat_of_formation_kcal", heat)?;
    }
    d.set_item("iterations", g.iterations)?;
    d.set_item("unrestricted", g.unrestricted)?;
    d.set_item("gradient_hartree_per_bohr", grad_au)?;
    d.set_item("gradient_ev_per_angstrom", grad_ev_ang)?;
    Ok(d.into())
}

#[pyfunction]
#[pyo3(signature = (numbers, positions, charge=0.0, multiplicity=1, reference="auto", method="mndo"))]
fn forces(
    py: Python<'_>,
    numbers: Vec<u8>,
    positions: Vec<Vec<f64>>,
    charge: f64,
    multiplicity: usize,
    reference: &str,
    method: &str,
) -> PyResult<PyObject> {
    let meth = parse_method(method)?;
    if !meth.supports_gradient() {
        return Err(PyValueError::new_err(meth.execution_error()));
    }
    let mol = build_molecule(&numbers, &positions, charge, multiplicity)?;
    let g = run_gradient(
        &mol,
        meth,
        &nddo_options(charge, multiplicity, parse_reference(reference)?),
    )
    .map_err(to_py_err)?;
    let f_au: Vec<[f64; 3]> = g
        .forces
        .iter()
        .map(|v| {
            [
                v.x * EV_TO_HARTREE,
                v.y * EV_TO_HARTREE,
                v.z * EV_TO_HARTREE,
            ]
        })
        .collect();
    let f_ev_ang: Vec<[f64; 3]> = g
        .forces
        .iter()
        .map(|v| {
            [
                v.x * ANGSTROM_TO_BOHR,
                v.y * ANGSTROM_TO_BOHR,
                v.z * ANGSTROM_TO_BOHR,
            ]
        })
        .collect();
    let d = PyDict::new(py);
    d.set_item("method", meth.as_str())?;
    d.set_item("energy_hartree", g.energy_ev * EV_TO_HARTREE)?;
    d.set_item("energy_ev", g.energy_ev)?;
    if let Some(heat) = g.heat_of_formation_kcal {
        d.set_item("heat_of_formation_kcal", heat)?;
    }
    d.set_item("iterations", g.iterations)?;
    d.set_item("unrestricted", g.unrestricted)?;
    d.set_item("forces_hartree_per_bohr", f_au)?;
    d.set_item("forces_ev_per_angstrom", f_ev_ang)?;
    Ok(d.into())
}

#[pyfunction]
#[pyo3(signature = (numbers, positions, charge=0.0, multiplicity=1, reference="auto", method="mndo", max_iter=200, gtol=1.0e-3, history=8))]
#[allow(clippy::too_many_arguments)]
fn optimize(
    py: Python<'_>,
    numbers: Vec<u8>,
    positions: Vec<Vec<f64>>,
    charge: f64,
    multiplicity: usize,
    reference: &str,
    method: &str,
    max_iter: usize,
    gtol: f64,
    history: usize,
) -> PyResult<PyObject> {
    let meth = parse_method(method)?;
    if !meth.supports_gradient() {
        return Err(PyValueError::new_err(meth.execution_error()));
    }
    let mol = build_molecule(&numbers, &positions, charge, multiplicity)?;
    let res = opt_geom(
        &mol,
        meth,
        &nddo_options(charge, multiplicity, parse_reference(reference)?),
        &OptOptions {
            max_iter,
            gtol,
            history,
        },
    )
    .map_err(to_py_err)?;
    let coords: Vec<[f64; 3]> = res
        .molecule
        .atoms
        .iter()
        .map(|a| {
            let p = a.position * BOHR_TO_ANGSTROM;
            [p.x, p.y, p.z]
        })
        .collect();
    let d = PyDict::new(py);
    d.set_item("method", meth.as_str())?;
    d.set_item("positions_angstrom", coords)?;
    d.set_item("energy_hartree", res.energy_ev * EV_TO_HARTREE)?;
    d.set_item("energy_ev", res.energy_ev)?;
    if let Some(heat) = res.heat_of_formation_kcal {
        d.set_item("heat_of_formation_kcal", heat)?;
    }
    d.set_item("converged", res.converged)?;
    d.set_item("iterations", res.iterations)?;
    Ok(d.into())
}

#[pyfunction]
#[pyo3(signature = (numbers, positions, charge=0.0, multiplicity=1, reference="auto", method="mndo"))]
fn frequencies(
    py: Python<'_>,
    numbers: Vec<u8>,
    positions: Vec<Vec<f64>>,
    charge: f64,
    multiplicity: usize,
    reference: &str,
    method: &str,
) -> PyResult<PyObject> {
    let meth = parse_method(method)?;
    if !meth.supports_hessian() {
        return Err(PyValueError::new_err(meth.execution_error()));
    }
    let mol = build_molecule(&numbers, &positions, charge, multiplicity)?;
    let h = run_hessian(
        &mol,
        meth,
        &nddo_options(charge, multiplicity, parse_reference(reference)?),
    )
    .map_err(to_py_err)?;
    let vib =
        crate::hessian::vibrational_analysis_from_hessian(&mol, h.hessian).map_err(to_py_err)?;
    let d = PyDict::new(py);
    d.set_item("method", meth.as_str())?;
    d.set_item("frequencies_cm", vib.frequencies_cm)?;
    d.set_item("eigenvalues", vib.eigenvalues)?;
    Ok(d.into())
}

#[pyfunction]
#[pyo3(signature = (numbers, positions, charge=0.0, multiplicity=1, reference="auto", method="mndo"))]
fn hessian(
    py: Python<'_>,
    numbers: Vec<u8>,
    positions: Vec<Vec<f64>>,
    charge: f64,
    multiplicity: usize,
    reference: &str,
    method: &str,
) -> PyResult<PyObject> {
    let meth = parse_method(method)?;
    if !meth.supports_hessian() {
        return Err(PyValueError::new_err(meth.execution_error()));
    }
    let mol = build_molecule(&numbers, &positions, charge, multiplicity)?;
    let result = run_hessian(
        &mol,
        meth,
        &nddo_options(charge, multiplicity, parse_reference(reference)?),
    )
    .map_err(to_py_err)?;
    let h = result.hessian;
    let mut rows = Vec::with_capacity(h.rows);
    for i in 0..h.rows {
        let mut row = Vec::with_capacity(h.cols);
        for j in 0..h.cols {
            row.push(h[(i, j)] * EV_TO_HARTREE);
        }
        rows.push(row);
    }
    let d = PyDict::new(py);
    d.set_item("method", meth.as_str())?;
    d.set_item("hessian_hartree_per_bohr2", rows)?;
    d.set_item("ndof", h.rows)?;
    Ok(d.into())
}

#[pyfunction]
#[pyo3(signature = (numbers, positions, charge=0.0, n_states=10, active_occupied=None, active_virtual=None, state_type="singlet", multiplicity=1, reference="auto"))]
#[allow(clippy::too_many_arguments)]
fn excited_state_gradients(
    py: Python<'_>,
    numbers: Vec<u8>,
    positions: Vec<Vec<f64>>,
    charge: f64,
    n_states: usize,
    active_occupied: Option<usize>,
    active_virtual: Option<usize>,
    state_type: &str,
    multiplicity: usize,
    reference: &str,
) -> PyResult<PyObject> {
    let parsed_reference = parse_reference(reference)?;
    let mol = build_molecule(&numbers, &positions, charge, multiplicity)?;
    let params = ZindoParameters::cached().map_err(to_py_err)?;
    let options = ZindoOptions {
        charge,
        multiplicity,
        reference: parsed_reference,
        n_states,
        active_occupied,
        active_virtual,
        ..ZindoOptions::default()
    };
    let states = PyList::empty(py);
    let mut ground_energy = None;
    let results = if parsed_reference == Reference::Uhf || multiplicity != 1 {
        vec![zindo_s_ucis_gradients(&mol, params, &options).map_err(to_py_err)?]
    } else {
        parse_cis_spins(state_type)?
            .into_iter()
            .map(|spin| zindo_s_cis_gradients(&mol, params, &options, spin).map_err(to_py_err))
            .collect::<PyResult<Vec<_>>>()?
    };
    for result in results {
        ground_energy = Some(result.ground.total_ev);
        for state in result.states {
            let d = PyDict::new(py);
            d.set_item("root", state.root)?;
            d.set_item("spin", state.spin.as_str())?;
            d.set_item("excitation_energy_ev", state.excitation_energy_ev)?;
            d.set_item("state_total_energy_ev", state.state_total_energy_ev)?;
            d.set_item(
                "excitation_gradient_hartree_per_bohr",
                state
                    .excitation_gradient
                    .iter()
                    .map(|v| {
                        [
                            v.x * EV_TO_HARTREE,
                            v.y * EV_TO_HARTREE,
                            v.z * EV_TO_HARTREE,
                        ]
                    })
                    .collect::<Vec<_>>(),
            )?;
            d.set_item(
                "excitation_gradient_ev_per_angstrom",
                state
                    .excitation_gradient
                    .iter()
                    .map(|v| {
                        [
                            v.x * ANGSTROM_TO_BOHR,
                            v.y * ANGSTROM_TO_BOHR,
                            v.z * ANGSTROM_TO_BOHR,
                        ]
                    })
                    .collect::<Vec<_>>(),
            )?;
            d.set_item(
                "state_gradient_hartree_per_bohr",
                state
                    .state_gradient
                    .iter()
                    .map(|v| {
                        [
                            v.x * EV_TO_HARTREE,
                            v.y * EV_TO_HARTREE,
                            v.z * EV_TO_HARTREE,
                        ]
                    })
                    .collect::<Vec<_>>(),
            )?;
            d.set_item(
                "state_gradient_ev_per_angstrom",
                state
                    .state_gradient
                    .iter()
                    .map(|v| {
                        [
                            v.x * ANGSTROM_TO_BOHR,
                            v.y * ANGSTROM_TO_BOHR,
                            v.z * ANGSTROM_TO_BOHR,
                        ]
                    })
                    .collect::<Vec<_>>(),
            )?;
            states.append(d)?;
        }
    }
    let out = PyDict::new(py);
    out.set_item("method", "ZINDO/S")?;
    out.set_item(
        "reference",
        if parsed_reference == Reference::Uhf || multiplicity != 1 {
            "UHF"
        } else {
            "RHF"
        },
    )?;
    out.set_item("ground_energy_ev", ground_energy.unwrap_or(f64::NAN))?;
    out.set_item("states", states)?;
    Ok(out.into())
}

#[pyfunction]
#[pyo3(signature = (numbers, positions, charge=0.0, n_states=10, active_occupied=None, active_virtual=None, state_type="singlet", multiplicity=1, reference="auto"))]
#[allow(clippy::too_many_arguments)]
fn excited_state_hessians(
    py: Python<'_>,
    numbers: Vec<u8>,
    positions: Vec<Vec<f64>>,
    charge: f64,
    n_states: usize,
    active_occupied: Option<usize>,
    active_virtual: Option<usize>,
    state_type: &str,
    multiplicity: usize,
    reference: &str,
) -> PyResult<PyObject> {
    let parsed_reference = parse_reference(reference)?;
    let mol = build_molecule(&numbers, &positions, charge, multiplicity)?;
    let params = ZindoParameters::cached().map_err(to_py_err)?;
    let options = ZindoOptions {
        charge,
        multiplicity,
        reference: parsed_reference,
        n_states,
        active_occupied,
        active_virtual,
        ..ZindoOptions::default()
    };
    let states = PyList::empty(py);
    let scale_ev_ang2 = ANGSTROM_TO_BOHR * ANGSTROM_TO_BOHR;
    let results = if parsed_reference == Reference::Uhf || multiplicity != 1 {
        vec![zindo_s_ucis_hessians(&mol, params, &options).map_err(to_py_err)?]
    } else {
        parse_cis_spins(state_type)?
            .into_iter()
            .map(|spin| zindo_s_cis_hessians(&mol, params, &options, spin).map_err(to_py_err))
            .collect::<PyResult<Vec<_>>>()?
    };
    for result in results {
        for state in result.states {
            let d = PyDict::new(py);
            d.set_item("root", state.root)?;
            d.set_item("spin", state.spin.as_str())?;
            d.set_item("excitation_energy_ev", state.excitation_energy_ev)?;
            d.set_item("state_total_energy_ev", state.state_total_energy_ev)?;
            d.set_item(
                "excitation_hessian_hartree_per_bohr2",
                matrix_rows(&state.excitation_hessian)
                    .into_iter()
                    .map(|row| {
                        row.into_iter()
                            .map(|v| v * EV_TO_HARTREE)
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>(),
            )?;
            d.set_item(
                "excitation_hessian_ev_per_angstrom2",
                matrix_rows(&state.excitation_hessian)
                    .into_iter()
                    .map(|row| {
                        row.into_iter()
                            .map(|v| v * scale_ev_ang2)
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>(),
            )?;
            d.set_item(
                "state_hessian_hartree_per_bohr2",
                matrix_rows(&state.state_hessian)
                    .into_iter()
                    .map(|row| {
                        row.into_iter()
                            .map(|v| v * EV_TO_HARTREE)
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>(),
            )?;
            d.set_item(
                "state_hessian_ev_per_angstrom2",
                matrix_rows(&state.state_hessian)
                    .into_iter()
                    .map(|row| {
                        row.into_iter()
                            .map(|v| v * scale_ev_ang2)
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>(),
            )?;
            states.append(d)?;
        }
    }
    let out = PyDict::new(py);
    out.set_item("method", "ZINDO/S")?;
    out.set_item(
        "reference",
        if parsed_reference == Reference::Uhf || multiplicity != 1 {
            "UHF"
        } else {
            "RHF"
        },
    )?;
    out.set_item("states", states)?;
    Ok(out.into())
}

/// ZINDO/S vertical excited-state properties (spin-adapted RHF-CIS or UCIS).
///
/// Accepted method strings: `zindo/s`, `zindos`, `zindo`, `indo/s`, `indos`,
/// `state_type` accepts `singlet`, `triplet`, or `both` for RHF; a UHF
/// reference selects the unrestricted spin-orbital CIS sector.
#[pyfunction]
#[pyo3(signature = (numbers, positions, charge=0.0, n_states=10, active_occupied=None, active_virtual=None, method="zindo/s", state_type="singlet", multiplicity=1, reference="auto"))]
#[allow(clippy::too_many_arguments)]
fn excited_states(
    py: Python<'_>,
    numbers: Vec<u8>,
    positions: Vec<Vec<f64>>,
    charge: f64,
    n_states: usize,
    active_occupied: Option<usize>,
    active_virtual: Option<usize>,
    method: &str,
    state_type: &str,
    multiplicity: usize,
    reference: &str,
) -> PyResult<PyObject> {
    let meth = parse_method(method)?;
    if meth != Method::ZindoS {
        return Err(PyValueError::new_err(format!(
            "excited_states() in v0.2.4 executes ZINDO/S only; {meth} is not an enabled excited-state engine"
        )));
    }
    let parsed_reference = parse_reference(reference)?;
    let mol = build_molecule(&numbers, &positions, charge, multiplicity)?;
    let params = ZindoParameters::cached().map_err(to_py_err)?;
    let opts = ZindoOptions {
        charge,
        multiplicity,
        reference: parsed_reference,
        n_states,
        active_occupied,
        active_virtual,
        ..ZindoOptions::default()
    };
    let unrestricted = parsed_reference == Reference::Uhf || multiplicity != 1;
    let requested = if unrestricted {
        "unrestricted".to_string()
    } else {
        state_type.trim().to_ascii_lowercase()
    };
    let spectra: Vec<ZindoSpectrum> = if unrestricted {
        vec![zindo_s_ucis(&mol, params, &opts).map_err(to_py_err)?]
    } else {
        parse_cis_spins(&requested)?
            .into_iter()
            .map(|spin| zindo_s_cis_spin(&mol, params, &opts, spin).map_err(to_py_err))
            .collect::<PyResult<Vec<_>>>()?
    };
    let first = spectra
        .first()
        .ok_or_else(|| PyValueError::new_err("no CIS spin sector requested"))?;
    let d = PyDict::new(py);
    d.set_item("method", "ZINDO/S")?;
    d.set_item("reference", if unrestricted { "UHF" } else { "RHF" })?;
    d.set_item("state_type", requested)?;
    d.set_item("ground_energy_hartree", first.ground.total_ev / ZINDO_AU2EV)?;
    d.set_item("ground_energy_ev", first.ground.total_ev)?;
    d.set_item("ground_charges", first.ground.charges.clone())?;
    d.set_item("ground_dipole_au", first.ground_dipole_au)?;
    d.set_item("ground_dipole_debye", first.ground_dipole_debye)?;
    d.set_item("mo_energies_ev", first.ground.mo_energies_ev.clone())?;
    d.set_item(
        "mo_energies_beta_ev",
        first.ground.mo_energies_beta_ev.clone(),
    )?;
    let states = PyList::empty(py);
    for spec in &spectra {
        for (idx, st) in spec.states.iter().enumerate() {
            let sd = PyDict::new(py);
            sd.set_item("state", idx + 1)?;
            sd.set_item("spin", st.spin.as_str())?;
            sd.set_item("spin_multiplicity", st.spin_multiplicity)?;
            sd.set_item("s2_expectation", st.s2_expectation)?;
            sd.set_item("energy_ev", st.energy_ev)?;
            sd.set_item("state_total_energy_ev", st.state_total_energy_ev)?;
            sd.set_item("state_total_energy_hartree", st.state_total_energy_hartree)?;
            sd.set_item("energy_cm1", st.energy_cm1)?;
            sd.set_item("wavelength_nm", st.wavelength_nm)?;
            sd.set_item("oscillator_strength", st.oscillator_strength)?;
            sd.set_item("transition_dipole_au", st.transition_dipole_au)?;
            sd.set_item("transition_dipole_debye", st.transition_dipole_debye)?;
            sd.set_item(
                "transition_dipole_magnitude_au",
                st.transition_dipole_magnitude_au,
            )?;
            sd.set_item("dipole_strength_au2", st.dipole_strength_au2)?;
            sd.set_item("permanent_dipole_au", st.permanent_dipole_au)?;
            sd.set_item("permanent_dipole_debye", st.permanent_dipole_debye)?;
            sd.set_item(
                "permanent_dipole_magnitude_debye",
                st.permanent_dipole_magnitude_debye,
            )?;
            sd.set_item("difference_dipole_au", st.difference_dipole_au)?;
            sd.set_item("difference_dipole_debye", st.difference_dipole_debye)?;
            sd.set_item(
                "difference_dipole_magnitude_debye",
                st.difference_dipole_magnitude_debye,
            )?;
            sd.set_item("charges", st.charges.clone())?;
            sd.set_item("hole_population", st.hole_population.clone())?;
            sd.set_item("electron_population", st.electron_population.clone())?;
            sd.set_item("hole_centroid_angstrom", st.hole_centroid_angstrom)?;
            sd.set_item("electron_centroid_angstrom", st.electron_centroid_angstrom)?;
            sd.set_item(
                "charge_transfer_distance_angstrom",
                st.charge_transfer_distance_angstrom,
            )?;
            let dom: Vec<(usize, usize, f64)> = st
                .dominant
                .iter()
                .map(|c| (c.occupied, c.virtual_orbital, c.coefficient))
                .collect();
            sd.set_item("dominant_configurations", dom)?;
            states.append(sd)?;
        }
    }
    d.set_item("states", states)?;
    Ok(d.into())
}

/// Alias emphasizing that excited_states() returns state properties, not only
/// an excitation-energy list.
#[pyfunction]
#[pyo3(signature = (numbers, positions, charge=0.0, n_states=10, active_occupied=None, active_virtual=None, method="zindo/s", state_type="singlet", multiplicity=1, reference="auto"))]
#[allow(clippy::too_many_arguments)]
fn excited_properties(
    py: Python<'_>,
    numbers: Vec<u8>,
    positions: Vec<Vec<f64>>,
    charge: f64,
    n_states: usize,
    active_occupied: Option<usize>,
    active_virtual: Option<usize>,
    method: &str,
    state_type: &str,
    multiplicity: usize,
    reference: &str,
) -> PyResult<PyObject> {
    excited_states(
        py,
        numbers,
        positions,
        charge,
        n_states,
        active_occupied,
        active_virtual,
        method,
        state_type,
        multiplicity,
        reference,
    )
}

/// Explicit UV-visible stick-spectrum entry point for RHF-CIS or UHF-UCIS.
#[pyfunction]
#[pyo3(signature = (numbers, positions, charge=0.0, n_states=10, active_occupied=None, active_virtual=None, method="zindo/s", multiplicity=1, reference="auto"))]
#[allow(clippy::too_many_arguments)]
fn uv_vis_spectrum(
    py: Python<'_>,
    numbers: Vec<u8>,
    positions: Vec<Vec<f64>>,
    charge: f64,
    n_states: usize,
    active_occupied: Option<usize>,
    active_virtual: Option<usize>,
    method: &str,
    multiplicity: usize,
    reference: &str,
) -> PyResult<PyObject> {
    excited_states(
        py,
        numbers,
        positions,
        charge,
        n_states,
        active_occupied,
        active_virtual,
        method,
        "singlet",
        multiplicity,
        reference,
    )
}

#[pyfunction]
fn parameter_dataset(name: &str) -> PyResult<String> {
    let key = name
        .trim()
        .to_ascii_lowercase()
        .replace(['-', '/', ' '], "_");
    let text = match key.as_str() {
        "mndo" | "openmopac_mndo" => data_tables::MNDO_PARAM_CSV,
        "mndo_pair" | "mndo_pairs" => data_tables::MNDO_PAIR_CSV,
        "mndod" | "mndo_d" | "openmopac_mndod" => data_tables::MNDOD_PARAM_CSV,
        "mndod_pair" | "mndod_pairs" | "mndo_d_pair" | "mndo_d_pairs" => data_tables::MNDOD_PAIR_CSV,
        "zindo_s" | "openmopac_zindo_s" => data_tables::ZINDO_S_PARAM_CSV,
        "mindo3" | "mindo_3" | "mopac7_mindo3" => data_tables::MINDO3_PARAM_CSV,
        "mindo3_pair" | "mindo3_pairs" | "mindo_3_pair" | "mindo_3_pairs" => data_tables::MINDO3_PAIR_PARAM_CSV,
        "molds_cndo2_indo" | "cndo2_indo" => data_tables::MOLDS_CNDO2_INDO_PARAM_CSV,
        "molds_zindo_s" => data_tables::MOLDS_ZINDO_S_PARAM_CSV,
        "catalog" | "legacy_parameter_catalog" => data_tables::LEGACY_PARAMETER_CATALOG_CSV,
        "manifest" | "legacy_parameter_manifest" => data_tables::LEGACY_PARAMETER_MANIFEST_SHA256,
        other => return Err(PyValueError::new_err(format!(
            "unknown parameter dataset {other:?}; call parameter_datasets() for available dataset IDs"
        ))),
    };
    Ok(text.to_string())
}

#[pyfunction]
fn parameter_datasets(py: Python<'_>) -> PyResult<PyObject> {
    let table = data_tables::CsvTable::parse(data_tables::LEGACY_PARAMETER_CATALOG_CSV)
        .ok_or_else(|| PyValueError::new_err("embedded legacy parameter catalog is invalid"))?;
    let out = PyList::empty(py);
    for row in &table.rows {
        let d = PyDict::new(py);
        for (idx, key) in table.header.iter().enumerate() {
            d.set_item(key.as_str(), row.get(idx).map(String::as_str).unwrap_or(""))?;
        }
        out.append(d)?;
    }
    Ok(out.into())
}

#[pyfunction]
fn available_methods(py: Python<'_>) -> PyResult<PyObject> {
    let out = PyList::empty(py);
    for method in Method::ALL {
        let d = PyDict::new(py);
        d.set_item("name", method.as_str())?;
        d.set_item("family", method.family())?;
        d.set_item("status", method.status().as_str())?;
        d.set_item("energy", method.supports_energy())?;
        d.set_item("gradient", method.supports_gradient())?;
        d.set_item("hessian", method.supports_hessian())?;
        d.set_item("uhf", method.supports_uhf())?;
        d.set_item("spectrum", method.supports_spectrum())?;
        d.set_item("excited_properties", method.supports_excited_properties())?;
        d.set_item("accepted_strings", method.accepted_strings().to_vec())?;
        d.set_item("apis", method.api_names().to_vec())?;
        d.set_item("note", method.implementation_note())?;
        out.append(d)?;
    }
    Ok(out.into())
}

/// Return an explicit API -> accepted method-string map for v0.2.4.
#[pyfunction]
fn api_methods(py: Python<'_>) -> PyResult<PyObject> {
    let apis = [
        "single_point",
        "gradient",
        "forces",
        "optimize",
        "hessian",
        "frequencies",
        "excited_states",
        "excited_properties",
        "uv_vis_spectrum",
        "excited_state_gradients",
        "excited_state_hessians",
    ];
    let out = PyDict::new(py);
    for api in apis {
        let rows = PyList::empty(py);
        for method in Method::ALL {
            let enabled = method.api_names().contains(&api);
            if !enabled {
                continue;
            }
            let d = PyDict::new(py);
            d.set_item("method", method.as_str())?;
            d.set_item("accepted_strings", method.accepted_strings().to_vec())?;
            d.set_item("uhf", method.supports_uhf())?;
            let excited_api = matches!(
                api,
                "excited_states"
                    | "excited_properties"
                    | "uv_vis_spectrum"
                    | "excited_state_gradients"
                    | "excited_state_hessians"
            );
            d.set_item("reference_argument", true)?;
            if method.supports_uhf() {
                d.set_item("reference_strings", vec!["auto", "rhf", "uhf"])?;
                d.set_item("fixed_reference", py.None())?;
            } else {
                d.set_item("reference_strings", Vec::<&str>::new())?;
                d.set_item("fixed_reference", py.None())?;
            }
            if api == "uv_vis_spectrum" && method == Method::ZindoS {
                d.set_item("state_type_strings", vec!["singlet", "unrestricted"])?;
            } else if excited_api && method == Method::ZindoS {
                d.set_item(
                    "state_type_strings",
                    vec!["singlet", "triplet", "both", "unrestricted"],
                )?;
            } else {
                d.set_item("state_type_strings", Vec::<&str>::new())?;
            }
            rows.append(d)?;
        }
        out.set_item(api, rows)?;
    }
    Ok(out.into())
}

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(single_point, m)?)?;
    m.add_function(wrap_pyfunction!(orbital_energies, m)?)?;
    m.add_function(wrap_pyfunction!(molden, m)?)?;
    m.add_function(wrap_pyfunction!(third_party_licenses, m)?)?;
    m.add_function(wrap_pyfunction!(gradient, m)?)?;
    m.add_function(wrap_pyfunction!(forces, m)?)?;
    m.add_function(wrap_pyfunction!(optimize, m)?)?;
    m.add_function(wrap_pyfunction!(frequencies, m)?)?;
    m.add_function(wrap_pyfunction!(hessian, m)?)?;
    m.add_function(wrap_pyfunction!(excited_states, m)?)?;
    m.add_function(wrap_pyfunction!(excited_state_gradients, m)?)?;
    m.add_function(wrap_pyfunction!(excited_state_hessians, m)?)?;
    m.add_function(wrap_pyfunction!(excited_properties, m)?)?;
    m.add_function(wrap_pyfunction!(uv_vis_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(api_methods, m)?)?;
    m.add_function(wrap_pyfunction!(parameter_dataset, m)?)?;
    m.add_function(wrap_pyfunction!(parameter_datasets, m)?)?;
    m.add_function(wrap_pyfunction!(available_methods, m)?)?;
    Ok(())
}
