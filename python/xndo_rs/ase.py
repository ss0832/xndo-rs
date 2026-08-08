# SPDX-License-Identifier: GPL-3.0-or-later
"""ASE calculators for native XNDO gradients and ZINDO/S spectroscopy."""
from __future__ import annotations
import numpy as np
try:
    from ase.calculators.calculator import Calculator, all_changes
    from ase.units import Bohr, Hartree
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "The xndo-rs ASE calculator requires ASE. Install with "
        "`pip install xndo-rs-python[ase]`."
    ) from exc
from . import native

_DEBYE_TO_E_ANGSTROM = 0.2081943


class XNDO(Calculator):
    """ASE calculator for all native analytic ground-state gradients.

    Exact ``method`` strings: MNDO = ``mndo``; MNDO/d = ``mndod``,
    ``mndo/d``, or ``mndo-d``. ``reference`` accepts ``auto``, ``rhf``, or
    ``uhf``.

    CNDO/2, INDO, MNDO, MNDO/d, MINDO/3, and ZINDO/S RHF/UHF ground states are
    accepted. Excited-state ZINDO/S gradients are exposed by the dedicated
    excited-state API.
    """
    implemented_properties = [
        "energy", "forces", "charges", "dipole", "heat_of_formation_kcal", "hessian"
    ]

    def __init__(self, charge: int = 0, multiplicity: int = 1,
                 reference: str = "auto", method: str = "mndo", **kwargs):
        super().__init__(**kwargs)
        key = str(method).lower().replace("_", "").replace("-", "")
        if key not in {
            "cndo", "cndo2", "cndo/2", "indo", "mndo", "mndod", "mndo/d",
            "mindo", "mindo3", "mindo/3", "zindo", "zindos", "zindo/s",
            "indos", "indo/s",
        }:
            raise ValueError(
                "xndo_rs.ase requires a native gradient method: CNDO/2, INDO, "
                "MNDO, MNDO/d, MINDO/3, or ZINDO/S"
            )
        self.charge = int(charge)
        self.multiplicity = int(multiplicity)
        self.reference = str(reference)
        self.method = str(method)

    def _args(self, atoms):
        if atoms is None:
            raise RuntimeError("XNDO has no Atoms object yet")
        return (atoms.get_atomic_numbers(), atoms.get_positions(), self.charge,
                self.multiplicity, self.reference, self.method)

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        sp = native.single_point(*self._args(self.atoms))
        self.results["energy"] = float(sp["energy_ev"])
        self.results["charges"] = np.asarray(sp["charges"], dtype=float)
        if "dipole_debye" in sp:
            self.results["dipole"] = np.asarray(sp["dipole_debye"], dtype=float) * _DEBYE_TO_E_ANGSTROM
        if "heat_of_formation_kcal" in sp:
            self.results["heat_of_formation_kcal"] = float(sp["heat_of_formation_kcal"])
        if "forces" in properties:
            fr = native.forces(*self._args(self.atoms))
            self.results["forces"] = np.asarray(fr["forces_ev_per_angstrom"], dtype=float)
        if "hessian" in properties:
            self.results["hessian"] = self.get_hessian(self.atoms)

    def get_hessian(self, atoms=None):
        atoms = atoms if atoms is not None else self.atoms
        h = native.hessian(*self._args(atoms))
        # Hartree/Bohr2 -> eV/A2
        return np.asarray(h["hessian_hartree_per_bohr2"], dtype=float) * Hartree / (Bohr * Bohr)

    def get_frequencies(self, atoms=None):
        atoms = atoms if atoms is not None else self.atoms
        return np.asarray(native.frequencies(*self._args(atoms))["frequencies_cm"], dtype=float)

    def get_gradient(self, atoms=None):
        atoms = atoms if atoms is not None else self.atoms
        return -np.asarray(self.get_forces(atoms), dtype=float)


class MNDO(XNDO):
    """Method-specific MNDO convenience calculator.

    This is equivalent to ``XNDO(method="mndo")`` and remains the exported
    method name for ``from xndo_rs.ase import *``.
    """

    def __init__(self, **kwargs):
        kwargs["method"] = "mndo"
        super().__init__(**kwargs)


class MNDOD(XNDO):
    def __init__(self, **kwargs):
        kwargs["method"] = "mndod"
        super().__init__(**kwargs)


class CNDO2(XNDO):
    def __init__(self, **kwargs):
        kwargs["method"] = "cndo2"
        super().__init__(**kwargs)


class INDO(XNDO):
    def __init__(self, **kwargs):
        kwargs["method"] = "indo"
        super().__init__(**kwargs)


class MINDO3(XNDO):
    def __init__(self, **kwargs):
        kwargs["method"] = "mindo3"
        super().__init__(**kwargs)


class ZINDOS(Calculator):
    """ASE RHF/UHF ground-state and RHF-CIS/UHF-UCIS calculator.

    ``get_uv_vis_spectrum`` returns RHF-CIS or UHF-UCIS stick spectra and
    related excited-state properties. ``forces`` follow the selected reference.
    """

    implemented_properties = ["energy", "forces", "charges", "dipole", "hessian"]

    def __init__(self, charge: int = 0, multiplicity: int = 1,
                 reference: str = "auto", n_states: int = 10,
                 active_occupied=None, active_virtual=None, **kwargs):
        super().__init__(**kwargs)
        self.charge = int(charge)
        self.multiplicity = int(multiplicity)
        self.reference = str(reference)
        self.n_states = int(n_states)
        self.active_occupied = active_occupied
        self.active_virtual = active_virtual

    def _geometry(self, atoms):
        if atoms is None:
            raise RuntimeError("ZINDOS has no Atoms object yet")
        return atoms.get_atomic_numbers(), atoms.get_positions()

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        numbers, positions = self._geometry(self.atoms)
        spectrum = native.uv_vis_spectrum(
            numbers,
            positions,
            charge=self.charge,
            n_states=self.n_states,
            active_occupied=self.active_occupied,
            active_virtual=self.active_virtual,
            multiplicity=self.multiplicity,
            reference=self.reference,
        )
        self.results["energy"] = float(spectrum["ground_energy_ev"])
        self.results["charges"] = np.asarray(spectrum["ground_charges"], dtype=float)
        self.results["dipole"] = (
            np.asarray(spectrum["ground_dipole_debye"], dtype=float)
            * _DEBYE_TO_E_ANGSTROM
        )
        self.results["uv_vis_spectrum"] = spectrum
        if "forces" in properties:
            force_result = native.forces(
                numbers,
                positions,
                charge=self.charge,
                multiplicity=self.multiplicity,
                reference=self.reference,
                method="zindo/s",
            )
            self.results["forces"] = np.asarray(
                force_result["forces_ev_per_angstrom"], dtype=float
            )
        if "hessian" in properties:
            self.results["hessian"] = self.get_hessian(self.atoms)

    def get_gradient(self, atoms=None):
        atoms = atoms if atoms is not None else self.atoms
        return -np.asarray(self.get_forces(atoms), dtype=float)

    def get_hessian(self, atoms=None):
        atoms = atoms if atoms is not None else self.atoms
        numbers, positions = self._geometry(atoms)
        result = native.hessian(
            numbers,
            positions,
            charge=self.charge,
            multiplicity=self.multiplicity,
            reference=self.reference,
            method="zindo/s",
        )
        return np.asarray(result["hessian_hartree_per_bohr2"], dtype=float) * Hartree / (Bohr * Bohr)

    def get_frequencies(self, atoms=None):
        atoms = atoms if atoms is not None else self.atoms
        numbers, positions = self._geometry(atoms)
        result = native.frequencies(
            numbers,
            positions,
            charge=self.charge,
            multiplicity=self.multiplicity,
            reference=self.reference,
            method="zindo/s",
        )
        return np.asarray(result["frequencies_cm"], dtype=float)

    def get_uv_vis_spectrum(self, atoms=None, n_states=None):
        atoms = atoms if atoms is not None else self.atoms
        requested = self.n_states if n_states is None else int(n_states)
        if (
            atoms is self.atoms
            and requested == self.n_states
            and "uv_vis_spectrum" in self.results
        ):
            return self.results["uv_vis_spectrum"]
        numbers, positions = self._geometry(atoms)
        return native.uv_vis_spectrum(
            numbers,
            positions,
            charge=self.charge,
            n_states=requested,
            active_occupied=self.active_occupied,
            active_virtual=self.active_virtual,
            multiplicity=self.multiplicity,
            reference=self.reference,
        )

    def get_excited_properties(self, atoms=None, n_states=None,
                               state_type="singlet"):
        atoms = atoms if atoms is not None else self.atoms
        numbers, positions = self._geometry(atoms)
        return native.excited_properties(
            numbers,
            positions,
            charge=self.charge,
            n_states=self.n_states if n_states is None else int(n_states),
            active_occupied=self.active_occupied,
            active_virtual=self.active_virtual,
            state_type=state_type,
            multiplicity=self.multiplicity,
            reference=self.reference,
        )

    def get_excited_state_gradients(self, atoms=None, n_states=None,
                                    state_type="singlet"):
        atoms = atoms if atoms is not None else self.atoms
        numbers, positions = self._geometry(atoms)
        return native.excited_state_gradients(
            numbers,
            positions,
            charge=self.charge,
            n_states=self.n_states if n_states is None else int(n_states),
            active_occupied=self.active_occupied,
            active_virtual=self.active_virtual,
            state_type=state_type,
            multiplicity=self.multiplicity,
            reference=self.reference,
        )

    def get_excited_state_hessians(self, atoms=None, n_states=None,
                                   state_type="singlet"):
        atoms = atoms if atoms is not None else self.atoms
        numbers, positions = self._geometry(atoms)
        return native.excited_state_hessians(
            numbers,
            positions,
            charge=self.charge,
            n_states=self.n_states if n_states is None else int(n_states),
            active_occupied=self.active_occupied,
            active_virtual=self.active_virtual,
            state_type=state_type,
            multiplicity=self.multiplicity,
            reference=self.reference,
        )


__all__ = ["CNDO2", "INDO", "MNDO", "MNDOD", "MINDO3", "ZINDOS"]
