# SPDX-License-Identifier: GPL-3.0-or-later
"""Native xndo-rs API.

Input coordinates are Angstrom. All native ground-state methods expose analytic
gradients and Hessians; every native ground-state method supports RHF/UHF.
ZINDO/S (``INDO/S``) exposes spin-adapted RHF-CIS and spin-orbital UCIS,
including transition/permanent dipoles, oscillator strengths, state charges,
and hole/electron descriptors. Call :func:`api_methods` for the exact accepted
method strings of every public API.
"""
from __future__ import annotations
from typing import Optional, Sequence
import numpy as np
from . import _native


def _as_lists(numbers, positions):
    numbers = [int(z) for z in np.asarray(numbers).reshape(-1)]
    positions = np.asarray(positions, dtype=float).reshape(len(numbers), 3).tolist()
    return numbers, positions


def single_point(numbers: Sequence[int], positions, charge: float = 0.0,
                 multiplicity: int = 1, reference: str = "auto",
                 method: str = "mndo") -> dict:
    """Fixed-geometry calculation.

    ``method`` strings:
      * CNDO/2: ``cndo2``, ``cndo/2``, ``cndo``
      * INDO (MolDS ground-state label): ``indo``
      * MNDO: ``mndo``
      * MNDO/d: ``mndod``, ``mndo/d``, ``mndo-d``
      * MINDO/3: ``mindo3``, ``mindo/3``, ``mindo``
      * ZINDO/S: ``zindo/s``, ``zindos``, ``zindo``, ``indo/s``, ``indos``

    ``reference`` accepts ``auto``, ``rhf`` or ``uhf`` for every native method.
    """
    n, p = _as_lists(numbers, positions)
    return _native.single_point(n, p, float(charge), int(multiplicity),
                                str(reference), str(method))


def orbital_energies(numbers: Sequence[int], positions, charge: float = 0.0,
                     multiplicity: int = 1, reference: str = "auto",
                     method: str = "mndo") -> dict:
    """Orbital energies, occupations and the frontier pair, in eV.

    ``mo_energies_ev`` and ``occupations`` are the restricted spectrum, or the
    alpha one when the reference is unrestricted. ``mo_energies_beta_ev`` and
    ``occupations_beta`` are the beta spectrum, and are ``None`` for a restricted
    reference. Also present: ``n_occ``, ``n_alpha``, ``n_beta``,
    ``homo_lumo_gap_ev``, and the per-channel ``homo_alpha_ev``,
    ``lumo_alpha_ev``, ``homo_beta_ev``, ``lumo_beta_ev``.

    ``homo_ev`` and ``lumo_ev`` are taken over *both* spin channels, so for an
    open-shell doublet the LUMO is usually the beta partner of the singly
    occupied orbital rather than the lowest unoccupied alpha orbital. The two
    channels are eigenvalues of different Fock operators, so neither spectrum
    alone gives the frontier pair. Use ``lumo_alpha_ev`` for the alpha diagram
    on its own. A frontier orbital that does not exist is ``None``.
    """
    n, p = _as_lists(numbers, positions)
    return _native.orbital_energies(n, p, float(charge), int(multiplicity),
                                    str(reference), str(method))


def molden(numbers: Sequence[int], positions, charge: float = 0.0,
           multiplicity: int = 1, reference: str = "auto",
           method: str = "mndo", coefficients: str = "lowdin") -> str:
    """A Molden wavefunction file, as a string.

    The MO coefficients are back-transformed with ``S^(-1/2)`` so that they are
    orthonormal over the Gaussians written in the file. Every engine here assumes
    an orthonormal AO basis, while a Molden file describes real Gaussians, which
    are not orthonormal; a reader handed the raw coefficients would compute a
    density that does not integrate to the electron count. Pass
    ``coefficients="raw"`` for the untransformed coefficients, to compare against
    programs that make that identification -- such a file says so in its title.

    The Slater basis is expanded as STO-6G (Stewart, *J. Chem. Phys.* **52**, 431
    (1970)), and d shells are written in Molden's ``[5D]`` order.
    """
    n, p = _as_lists(numbers, positions)
    return _native.molden(n, p, float(charge), int(multiplicity), str(reference),
                          str(method), str(coefficients))


def third_party_licenses() -> list[dict]:
    """The third-party licence and attribution documents this build embeds.

    One dict per document, with ``path``, ``role`` and the verbatim ``text``.
    These are the notices Apache-2.0 4(c) and GPL-3.0 5(a) require to be carried,
    so they travel with the installed wheel and not only with the source tree.
    """
    return _native.third_party_licenses()


def gradient(numbers: Sequence[int], positions, charge: float = 0.0,
             multiplicity: int = 1, reference: str = "auto",
             method: str = "mndo") -> dict:
    """Analytic Cartesian gradient for any native method."""
    n, p = _as_lists(numbers, positions)
    return _native.gradient(n, p, float(charge), int(multiplicity),
                            str(reference), str(method))


def forces(numbers: Sequence[int], positions, charge: float = 0.0,
           multiplicity: int = 1, reference: str = "auto",
           method: str = "mndo") -> dict:
    """Analytic forces (= -gradient) for any native method."""
    n, p = _as_lists(numbers, positions)
    return _native.forces(n, p, float(charge), int(multiplicity),
                          str(reference), str(method))


def excited_state_gradients(numbers: Sequence[int], positions, charge: float = 0.0,
                            n_states: int = 10, active_occupied: Optional[int] = None,
                            active_virtual: Optional[int] = None,
                            state_type: str = "singlet", multiplicity: int = 1,
                            reference: str = "auto") -> dict:
    """Analytic RHF-CIS or UHF-UCIS excitation and absolute-state gradients."""
    n, p = _as_lists(numbers, positions)
    return _native.excited_state_gradients(
        n, p, float(charge), int(n_states), active_occupied, active_virtual,
        str(state_type), int(multiplicity), str(reference)
    )


def excited_state_hessians(numbers: Sequence[int], positions, charge: float = 0.0,
                           n_states: int = 10, active_occupied: Optional[int] = None,
                           active_virtual: Optional[int] = None,
                           state_type: str = "singlet", multiplicity: int = 1,
                           reference: str = "auto") -> dict:
    """Analytic RHF-CIS or UHF-UCIS excitation and absolute-state Hessians."""
    n, p = _as_lists(numbers, positions)
    return _native.excited_state_hessians(
        n, p, float(charge), int(n_states), active_occupied, active_virtual,
        str(state_type), int(multiplicity), str(reference)
    )


def optimize(numbers: Sequence[int], positions, charge: float = 0.0,
             multiplicity: int = 1, reference: str = "auto",
             method: str = "mndo", max_iter: int = 200,
             gtol: float = 1.0e-3, history: int = 8) -> dict:
    """L-BFGS geometry optimization for any native ground-state method.

    ``gtol`` is the maximum Cartesian gradient-component threshold in eV/Bohr.
    ``max_iter`` and ``history`` control the iteration cap and L-BFGS memory.
    """
    n, p = _as_lists(numbers, positions)
    return _native.optimize(n, p, float(charge), int(multiplicity),
                            str(reference), str(method), int(max_iter),
                            float(gtol), int(history))


def frequencies(numbers: Sequence[int], positions, charge: float = 0.0,
                multiplicity: int = 1, reference: str = "auto",
                method: str = "mndo") -> dict:
    """Harmonic frequencies from an analytic native-method Hessian."""
    n, p = _as_lists(numbers, positions)
    return _native.frequencies(n, p, float(charge), int(multiplicity),
                               str(reference), str(method))


def hessian(numbers: Sequence[int], positions, charge: float = 0.0,
            multiplicity: int = 1, reference: str = "auto",
            method: str = "mndo") -> dict:
    """Analytic Cartesian Hessian in Hartree/Bohr2 for any native method.

    Every native method supports RHF/UHF. Open-shell diatomic state-specific
    Hessians require an explicitly chosen state-averaged or diabatic model.
    """
    n, p = _as_lists(numbers, positions)
    return _native.hessian(n, p, float(charge), int(multiplicity),
                           str(reference), str(method))


def excited_states(numbers: Sequence[int], positions, charge: float = 0.0,
                   n_states: int = 10, active_occupied: Optional[int] = None,
                   active_virtual: Optional[int] = None,
                   method: str = "zindo/s", state_type: str = "singlet",
                   multiplicity: int = 1, reference: str = "auto") -> dict:
    """ZINDO/S vertical excited-state calculation.

    ``method`` accepts ``zindo/s``, ``zindos``, ``zindo``, ``indo/s`` and
    ``indos``. Bare ``indo`` is the distinct ground-state INDO method and is
    deliberately not accepted here. On RHF references, ``state_type`` is
    ``singlet``, ``triplet`` or ``both``; UHF references use UCIS.

    Each state returns excitation and absolute-state energies, spin
    multiplicity and <S2>, wavelength, oscillator strength, transition dipole,
    CIS-state permanent and excited-minus-ground difference dipoles,
    state-specific NDO charges, dominant configurations, hole/electron
    populations and centroids, and the hole-to-electron centroid (CT) distance.

    The numerical engine is the ordinary 1/4-AO s/p ZINDO/S branch.
    Transition-metal 9-AO d-shell ZINDO/S and MSINDO-sCIS/UCIS remain gated
    until their full method-specific Hamiltonians are implemented.
    """
    n, p = _as_lists(numbers, positions)
    return _native.excited_states(n, p, float(charge), int(n_states),
                                  active_occupied, active_virtual,
                                  str(method), str(state_type), int(multiplicity),
                                  str(reference))


def excited_properties(numbers: Sequence[int], positions, charge: float = 0.0,
                       n_states: int = 10, active_occupied: Optional[int] = None,
                       active_virtual: Optional[int] = None,
                       method: str = "zindo/s", state_type: str = "singlet",
                       multiplicity: int = 1, reference: str = "auto") -> dict:
    """Alias of :func:`excited_states`, emphasizing the property-rich output."""
    n, p = _as_lists(numbers, positions)
    return _native.excited_properties(n, p, float(charge), int(n_states),
                                      active_occupied, active_virtual,
                                      str(method), str(state_type), int(multiplicity),
                                      str(reference))


def uv_vis_spectrum(numbers: Sequence[int], positions, charge: float = 0.0,
                    n_states: int = 10, active_occupied: Optional[int] = None,
                    active_virtual: Optional[int] = None,
                    method: str = "zindo/s", multiplicity: int = 1,
                    reference: str = "auto") -> dict:
    """Return a ZINDO/S UV-visible stick spectrum and state properties.

    Each stick includes excitation energy, wavenumber, wavelength, oscillator
    and dipole strengths, transition/permanent/difference dipoles, atomic
    charges, dominant configurations, and hole/electron descriptors. No line
    broadening or solvent shift is applied. RHF uses singlet CIS sticks; UHF
    uses spin-orbital UCIS sticks.
    """
    n, p = _as_lists(numbers, positions)
    return _native.uv_vis_spectrum(n, p, float(charge), int(n_states),
                                   active_occupied, active_virtual, str(method),
                                   int(multiplicity), str(reference))


def available_methods() -> list[dict]:
    """Return every registered method, implementation status, and capabilities."""
    return _native.available_methods()


def api_methods() -> dict[str, list[dict]]:
    """Return the exact method strings accepted by each computational API."""
    return _native.api_methods()


def parameter_datasets() -> list[dict]:
    """Return metadata for every bundled reusable legacy parameter dataset."""
    return _native.parameter_datasets()


def parameter_dataset(name: str) -> str:
    """Return one bundled parameter dataset as its provenance-preserving CSV text."""
    return _native.parameter_dataset(str(name))
