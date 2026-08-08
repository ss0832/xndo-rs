# SPDX-License-Identifier: GPL-3.0-or-later
import pytest

xndo_rs = pytest.importorskip("xndo_rs")

WATER_Z = [8, 1, 1]
WATER_XYZ = [[0.0, 0.0, 0.0], [0.9584, 0.0, 0.0], [-0.2400, 0.9278, 0.0]]


def test_mndo_single_point_surface():
    r = xndo_rs.single_point(WATER_Z, WATER_XYZ, method="mndo")
    assert r["method"] == "MNDO"
    assert r["converged"]
    assert "energy_hartree" in r and "heat_of_formation_kcal" in r


def test_mndod_gradient_surface():
    r = xndo_rs.gradient(WATER_Z, WATER_XYZ, method="mndod")
    assert r["method"] == "MNDO/d"
    assert len(r["gradient_hartree_per_bohr"]) == 3


def test_zindo_spectrum_surface():
    r = xndo_rs.excited_states(WATER_Z, WATER_XYZ, n_states=3)
    assert r["method"] == "ZINDO/S"
    assert "states" in r


def test_explicit_uv_vis_spectrum_surface():
    r = xndo_rs.uv_vis_spectrum(WATER_Z, WATER_XYZ, n_states=3)
    assert r["method"] == "ZINDO/S"
    assert r["state_type"] == "singlet"
    assert len(r["states"]) == 3
    assert all(state["wavelength_nm"] > 0.0 for state in r["states"])
    assert all(state["oscillator_strength"] >= 0.0 for state in r["states"])


def test_all_native_ground_derivative_surfaces():
    for method in ["cndo2", "indo", "mindo3", "zindo/s"]:
        gradient = xndo_rs.gradient(WATER_Z, WATER_XYZ, method=method)
        hessian = xndo_rs.hessian(WATER_Z, WATER_XYZ, method=method)
        assert len(gradient["gradient_hartree_per_bohr"]) == 3
        assert len(hessian["hessian_hartree_per_bohr2"]) == 9


def test_zindo_cis_derivative_surfaces():
    gradients = xndo_rs.excited_state_gradients(
        WATER_Z, WATER_XYZ, n_states=1, state_type="both"
    )
    hessians = xndo_rs.excited_state_hessians(
        WATER_Z, WATER_XYZ, n_states=1, state_type="both"
    )
    assert {state["spin"] for state in gradients["states"]} == {"singlet", "triplet"}
    assert {state["spin"] for state in hessians["states"]} == {"singlet", "triplet"}
    for state in gradients["states"]:
        assert len(state["state_gradient_hartree_per_bohr"]) == 3
    for state in hessians["states"]:
        assert len(state["state_hessian_hartree_per_bohr2"]) == 9


def test_new_uhf_ground_derivative_surfaces():
    numbers = [6, 1, 1, 1]
    positions = [
        [0.0, 0.0, 0.05],
        [1.07, 0.0, -0.02],
        [-0.535, 0.927, 0.01],
        [-0.515, -0.942, -0.04],
    ]
    for method in ["cndo2", "indo", "mindo3", "zindo/s"]:
        gradient = xndo_rs.gradient(
            numbers, positions, multiplicity=2, reference="uhf", method=method
        )
        hessian = xndo_rs.hessian(
            numbers, positions, multiplicity=2, reference="uhf", method=method
        )
        assert gradient["unrestricted"] is True
        assert len(gradient["gradient_hartree_per_bohr"]) == 4
        assert len(hessian["hessian_hartree_per_bohr2"]) == 12


def test_zindo_uhf_ucis_native_surfaces():
    numbers = [6, 1, 1, 1]
    positions = [
        [0.07, -0.02, 0.05],
        [0.98, 0.10, -0.06],
        [-0.38, 0.82, 0.12],
        [-0.34, -0.77, -0.18],
    ]
    kwargs = {
        "multiplicity": 2,
        "reference": "uhf",
        "n_states": 1,
        "active_occupied": 1,
        "active_virtual": 1,
    }
    ground = xndo_rs.single_point(
        numbers, positions, multiplicity=2, reference="uhf", method="zindo/s"
    )
    spectrum = xndo_rs.excited_states(numbers, positions, **kwargs)
    gradient = xndo_rs.excited_state_gradients(numbers, positions, **kwargs)
    hessian = xndo_rs.excited_state_hessians(numbers, positions, **kwargs)
    assert ground["unrestricted"] is True
    assert ground["spin_density"] is not None
    assert spectrum["reference"] == "UHF"
    assert spectrum["state_type"] == "unrestricted"
    assert spectrum["states"][0]["spin"] == "unrestricted"
    assert gradient["states"][0]["spin"] == "unrestricted"
    assert hessian["states"][0]["spin"] == "unrestricted"
    assert len(gradient["states"][0]["state_gradient_hartree_per_bohr"]) == 4
    assert len(hessian["states"][0]["state_hessian_hartree_per_bohr2"]) == 12


def test_oh_radical_is_retained_with_explicit_degenerate_derivative_boundary():
    numbers = [8, 1]
    positions = [[0.11, -0.07, 0.03], [1.01, 0.12, -0.08]]
    kwargs = {"multiplicity": 2, "reference": "uhf"}
    gradient = xndo_rs.gradient(numbers, positions, method="zindo/s", **kwargs)
    spectrum = xndo_rs.uv_vis_spectrum(numbers, positions, n_states=1, **kwargs)
    assert gradient["unrestricted"] is True
    assert len(gradient["gradient_hartree_per_bohr"]) == 2
    assert spectrum["state_type"] == "unrestricted"
    assert len(spectrum["states"]) == 1
    with pytest.raises(ValueError, match="not unique"):
        xndo_rs.hessian(numbers, positions, method="zindo/s", **kwargs)
    with pytest.raises(ValueError, match="not unique"):
        xndo_rs.excited_state_gradients(numbers, positions, n_states=1, **kwargs)
    with pytest.raises(ValueError, match="not unique"):
        xndo_rs.excited_state_hessians(numbers, positions, n_states=1, **kwargs)


def test_mindo3_single_point_surface():
    r = xndo_rs.single_point(WATER_Z, WATER_XYZ, method="mindo3")
    assert r["method"] == "MINDO/3"
    assert "heat_of_formation_kcal" in r


def test_cndo2_single_point_surface():
    r = xndo_rs.single_point(WATER_Z, WATER_XYZ, method="cndo2")
    assert r["method"] == "CNDO/2"
    assert r["converged"]
    assert r["unrestricted"] is False
    assert len(r["charges"]) == 3


def test_bare_indo_is_ground_state_indo_not_zindo_s():
    r = xndo_rs.single_point(WATER_Z, WATER_XYZ, method="indo")
    assert r["method"] == "INDO"
    assert r["converged"]
    assert "states" not in r
    z = xndo_rs.single_point(WATER_Z, WATER_XYZ, method="indo/s")
    assert z["method"] == "ZINDO/S"

def test_noncanonical_cndo3_is_explicitly_rejected():
    with pytest.raises(ValueError, match="canonical"):
        xndo_rs.single_point(WATER_Z, WATER_XYZ, method="cndo3")


def test_mndo_uhf_surface():
    z = [6, 1, 1, 1]
    xyz = [[0.0, 0.0, 0.05], [1.07, 0.0, -0.02], [-0.535, 0.927, 0.01], [-0.515, -0.942, -0.04]]
    r = xndo_rs.single_point(z, xyz, multiplicity=2, reference="uhf", method="mndo")
    assert r["unrestricted"] is True
    assert r["spin_density"] is not None


def test_zindo_excited_properties_surface_is_property_rich():
    r = xndo_rs.excited_properties(WATER_Z, WATER_XYZ, n_states=2, state_type="both")
    assert r["method"] == "ZINDO/S"
    spins = {state["spin"] for state in r["states"]}
    assert "singlet" in spins
    assert "triplet" in spins
    for state in r["states"]:
        for key in [
            "energy_ev", "energy_cm1", "wavelength_nm", "oscillator_strength",
            "state_total_energy_ev", "state_total_energy_hartree",
            "spin_multiplicity", "s2_expectation",
            "transition_dipole_au", "transition_dipole_debye",
            "permanent_dipole_debye", "difference_dipole_debye",
            "difference_dipole_magnitude_debye", "charges", "hole_population",
            "electron_population", "hole_centroid_angstrom",
            "electron_centroid_angstrom", "charge_transfer_distance_angstrom",
            "dominant_configurations",
        ]:
            assert key in state
        assert abs(sum(state["hole_population"]) - 1.0) < 1.0e-6
        assert abs(sum(state["electron_population"]) - 1.0) < 1.0e-6
        if state["spin"] == "triplet":
            assert state["oscillator_strength"] == 0.0
            assert state["spin_multiplicity"] == 3
            assert state["s2_expectation"] == 2.0
        else:
            assert state["spin_multiplicity"] == 1
            assert state["s2_expectation"] == 0.0
