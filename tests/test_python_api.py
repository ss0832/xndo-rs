# SPDX-License-Identifier: GPL-3.0-or-later
import pytest
import numpy as np

xndo_rs = pytest.importorskip("xndo_rs")


def test_top_level_exports():
    for name in ["single_point", "orbital_energies", "molden", "third_party_licenses", "gradient", "forces", "optimize", "frequencies", "hessian", "excited_states", "excited_properties", "uv_vis_spectrum", "excited_state_gradients", "excited_state_hessians", "available_methods", "api_methods", "parameter_datasets", "parameter_dataset"]:
        assert name in xndo_rs.__all__
        assert callable(getattr(xndo_rs, name))


def test_every_native_function_is_wired_through_all_three_layers():
    """The extension, the shim and the package must export the same set.

    There are three layers -- the pyo3 module `_native`, the `native.py` shim
    that adapts arrays and documents the arguments, and the package `__init__`
    that re-exports -- and adding a function to the first without the other two
    is invisible to a hand-written list of names. It happened: `_native`
    registered `orbital_energies`, `molden` and `third_party_licenses`, the
    package `__init__` re-exported all three, the shim defined none of them, and
    `import xndo_rs` raised AttributeError at line 6 for anyone who installed
    the wheel. Deriving the expectation from `_native` is what makes the next
    one fail here instead.
    """
    from xndo_rs import native
    exported = {
        name for name in dir(native._native)
        if not name.startswith("_") and callable(getattr(native._native, name))
    }
    assert exported, "no functions found in the compiled module"
    missing_shim = sorted(n for n in exported if not callable(getattr(native, n, None)))
    assert not missing_shim, f"registered in _native but absent from native.py: {missing_shim}"
    missing_pkg = sorted(n for n in exported if n not in xndo_rs.__all__)
    assert not missing_pkg, f"in native.py but not re-exported by the package: {missing_pkg}"


def test_no_exported_name_is_dangling():
    for name in xndo_rs.__all__:
        assert hasattr(xndo_rs, name), f"__all__ names {name}, which does not resolve"


def test_ase_classes_if_available():
    pytest.importorskip("ase")
    import xndo_rs.ase as xndo_ase
    from xndo_rs.ase import CNDO2, INDO, MINDO3, XNDO, MNDO, MNDOD, ZINDOS
    assert "XNDO" not in xndo_ase.__all__
    assert "MNDO" in xndo_ase.__all__
    assert XNDO(method="mndo").method == "mndo"
    assert MNDO().method == "mndo"
    assert MNDOD().method == "mndod"
    assert XNDO(method="mndo/d").method == "mndo/d"
    assert XNDO(method="mndo-d").method == "mndo-d"
    assert XNDO(method="zindo/s").method == "zindo/s"
    assert CNDO2().method == "cndo2"
    assert INDO().method == "indo"
    assert MINDO3().method == "mindo3"
    assert ZINDOS(n_states=3).n_states == 3


def test_ase_runtime_energy_forces_hessian_and_frequencies():
    pytest.importorskip("ase")
    from ase import Atoms
    from xndo_rs.ase import MNDO

    atoms = Atoms(
        numbers=[8, 1, 1],
        positions=[[0.0, 0.0, 0.0], [0.9584, 0.0, 0.0], [-0.2400, 0.9278, 0.0]],
    )
    atoms.calc = MNDO()
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    hessian = atoms.calc.get_hessian(atoms)
    property_hessian = atoms.calc.get_property("hessian", atoms)
    frequencies = atoms.calc.get_frequencies(atoms)

    native = xndo_rs.single_point(atoms.numbers, atoms.positions, method="mndo")
    assert energy == pytest.approx(native["energy_ev"], abs=1.0e-10)
    assert forces.shape == (3, 3) and np.isfinite(forces).all()
    assert hessian.shape == (9, 9) and np.isfinite(hessian).all()
    assert property_hessian == pytest.approx(hessian)
    assert frequencies.shape == (9,) and np.isfinite(frequencies).all()


def test_ase_zindos_uv_vis_runtime():
    pytest.importorskip("ase")
    from ase import Atoms
    from xndo_rs.ase import ZINDOS

    atoms = Atoms(
        numbers=[8, 1, 1],
        positions=[[0.0, 0.0, 0.0], [0.9584, 0.0, 0.0], [-0.2400, 0.9278, 0.0]],
    )
    atoms.calc = ZINDOS(n_states=3)
    assert np.isfinite(atoms.get_potential_energy())
    assert atoms.get_charges().shape == (3,)
    assert atoms.get_dipole_moment().shape == (3,)
    spectrum = atoms.calc.get_uv_vis_spectrum(atoms)
    assert spectrum["state_type"] == "singlet"
    assert len(spectrum["states"]) == 3
    for state in spectrum["states"]:
        assert state["wavelength_nm"] > 0.0
        assert state["oscillator_strength"] >= 0.0


def test_ase_zindos_uhf_ucis_runtime():
    pytest.importorskip("ase")
    from ase import Atoms
    from xndo_rs.ase import ZINDOS

    atoms = Atoms(
        numbers=[6, 1, 1, 1],
        positions=[
            [0.07, -0.02, 0.05],
            [0.98, 0.10, -0.06],
            [-0.38, 0.82, 0.12],
            [-0.34, -0.77, -0.18],
        ],
    )
    atoms.calc = ZINDOS(
        multiplicity=2,
        reference="uhf",
        n_states=1,
        active_occupied=1,
        active_virtual=1,
    )
    assert np.isfinite(atoms.get_potential_energy())
    assert np.isfinite(atoms.get_forces()).all()
    assert np.isfinite(atoms.calc.get_hessian(atoms)).all()
    spectrum = atoms.calc.get_uv_vis_spectrum(atoms)
    gradient = atoms.calc.get_excited_state_gradients(atoms)
    assert spectrum["reference"] == "UHF"
    assert spectrum["states"][0]["spin"] == "unrestricted"
    assert gradient["states"][0]["spin"] == "unrestricted"


def test_method_registry_surface():
    methods = {row["name"]: row for row in xndo_rs.available_methods()}
    assert methods["MNDO"]["status"] == "native"
    assert methods["MINDO/3"]["status"] == "native"
    assert methods["CNDO/2"]["status"] == "native"
    assert methods["INDO"]["status"] == "native"
    assert methods["INDO"]["accepted_strings"] == ["indo"]
    assert methods["CNDO/3"]["status"] == "noncanonical"


def test_parameter_dataset_surface():
    catalog = {row["dataset_id"]: row for row in xndo_rs.parameter_datasets()}
    assert "mopac7_mindo3" in catalog
    assert catalog["molds_cndo2_indo"]["upstream_license"] == "GPL-3.0-or-later"
    text = xndo_rs.parameter_dataset("mindo3_pair")
    assert "15,17,0.277322,1.543720" in text
    assert "17,17,0.258969,1.792125" in text


def test_api_method_strings_are_discoverable():
    apis = xndo_rs.api_methods()
    assert {row["method"] for row in apis["single_point"]} >= {"CNDO/2", "INDO", "MNDO", "MNDO/d", "MINDO/3", "ZINDO/S"}
    assert {row["method"] for row in apis["optimize"]} == {"CNDO/2", "INDO", "MNDO", "MNDO/d", "MINDO/3", "ZINDO/S"}
    z = next(row for row in apis["excited_properties"] if row["method"] == "ZINDO/S")
    for alias in ["zindo/s", "zindos", "zindo", "indo/s", "indos"]:
        assert alias in z["accepted_strings"]
    assert "indo" not in z["accepted_strings"]
    assert z["reference_argument"] is True
    assert z["reference_strings"] == ["auto", "rhf", "uhf"]
    assert z["fixed_reference"] is None
    assert z["state_type_strings"] == ["singlet", "triplet", "both", "unrestricted"]
    uv = next(row for row in apis["uv_vis_spectrum"] if row["method"] == "ZINDO/S")
    assert uv["state_type_strings"] == ["singlet", "unrestricted"]
    cndo = next(row for row in apis["single_point"] if row["method"] == "CNDO/2")
    indo = next(row for row in apis["single_point"] if row["method"] == "INDO")
    assert cndo["accepted_strings"] == ["cndo2", "cndo/2", "cndo"]
    assert indo["accepted_strings"] == ["indo"]
    assert cndo["uhf"] is True and indo["uhf"] is True
    assert cndo["reference_strings"] == ["auto", "rhf", "uhf"]
    assert indo["reference_strings"] == ["auto", "rhf", "uhf"]
    zsp = next(row for row in apis["single_point"] if row["method"] == "ZINDO/S")
    assert zsp["reference_argument"] is True
    assert zsp["reference_strings"] == ["auto", "rhf", "uhf"]
    assert zsp["fixed_reference"] is None
    mndo = next(row for row in apis["hessian"] if row["method"] == "MNDO")
    assert mndo["uhf"] is True


@pytest.mark.parametrize("method", ["cndo2", "indo", "mndo", "mndod", "mindo3", "zindo/s"])
def test_optimize_runs_for_every_native_method(method):
    numbers = [8, 1, 1]
    positions = [[0.0, 0.0, 0.0], [1.05, 0.0, 0.0], [-0.30, 1.02, 0.0]]
    initial = xndo_rs.single_point(numbers, positions, method=method)
    result = xndo_rs.optimize(numbers, positions, method=method, gtol=2.0e-3)
    expected_name = next(
        row["method"] for row in xndo_rs.api_methods()["optimize"]
        if method in row["accepted_strings"]
    )
    assert result["method"] == expected_name
    assert result["converged"] is True
    assert result["energy_ev"] < initial["energy_ev"]
    assert np.asarray(result["positions_angstrom"]).shape == (3, 3)


def test_optimize_rejects_invalid_convergence_settings():
    numbers = [1, 1]
    positions = [[0.0, 0.0, 0.0], [0.8, 0.0, 0.0]]
    with pytest.raises(ValueError, match="gtol"):
        xndo_rs.optimize(numbers, positions, gtol=0.0)
    with pytest.raises(ValueError, match="max_iter"):
        xndo_rs.optimize(numbers, positions, max_iter=0)
    with pytest.raises(ValueError, match="history"):
        xndo_rs.optimize(numbers, positions, history=0)
