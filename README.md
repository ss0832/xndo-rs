# xndo-rs 

`xndo-rs` is a Rust library with Python-native, ASE, and CLI
interfaces for selected NDO-family semiempirical methods.

## Implemented methods

| Method | Ground state | Derivatives | Excited states | Parameter source |
| --- | --- | --- | --- | --- |
| CNDO/2 | RHF/UHF | analytic gradient, Hessian, optimization, frequencies | no | MolDS 0.3.1 |
| INDO | RHF/UHF | analytic gradient, Hessian, optimization, frequencies | no | MolDS 0.3.1 |
| MNDO | RHF/UHF | analytic gradient, Hessian, optimization, frequencies | no | OpenMOPAC 23.2.5 |
| MNDO/d | RHF/UHF | analytic gradient, Hessian, optimization, frequencies, including open-shell d AOs | no | OpenMOPAC 23.2.5 |
| MINDO/3 | RHF/UHF | analytic gradient, Hessian, optimization, frequencies | no | MOPAC7 |
| ZINDO/S | RHF/UHF | analytic ground optimization/derivatives and state gradients/Hessians | RHF singlet/triplet CIS; UHF spin-orbital UCIS; UV-vis properties | OpenMOPAC 23.2.5 |

Every method also reports orbital energies, occupations and the HOMO/LUMO pair,
and can write a Molden wavefunction file. Molden coefficients are
back-transformed with `S^(-1/2)`: the engines assume an orthonormal AO basis and
a Molden file does not describe one, so the raw coefficients would give a
reader a density that does not integrate to the electron count.

Each method is checked against the program its own parameters came from --
OpenMOPAC 23.2.5, MolDS 0.3.1, MOPAC7 1.15 -- over 2415 independent scalar
comparisons with no known disagreements. `VALIDATION_STATUS.md` has the table
and `tests/data/ORACLE_NOTES.md` the detail.

**No correlated method is implemented**: no MP2, no coupled cluster, no CI, no
multireference. `docs/scope.md` lists what is out, and
`docs/v0.3.0-delivery.md` records what this release delivers against what was
planned, together with its known problems.

CNDO/1, INDO/1, INDO/2, ZINDO/1, ZINDO/2, MINDO/1, MINDO/2,
SINDO1, and MSINDO are registered but deliberately fail at execution because
no complete, compatible, provenance-audited parameter and Hamiltonian set is
bundled. CNDO/3, INDO/3, and ZINDO/3 are non-canonical compatibility labels.

For open-shell diatomics such as OH, the UHF ground-state energy and analytic
gradient and the vertical UCIS spectrum are available and tested. A
state-specific UHF Hessian or UCIS derivative is not unique on the degenerate
electronic surface without a state-averaged or diabatic definition, so those
requests return an explicit error instead of a branch-dependent value.
The same safeguard applies when requested CIS/UCIS roots are separated by less
than 0.0001 eV; vertical spectra remain available, while state-specific
derivatives require a state-averaged or diabatic treatment.

## Build and test

Install the published Python package with optional ASE support:

```text
pip install "xndo-rs-python[ase]"
```

Build local PyPI artifacts with:

```text
cargo check --all-targets --all-features
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
maturin build --release --features python
maturin sdist
twine check target/wheels/*
pip install --force-reinstall --no-deps target/wheels/*.whl
python -m pytest tests/test_python.py tests/test_python_api.py
```

Install the wheel before running the Python tests. They open with
`pytest.importorskip("xndo_rs")`, so against a source tree with nothing
installed they all skip and report success without having run.

A release build records the directories it was built from -- the project, the
cargo registry and the rustup toolchains -- inside the binary, which for a
published artifact means the maintainer's home directory. `[profile.release]`
sets `strip = "symbols"`; supply the rest at build time. `CARGO_HOME` and
`RUSTUP_HOME` are the canonical variables for the other two; set them if your
environment has not already.

```text
RUSTFLAGS="--remap-path-prefix=$PWD=. --remap-path-prefix=$CARGO_HOME=/cargo --remap-path-prefix=$RUSTUP_HOME=/rustup" maturin build --release --features python
```

Then check, rather than assume:

```text
python tools/privacy_scan.py --artifacts target/wheels/*.whl target/wheels/*.tar.gz target/release/xndo_rs_cli*
```

The scanner reads inside archives and binaries. It permits an e-mail address
only where the same address appears in this tree's licence and notice documents,
because those are the attributions Apache-2.0 4(c) and the GPL require to be
carried and `src/licenses.rs` embeds them in every binary. One known finding it
will report is maturin's generated `.dist-info/sboms/*.cyclonedx.json`, which
records absolute source paths; remove it, or do not publish that wheel.

The release profile uses fat LTO and one code-generation unit. Embedded tables
are parsed once per process. The ZDO Hessian engines solve all Cartesian CPHF
right-hand sides with one factorization. Memory guards are available through
`NddoOptions`; `XNDO_MAX_PAIR_CACHE_MB` can override the two-electron
pair-cache limit.

## Rust API

```rust
use xndo_rs::{run_gradient, run_hessian, Method, Molecule, NddoOptions};

let xyz = "3\nwater\nO 0 0 0\nH 0.9584 0 0\nH -0.2400 0.9278 0\n";
let molecule = Molecule::from_xyz_str(xyz, 0.0)?;
let gradient = run_gradient(&molecule, Method::Cndo2, &NddoOptions::default())?;
let hessian = run_hessian(&molecule, Method::Cndo2, &NddoOptions::default())?;
println!("{} eV; {} Cartesian coordinates", gradient.energy_ev, hessian.hessian.rows);
# Ok::<(), xndo_rs::XndoError>(())
```

Use `run_method` for the unified fixed-geometry dispatcher and
`methods_for_api` for programmatic capability discovery.

## Python and ASE

```python
import xndo_rs

z = [8, 1, 1]
xyz = [[0.0, 0.0, 0.0], [0.9584, 0.0, 0.0], [-0.2400, 0.9278, 0.0]]
result = xndo_rs.single_point(z, xyz, method="mndo")
gradient = xndo_rs.gradient(z, xyz, method="cndo2")
hessian = xndo_rs.hessian(z, xyz, method="mindo3")
optimized = xndo_rs.optimize(z, xyz, method="indo", gtol=2.0e-3)
states = xndo_rs.excited_states(z, xyz, method="zindo/s", n_states=5)
uv_vis = xndo_rs.uv_vis_spectrum(z, xyz, n_states=10)
state_gradients = xndo_rs.excited_state_gradients(z, xyz, n_states=5)
state_hessians = xndo_rs.excited_state_hessians(z, xyz, n_states=5)
```

```python
from ase import Atoms
from xndo_rs.ase import XNDO

atoms = Atoms(numbers=z, positions=xyz)
atoms.calc = XNDO(method="mndo")
print(atoms.get_potential_energy())
print(atoms.get_forces())
```

Python coordinates and ASE positions are Angstrom. Returned dictionary keys
carry units in their names. ASE uses eV and Angstrom conventions.

ref.: A. H. Larsen
et al., "The atomic simulation environment - a Python library for working with
atoms," *J. Phys.: Condens. Matter* **29**, 273002 (2017),
[doi:10.1088/1361-648X/aa680e](https://doi.org/10.1088/1361-648X/aa680e).

## CLI

```text
xndo_rs_cli methods
xndo_rs_cli energy examples/water.xyz --method mndo
xndo_rs_cli gradient examples/water.xyz --method mndod
xndo_rs_cli hessian examples/water.xyz --method cndo2
xndo_rs_cli optimize examples/water.xyz --method mindo3 --opt-gtol 2e-3
xndo_rs_cli uv-vis examples/water.xyz --method zindo/s --states 10
xndo_rs_cli excited-gradient examples/water.xyz --method zindo/s --states 5
xndo_rs_cli excited-hessian examples/water.xyz --method zindo/s --states 5
```

See `docs/methods.md`, `docs/parameter-provenance.md`,
`docs/parameter-search.md`, and `docs/oracles.md` for exact scope and evidence.

## License

The library is GPL-3.0-or-later. Retained upstream notices and data licenses
are recorded in `THIRD_PARTY_NOTICES.md` and `third_party/`.
