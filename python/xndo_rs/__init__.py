# SPDX-License-Identifier: GPL-3.0-or-later
"""xndo-rs Python API for legacy NDO-family semiempirical methods."""
from . import native

single_point = native.single_point
gradient = native.gradient
forces = native.forces
optimize = native.optimize
frequencies = native.frequencies
hessian = native.hessian
excited_states = native.excited_states
excited_properties = native.excited_properties
uv_vis_spectrum = native.uv_vis_spectrum
excited_state_gradients = native.excited_state_gradients
excited_state_hessians = native.excited_state_hessians
available_methods = native.available_methods
api_methods = native.api_methods
parameter_datasets = native.parameter_datasets
parameter_dataset = native.parameter_dataset

__all__ = [
    "native", "single_point", "gradient", "forces", "optimize",
    "frequencies", "hessian", "excited_states", "excited_properties",
    "uv_vis_spectrum", "excited_state_gradients", "excited_state_hessians",
    "available_methods", "api_methods", "parameter_datasets", "parameter_dataset",
]
__version__ = "0.2.4"
