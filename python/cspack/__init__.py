"""Small NumPy wrapper for selected CSpack routines."""

import importlib

__all__ = [
    "CSpackLibraryError",
    "KernelType",
    "NeutrinoKernel",
    "PhotonKernel",
    "ThermalPhotonKernel",
    "ensure_library",
    "integral_weights",
    "kernel",
    "kernel_bounds",
    "library_path",
    "neutrino_kernel",
    "neutrino_kernel_norm",
    "omega_crit",
    "omega_max",
    "omega_min",
    "planck_occupation",
    "scattering_matrix_explicit",
    "scattering_matrix_kr",
    "thermal_kernel",
    "thermal_neutrino_kernel",
    "thermal_photon_fd_kernel",
]


def __getattr__(name: str):
    if name not in __all__:
        raise AttributeError(f"module 'cspack' has no attribute {name!r}")

    core = importlib.import_module(".core", __name__)
    return getattr(core, name)
