from __future__ import annotations

import ctypes
import os
import sys
from enum import Enum
from pathlib import Path
from typing import Callable, Optional, Union

import numpy as np

ArrayLike = Union[float, list[float], tuple[float, ...], np.ndarray]


class CSpackLibraryError(RuntimeError):
    """Raised when the CSpack shared library cannot be loaded."""


class PhotonKernel(str, Enum):
    EXACT = "exact"
    RECOIL = "recoil"
    DOPPLER = "doppler"
    UR = "ur"


class ThermalPhotonKernel(str, Enum):
    EXACT = "exact"
    EXACT_SS_C = "exact+SS_C"
    RECOIL = "recoil"
    DOPPLER = "doppler"
    UR = "ur"
    SS_K = "SS_K"
    SS_C = "SS_C"


class NeutrinoKernel(str, Enum):
    NUE_E = "nue_e"
    NUE_P = "nue_p"
    NUE_EP = "nue_ep"
    NUMU_E = "numu_e"
    NUMU_P = "numu_p"
    NUTAU_E = "nutau_e"
    NUTAU_P = "nutau_p"


class KernelType(str, Enum):
    EXACT = "exact"
    SS_K = "SS_K"
    SS_C = "SS_C"


_LIB: Optional[ctypes.CDLL] = None
_DOUBLE_PTR = ctypes.POINTER(ctypes.c_double)
_DOUBLE_PTR_PTR = ctypes.POINTER(_DOUBLE_PTR)
_NDARRAY_1D = np.ctypeslib.ndpointer(
    dtype=np.float64,
    ndim=1,
    flags=("C_CONTIGUOUS", "ALIGNED"),
)


def _library_filename() -> str:
    if sys.platform == "darwin":
        return "libcspack_py.dylib"
    if os.name == "nt":
        return "cspack_py.dll"
    return "libcspack_py.so"


def _package_library_path() -> Path:
    return Path(__file__).resolve().parent / "lib" / _library_filename()


def library_path() -> Optional[Path]:
    """Return the shared-library path that will be loaded, if it exists."""
    override = os.environ.get("CSPACK_LIBRARY")
    if override:
        path = Path(override).expanduser()
        return path if path.exists() else None

    path = _package_library_path()
    return path if path.exists() else None


def ensure_library(force: bool = False, **build_options: object) -> Path:
    """Build the local CSpack shared library if needed and return its path.

    Extra keyword arguments are forwarded to ``cspack.build.build_shared_library``.
    """
    global _LIB

    path = library_path()
    if path is not None and not force:
        return path

    from .build import build_shared_library

    path = build_shared_library(**build_options)
    _LIB = None
    return path


def _load_library() -> ctypes.CDLL:
    global _LIB

    if _LIB is not None:
        return _LIB

    path = library_path()
    if path is None:
        raise CSpackLibraryError(
            "Could not find the CSpack shared library. Run "
            "`python -m cspack.build` from the python directory, or set "
            "CSPACK_LIBRARY to an existing shared library path."
        )

    try:
        lib = ctypes.CDLL(str(path))
    except OSError as exc:
        raise CSpackLibraryError(f"Failed to load CSpack shared library at {path}: {exc}") from exc

    lib.nbb_func_C.argtypes = [ctypes.c_double]
    lib.nbb_func_C.restype = ctypes.c_double

    lib.Integral_weights.argtypes = [_NDARRAY_1D, _NDARRAY_1D, ctypes.c_int]
    lib.Integral_weights.restype = None

    lib.compute_scattering_matrix_explicit.argtypes = [
        _NDARRAY_1D,
        _NDARRAY_1D,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_double,
        _DOUBLE_PTR_PTR,
    ]
    lib.compute_scattering_matrix_explicit.restype = None

    lib.compute_scattering_matrix_KR.argtypes = [
        _NDARRAY_1D,
        _NDARRAY_1D,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        _DOUBLE_PTR_PTR,
    ]
    lib.compute_scattering_matrix_KR.restype = None

    lib.cspack_omega_min.argtypes = [
        ctypes.c_double,
        ctypes.c_double,
    ]
    lib.cspack_omega_min.restype = ctypes.c_double

    lib.cspack_omega_max.argtypes = [
        ctypes.c_double,
        ctypes.c_double,
    ]
    lib.cspack_omega_max.restype = ctypes.c_double

    lib.cspack_omega_crit.argtypes = [
        ctypes.c_double,
        ctypes.c_double,
    ]
    lib.cspack_omega_crit.restype = ctypes.c_double

    lib.cspack_kernel.argtypes = [
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
    ]
    lib.cspack_kernel.restype = ctypes.c_double

    lib.cspack_thermal_kernel.argtypes = [
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
    ]
    lib.cspack_thermal_kernel.restype = ctypes.c_double

    lib.cspack_neutrino_kernel.argtypes = [
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
    ]
    lib.cspack_neutrino_kernel.restype = ctypes.c_double

    lib.cspack_thermal_photon_fd_kernel.argtypes = [
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_int,
    ]
    lib.cspack_thermal_photon_fd_kernel.restype = ctypes.c_double

    lib.cspack_thermal_neutrino_kernel.argtypes = [
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_int,
    ]
    lib.cspack_thermal_neutrino_kernel.restype = ctypes.c_double

    lib.cspack_neutrino_kernel_norm.argtypes = [ctypes.c_int]
    lib.cspack_neutrino_kernel_norm.restype = ctypes.c_double

    _LIB = lib
    return lib


def _as_grid(x: ArrayLike) -> np.ndarray:
    arr = np.ascontiguousarray(np.asarray(x, dtype=np.float64))
    if arr.ndim != 1:
        raise ValueError("x grid must be one-dimensional")
    if arr.size < 5:
        raise ValueError("x grid must contain at least five points")
    if not np.all(np.isfinite(arr)):
        raise ValueError("x grid contains non-finite values")
    if not np.all(np.diff(arr) > 0.0):
        raise ValueError("x grid must be strictly increasing")
    return arr


def _as_weights(weights: Optional[ArrayLike], x: np.ndarray) -> np.ndarray:
    if weights is None:
        return integral_weights(x)

    arr = np.ascontiguousarray(np.asarray(weights, dtype=np.float64))
    if arr.ndim != 1 or arr.size != x.size:
        raise ValueError("weights must be a one-dimensional array with the same length as x")
    if not np.all(np.isfinite(arr)):
        raise ValueError("weights contain non-finite values")
    return arr


def _positive_float(value: float, name: str) -> float:
    out = float(value)
    if not np.isfinite(out) or out <= 0.0:
        raise ValueError(f"{name} must be a positive finite value")
    return out


def _enum_code(value: Union[str, Enum], mapping: dict[str, int], label: str) -> int:
    text = value.value if isinstance(value, Enum) else str(value)
    normalized = text.strip().lower().replace("-", "_")
    try:
        return mapping[normalized]
    except KeyError as exc:
        names = ", ".join(repr(name) for name in mapping)
        raise ValueError(f"{label} must be one of {names}") from exc


def _matrix_kernel_code(kernel: Union[str, KernelType]) -> int:
    return _enum_code(
        kernel,
        {
            "exact": 0,
            "ss_k": 1,
            "ss_c": 2,
        },
        "kernel",
    )


def _photon_kernel_code(kernel: Union[str, PhotonKernel]) -> int:
    return _enum_code(
        kernel,
        {
            "exact": 0,
            "recoil": 1,
            "doppler": 2,
            "ur": 3,
        },
        "kernel",
    )


def _thermal_photon_kernel_code(kernel: Union[str, ThermalPhotonKernel]) -> int:
    return _enum_code(
        kernel,
        {
            "exact": 0,
            "exact+ss_c": 1,
            "exact_ss_c": 1,
            "recoil": 2,
            "doppler": 3,
            "ur": 4,
            "ss_k": 5,
            "ss_c": 6,
        },
        "kernel",
    )


def _neutrino_kernel_code(kernel: Union[str, NeutrinoKernel]) -> int:
    return _enum_code(
        kernel,
        {
            "nue_e": 0,
            "nue_p": 1,
            "nue_ep": 2,
            "numu_e": 3,
            "numu_p": 4,
            "nutau_e": 5,
            "nutau_p": 6,
        },
        "kernel",
    )


def _matrix_row_pointers(matrix: np.ndarray) -> ctypes.Array[_DOUBLE_PTR]:
    row_type = _DOUBLE_PTR * matrix.shape[0]
    return row_type(*[matrix[i].ctypes.data_as(_DOUBLE_PTR) for i in range(matrix.shape[0])])


def _finite_array(value: ArrayLike, name: str) -> np.ndarray:
    arr = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _positive_array(value: ArrayLike, name: str) -> np.ndarray:
    arr = _finite_array(value, name)
    if not np.all(arr > 0.0):
        raise ValueError(f"{name} must contain positive values")
    return arr


def _broadcast_kernel_call(
    function: Callable[..., float],
    kernel_code: int,
    *values: ArrayLike,
    converters: Optional[tuple[Callable[[float], object], ...]] = None,
) -> Union[float, np.ndarray]:
    if converters is None:
        converters = tuple(float for _ in values)
    if len(converters) != len(values):
        raise ValueError("converters must match the number of values")

    arrays = np.broadcast_arrays(*[np.asarray(value, dtype=np.float64) for value in values])
    if all(array.ndim == 0 for array in arrays):
        args = [converter(float(array)) for converter, array in zip(converters, arrays)]
        return float(function(kernel_code, *args))

    output = np.empty(arrays[0].shape, dtype=np.float64)
    flat_arrays = [np.ascontiguousarray(array).ravel() for array in arrays]
    flat_output = output.ravel()
    for index, items in enumerate(zip(*flat_arrays)):
        args = [converter(float(item)) for converter, item in zip(converters, items)]
        flat_output[index] = function(kernel_code, *args)
    return output


def _broadcast_call(
    function: Callable[..., float],
    *values: ArrayLike,
    converters: Optional[tuple[Callable[[float], object], ...]] = None,
) -> Union[float, np.ndarray]:
    if converters is None:
        converters = tuple(float for _ in values)
    if len(converters) != len(values):
        raise ValueError("converters must match the number of values")

    arrays = np.broadcast_arrays(*[np.asarray(value, dtype=np.float64) for value in values])
    if all(array.ndim == 0 for array in arrays):
        args = [converter(float(array)) for converter, array in zip(converters, arrays)]
        return float(function(*args))

    output = np.empty(arrays[0].shape, dtype=np.float64)
    flat_arrays = [np.ascontiguousarray(array).ravel() for array in arrays]
    flat_output = output.ravel()
    for index, items in enumerate(zip(*flat_arrays)):
        args = [converter(float(item)) for converter, item in zip(converters, items)]
        flat_output[index] = function(*args)
    return output


def planck_occupation(x: ArrayLike) -> Union[float, np.ndarray]:
    """Return ``1 / (exp(x) - 1)`` using CSpack's stable implementation."""
    lib = _load_library()
    arr = np.asarray(x, dtype=np.float64)

    if arr.ndim == 0:
        return float(lib.nbb_func_C(float(arr)))

    flat = np.ascontiguousarray(arr.ravel())
    out = np.empty_like(flat)
    for index, value in enumerate(flat):
        out[index] = lib.nbb_func_C(float(value))
    return out.reshape(arr.shape)


def omega_min(omega0: ArrayLike, p0: ArrayLike) -> Union[float, np.ndarray]:
    """Return the lower outgoing-energy support for a fixed-momentum kernel."""
    lib = _load_library()
    _positive_array(omega0, "omega0")
    _positive_array(p0, "p0")
    return _broadcast_call(lib.cspack_omega_min, omega0, p0)


def omega_max(omega0: ArrayLike, p0: ArrayLike) -> Union[float, np.ndarray]:
    """Return the upper outgoing-energy support for a fixed-momentum kernel."""
    lib = _load_library()
    _positive_array(omega0, "omega0")
    _positive_array(p0, "p0")
    return _broadcast_call(lib.cspack_omega_max, omega0, p0)


def omega_crit(omega0: ArrayLike, p0: ArrayLike) -> Union[float, np.ndarray]:
    """Return the critical outgoing energy separating exact-kernel zones."""
    lib = _load_library()
    _positive_array(omega0, "omega0")
    _positive_array(p0, "p0")
    return _broadcast_call(lib.cspack_omega_crit, omega0, p0)


def kernel_bounds(omega0: ArrayLike, p0: ArrayLike) -> tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
    """Return ``(omega_min, omega_max)`` for a fixed-momentum kernel."""
    return omega_min(omega0, p0), omega_max(omega0, p0)


def kernel(
    omega0: ArrayLike,
    p0: ArrayLike,
    omega: ArrayLike,
    *,
    kind: Union[str, PhotonKernel] = PhotonKernel.EXACT,
) -> Union[float, np.ndarray]:
    """Evaluate a photon redistribution kernel.

    ``omega0``, ``p0``, and ``omega`` follow CSpack's dimensionless convention.
    Inputs can be scalars or NumPy-broadcastable arrays.
    """
    lib = _load_library()
    _positive_array(omega0, "omega0")
    _positive_array(p0, "p0")
    _positive_array(omega, "omega")
    return _broadcast_kernel_call(lib.cspack_kernel, _photon_kernel_code(kind), omega0, p0, omega)


def thermal_kernel(
    omega0: ArrayLike,
    omega: ArrayLike,
    theta: ArrayLike,
    *,
    kind: Union[str, ThermalPhotonKernel] = ThermalPhotonKernel.EXACT,
) -> Union[float, np.ndarray]:
    """Evaluate a thermally averaged photon redistribution kernel."""
    lib = _load_library()
    _positive_array(omega0, "omega0")
    _positive_array(omega, "omega")
    _positive_array(theta, "theta")
    return _broadcast_kernel_call(
        lib.cspack_thermal_kernel,
        _thermal_photon_kernel_code(kind),
        omega0,
        omega,
        theta,
    )


def neutrino_kernel(
    omega0: ArrayLike,
    p0: ArrayLike,
    omega: ArrayLike,
    *,
    kind: Union[str, NeutrinoKernel] = NeutrinoKernel.NUE_EP,
) -> Union[float, np.ndarray]:
    """Evaluate an exact neutrino redistribution kernel."""
    lib = _load_library()
    _positive_array(omega0, "omega0")
    _positive_array(p0, "p0")
    _positive_array(omega, "omega")
    return _broadcast_kernel_call(
        lib.cspack_neutrino_kernel,
        _neutrino_kernel_code(kind),
        omega0,
        p0,
        omega,
    )


def thermal_neutrino_kernel(
    omega0: ArrayLike,
    omega: ArrayLike,
    theta: ArrayLike,
    *,
    mue: float = 0.0,
    add_final_state: int = 0,
    kind: Union[str, NeutrinoKernel] = NeutrinoKernel.NUE_EP,
) -> Union[float, np.ndarray]:
    """Evaluate a Fermi-Dirac averaged neutrino redistribution kernel."""
    lib = _load_library()
    _positive_array(omega0, "omega0")
    _positive_array(omega, "omega")
    _positive_array(theta, "theta")
    mue_value = float(mue)
    if not np.isfinite(mue_value):
        raise ValueError("mue must be finite")
    add_final_state_value = int(add_final_state)
    if add_final_state_value not in (0, 1, 2, 3):
        raise ValueError("add_final_state must be 0, 1, 2, or 3")

    base = _broadcast_kernel_call(
        lib.cspack_thermal_neutrino_kernel,
        _neutrino_kernel_code(kind),
        omega0,
        omega,
        theta,
        np.asarray(mue_value),
        np.asarray(float(add_final_state_value)),
        converters=(float, float, float, float, int),
    )
    return base


def thermal_photon_fd_kernel(
    omega0: ArrayLike,
    omega: ArrayLike,
    theta: ArrayLike,
    *,
    mue: float = -1.0e2,
    add_final_state: int = 0,
    kind: Union[str, PhotonKernel] = PhotonKernel.EXACT,
) -> Union[float, np.ndarray]:
    """Evaluate a Fermi-Dirac averaged photon-electron redistribution kernel."""
    lib = _load_library()
    _positive_array(omega0, "omega0")
    _positive_array(omega, "omega")
    _positive_array(theta, "theta")
    mue_value = float(mue)
    if not np.isfinite(mue_value):
        raise ValueError("mue must be finite")
    add_final_state_value = int(add_final_state)
    if add_final_state_value not in (0, 1, 2, 3):
        raise ValueError("add_final_state must be 0, 1, 2, or 3")

    return _broadcast_kernel_call(
        lib.cspack_thermal_photon_fd_kernel,
        _photon_kernel_code(kind),
        omega0,
        omega,
        theta,
        np.asarray(mue_value),
        np.asarray(float(add_final_state_value)),
        converters=(float, float, float, float, int),
    )


def neutrino_kernel_norm(kind: Union[str, NeutrinoKernel] = NeutrinoKernel.NUE_EP) -> float:
    """Return the weak-coupling normalization relative to ``nue_e``."""
    lib = _load_library()
    return float(lib.cspack_neutrino_kernel_norm(ctypes.c_int(_neutrino_kernel_code(kind))))


def integral_weights(x: ArrayLike) -> np.ndarray:
    """Compute CSpack quadrature weights for a strictly increasing x grid."""
    lib = _load_library()
    grid = _as_grid(x)
    weights = np.empty_like(grid)
    lib.Integral_weights(grid, weights, ctypes.c_int(grid.size))
    return weights


def scattering_matrix_explicit(
    x: ArrayLike,
    theta: float,
    *,
    weights: Optional[ArrayLike] = None,
    kernel: Union[str, KernelType] = KernelType.EXACT,
    epsilon: float = 1.0e-8,
) -> np.ndarray:
    """Compute an explicit photon scattering matrix.

    The returned matrix has shape ``(len(x), len(x))`` and follows CSpack's
    ``Msc[i, j] = w_j P_ij theta`` convention.
    """
    lib = _load_library()
    grid = _as_grid(x)
    grid_weights = _as_weights(weights, grid)
    theta_value = _positive_float(theta, "theta")
    epsilon_value = _positive_float(epsilon, "epsilon")
    matrix = np.zeros((grid.size, grid.size), dtype=np.float64, order="C")
    rows = _matrix_row_pointers(matrix)

    lib.compute_scattering_matrix_explicit(
        grid,
        grid_weights,
        ctypes.c_int(grid.size),
        ctypes.c_double(theta_value),
        ctypes.c_int(_matrix_kernel_code(kernel)),
        ctypes.c_double(epsilon_value),
        rows,
    )
    return matrix


def scattering_matrix_kr(
    x: ArrayLike,
    theta: float,
    *,
    weights: Optional[ArrayLike] = None,
    n_kernel: int = 80,
    epsilon: float = 1.0e-8,
) -> np.ndarray:
    """Compute a photon scattering matrix using CSpack's kernel representation."""
    lib = _load_library()
    grid = _as_grid(x)
    grid_weights = _as_weights(weights, grid)
    theta_value = _positive_float(theta, "theta")
    epsilon_value = _positive_float(epsilon, "epsilon")
    n_kernel_value = int(n_kernel)
    if n_kernel_value < 5:
        raise ValueError("n_kernel must be at least 5")

    matrix = np.zeros((grid.size, grid.size), dtype=np.float64, order="C")
    rows = _matrix_row_pointers(matrix)

    lib.compute_scattering_matrix_KR(
        grid,
        grid_weights,
        ctypes.c_int(grid.size),
        ctypes.c_int(n_kernel_value),
        ctypes.c_double(theta_value),
        ctypes.c_double(epsilon_value),
        rows,
    )
    return matrix
