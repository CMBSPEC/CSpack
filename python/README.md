# CSpack Python Wrapper

This folder contains a small Python wrapper around CSpack. It uses the stable C
ABI declared in `../CSpack/CSpack_c.h` where possible, plus a tiny local bridge
for direct kernel calls. It intentionally exposes only a useful subset:

- `kernel(omega0, p0, omega, kind=...)`
- `thermal_kernel(omega0, omega, theta, kind=...)`
- `thermal_photon_fd_kernel(omega0, omega, theta, kind=...)`
- `neutrino_kernel(omega0, p0, omega, kind=...)`
- `thermal_neutrino_kernel(omega0, omega, theta, kind=...)`
- `omega_min(omega0, p0)`, `omega_max(omega0, p0)`, `omega_crit(omega0, p0)`
- `kernel_bounds(omega0, p0)`
- `planck_occupation(x)`
- `integral_weights(x_grid)`
- `scattering_matrix_explicit(...)` for `exact`, `SS_K`, and `SS_C`
- `scattering_matrix_kr(...)` using the CSpack kernel representation

The wrapper uses `ctypes`, so no Python extension module is compiled. A local
shared library is built from the existing C++ sources.

## Use

From the main CSpack directory:

```sh
make py
make notebook
```

That builds the local shared library and opens
`python/notebooks/CSpack_examples.ipynb` in JupyterLab.
The notebook includes cached PNG figures, so the examples are visible
immediately when it opens.

## Requirements

On macOS/Homebrew, the project expects:

```sh
brew install gcc gsl boost jupyterlab python-matplotlib
```

The Makefile prefers `/opt/homebrew/bin/python3` and `/opt/homebrew/bin/jupyter`
when they exist. The build helper uses the newest `/opt/homebrew/bin/g++-*`
compiler it can find, then falls back to `g++`. If needed, you can override
paths:

```sh
CXX=/path/to/g++ python -m cspack.build
python -m cspack.build --gsl-include /opt/homebrew/include --gsl-lib /opt/homebrew/lib
```

To disable OpenMP:

```sh
python -m cspack.build --no-openmp
```

## Quick Example

```python
import numpy as np
import cspack

cspack.ensure_library()

omega = np.logspace(-4, -2, 100)
kernel_values = cspack.thermal_kernel(
    omega0=1.0e-3,
    omega=omega,
    theta=1.0e-3,
    kind="SS_C",
)

omega_min, omega_max = cspack.kernel_bounds(omega0=0.1, p0=0.5)
omega = np.linspace(omega_min, omega_max, 100)
photon_kernel = cspack.kernel(omega0=0.1, p0=0.5, omega=omega)

x = np.logspace(-3, 1, 64)
weights = cspack.integral_weights(x)
matrix = cspack.scattering_matrix_explicit(
    x,
    theta=0.1,
    weights=weights,
    kernel="exact",
    epsilon=1.0e-8,
)
```

## Notebook

Run:

```sh
make notebook
```

The notebook builds/loads the local shared library and plots photon and neutrino
kernels. It includes paper-style fixed-momentum and thermally averaged kernel
figures, photon full/recoil/Doppler comparisons, and separate neutrino
scattering-channel examples.

Inside JupyterLab, run the `Setup` cell once, then run any example or plotting
cell with `Shift-Enter`. The Markdown images are cached previews; the cells
below them recompute the same plots from the Python wrapper. The plot helpers
accept arguments such as `omega0`, `p0`, `p_values`, `omega1_values`, and
`theta`, so you can change examples directly in the notebook cells.

To regenerate the cached notebook figures without opening JupyterLab:

```sh
make -C python figures
```
