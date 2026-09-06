from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Iterable


NOTEBOOK_DIR = Path(__file__).resolve().parent
PYTHON_DIR = NOTEBOOK_DIR.parent
REPO_ROOT = PYTHON_DIR.parent
FIGURE_DIR = NOTEBOOK_DIR / "figures"
OUTPUT_DIR = REPO_ROOT / "outputs"

os.environ.setdefault("MPLCONFIGDIR", str(PYTHON_DIR / ".matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

if "ipykernel" not in sys.modules:
    import matplotlib

    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

sys.path.insert(0, str(PYTHON_DIR))
import cspack  # noqa: E402


P_VALUES = [1.0e-2, 5.0e-2, 1.0e-1, 2.0e-1, 5.0e-1]
OMEGA1_VALUES = [1.0e-3, 1.0e-2, 1.0e-1, 1.0, 10.0]
OMEGA1_LABELS = [
    r"$\omega_1=10^{-3}$",
    r"$\omega_1=10^{-2}$",
    r"$\omega_1=10^{-1}$",
    r"$\omega_1=1$",
    r"$\omega_1=10$",
]
THERMAL_COLORS = ["black", "#2369bd", "#7b3294", "#c9342c", "#2369bd"]
THERMAL_STYLES = ["-", "--", "-", "--", "-"]


def _table(name: str) -> np.ndarray:
    path = OUTPUT_DIR / name
    if not path.exists():
        raise FileNotFoundError(path)
    return np.loadtxt(path, comments="#")


def _positive(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return np.where(np.isfinite(arr) & (arr > 0.0), arr, np.nan)


def _as_float_list(values: Iterable[float] | None, default: list[float]) -> list[float]:
    if values is None:
        return list(default)
    return [float(value) for value in values]


def _same_values(values: list[float], default: list[float]) -> bool:
    return len(values) == len(default) and bool(np.allclose(values, default))


def _kernel_bounds(omega0: float, p0: float) -> tuple[float, float]:
    omega_min, omega_max = cspack.kernel_bounds(omega0, p0)
    omega_min = float(omega_min)
    omega_max = float(omega_max)
    if not np.isfinite(omega_min) or not np.isfinite(omega_max) or omega_max <= omega_min:
        raise ValueError(f"Invalid kernel support for omega0={omega0:g}, p0={p0:g}")
    return omega_min, omega_max


def _checked_omega_grid(omega: Iterable[float], omega_min: float, omega_max: float) -> np.ndarray:
    grid = np.asarray(list(omega), dtype=float)
    keep = np.isfinite(grid) & (grid >= omega_min) & (grid <= omega_max)
    grid = np.unique(np.sort(grid[keep]))
    if grid.size < 2:
        raise ValueError(
            f"At least two omega values must lie inside [{omega_min:g}, {omega_max:g}]"
        )
    return grid


def _positive_omega_grid(omega: Iterable[float]) -> np.ndarray:
    grid = np.asarray(list(omega), dtype=float)
    grid = np.unique(np.sort(grid[np.isfinite(grid) & (grid > 0.0)]))
    if grid.size < 2:
        raise ValueError("At least two positive finite omega values are required")
    return grid


def _kernel_grid(
    omega0: float,
    p0: float,
    n_omega: int,
    omega: Iterable[float] | None = None,
) -> np.ndarray:
    omega_min, omega_max = _kernel_bounds(omega0, p0)
    if omega is not None:
        return _checked_omega_grid(omega, omega_min, omega_max)
    return np.linspace(omega_min, omega_max, int(n_omega))


def _combined_kernel_grid(
    omega0: float,
    p_values: list[float],
    n_omega: int,
    omega: Iterable[float] | None = None,
) -> np.ndarray:
    bounds = np.array([_kernel_bounds(omega0, p0) for p0 in p_values])
    omega_min = float(bounds[:, 0].min())
    omega_max = float(bounds[:, 1].max())
    if omega is not None:
        return _checked_omega_grid(omega, omega_min, omega_max)
    return np.linspace(omega_min, omega_max, int(n_omega))


def _auto_log_ylim(ax: plt.Axes, *curves: np.ndarray, floor_fraction: float = 1.0e-4) -> None:
    values = np.concatenate([np.ravel(curve) for curve in curves])
    values = values[np.isfinite(values) & (values > 0.0)]
    if values.size == 0:
        return

    upper = float(values.max())
    lower = max(float(values.min()), upper * floor_fraction)
    if lower >= upper:
        lower = upper / 10.0
    ax.set_ylim(lower / 1.4, upper * 1.4)


def _format_axes(ax: plt.Axes) -> None:
    ax.tick_params(direction="in", which="both", top=True, right=True)
    ax.grid(True, which="major", color="0.82", linewidth=0.7)
    ax.grid(True, which="minor", color="0.9", linewidth=0.45)


def _save(fig: plt.Figure, filename: str | None, show: bool) -> plt.Figure:
    fig.tight_layout()
    if filename is not None:
        FIGURE_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIGURE_DIR / filename, dpi=180, bbox_inches="tight")
    if not show:
        plt.close(fig)
    return fig


def _from_tables_or_fixed_kernel_tables() -> tuple[np.ndarray, np.ndarray]:
    try:
        return _table("kernel_ph.dat"), _table("kernel_nu.dat")
    except FileNotFoundError:
        omega0 = 0.1
        omega = np.linspace(0.03399463336, 0.1977954475, 1000)
        ph = np.column_stack(
            [omega]
            + [cspack.kernel(omega0, p0, omega, kind="exact") for p0 in P_VALUES]
        )
        nu = np.column_stack(
            [omega]
            + [cspack.neutrino_kernel(omega0, p0, omega, kind="nue_e") for p0 in P_VALUES]
        )
        return ph, nu


def _from_tables_or_thermal(
    base_name: str,
    final_name: str,
    *,
    photon: bool,
    omega1_values: list[float],
    theta: float,
    mue: float,
    n_omega: int,
    omega: Iterable[float] | None,
    omega_range: tuple[float, float] | None,
    use_tables: bool,
) -> tuple[np.ndarray, np.ndarray]:
    use_reference_tables = (
        use_tables
        and np.isclose(theta, 0.1)
        and np.isclose(mue, -1.0e2)
        and _same_values(omega1_values, OMEGA1_VALUES)
        and omega is None
        and omega_range is None
    )
    if use_reference_tables:
        try:
            return _table(base_name), _table(final_name)
        except FileNotFoundError:
            pass

    if omega is not None:
        omega = _positive_omega_grid(omega)
    elif omega_range is None:
        omega_low = min(omega1_values) * 1.0e-3
        omega_high = max(omega1_values) * 1.0e3
        omega = np.logspace(np.log10(float(omega_low)), np.log10(float(omega_high)), int(n_omega))
    else:
        omega_low, omega_high = omega_range
        omega = np.logspace(np.log10(float(omega_low)), np.log10(float(omega_high)), int(n_omega))
    if photon:
        plain = np.column_stack(
            [omega]
            + [
                cspack.thermal_photon_fd_kernel(
                    omega0, omega, theta, mue=mue, add_final_state=0, kind="exact"
                )
                for omega0 in omega1_values
            ]
        )
        final = np.column_stack(
            [omega]
            + [
                cspack.thermal_photon_fd_kernel(
                    omega0, omega, theta, mue=mue, add_final_state=3, kind="exact"
                )
                for omega0 in omega1_values
            ]
        )
    else:
        plain = np.column_stack(
            [omega]
            + [
                cspack.thermal_neutrino_kernel(
                    omega0, omega, theta, mue=mue, add_final_state=0, kind="nue_e"
                )
                for omega0 in omega1_values
            ]
        )
        final = np.column_stack(
            [omega]
            + [
                cspack.thermal_neutrino_kernel(
                    omega0, omega, theta, mue=mue, add_final_state=2, kind="nue_e"
                )
                for omega0 in omega1_values
            ]
        )
    return plain, final


def plot_fixed_photon_neutrino(
    filename: str | None = "fixed_photon_neutrino_kernels.png",
    *,
    show: bool = False,
    omega0: float = 0.1,
    p_values: Iterable[float] | None = None,
    n_omega: int = 1000,
    omega: Iterable[float] | None = None,
    use_tables: bool = True,
) -> plt.Figure:
    p_values = _as_float_list(p_values, P_VALUES)
    use_reference_tables = (
        use_tables
        and np.isclose(omega0, 0.1)
        and _same_values(p_values, P_VALUES)
        and omega is None
    )
    if use_reference_tables:
        photon, neutrino = _from_tables_or_fixed_kernel_tables()
    else:
        omega = _combined_kernel_grid(omega0, p_values, n_omega, omega=omega)
        photon = np.column_stack(
            [omega]
            + [cspack.kernel(omega0, p0, omega, kind="exact") for p0 in p_values]
        )
        neutrino = np.column_stack(
            [omega]
            + [cspack.neutrino_kernel(omega0, p0, omega, kind="nue_e") for p0 in p_values]
        )
    omega = photon[:, 0]

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    plotted_curves = []
    for index, p0 in enumerate(p_values):
        neutrino_values = _positive(neutrino[:, index + 1])
        photon_values = _positive(photon[:, index + 1])
        plotted_curves.extend([neutrino_values, photon_values])
        ax.semilogy(
            omega,
            neutrino_values,
            color="black",
            linewidth=1.35,
            alpha=0.52 + 0.4 * index / max(len(p_values) - 1, 1),
            label="neutrino-electron" if index == 0 else None,
        )
        ax.semilogy(
            omega,
            photon_values,
            color="#2369bd",
            linewidth=1.35,
            alpha=0.52 + 0.4 * index / max(len(p_values) - 1, 1),
            label="photon-electron" if index == 0 else None,
        )
    ax.set_xlim(0.0 if use_reference_tables else float(omega.min()), float(omega.max()))
    if use_reference_tables:
        ax.set_ylim(1.0e-1, 1.0e2)
    else:
        _auto_log_ylim(ax, *plotted_curves)
    ax.set_xlabel(r"$\omega_3$")
    ax.set_ylabel("Kernel")
    ax.set_title(fr"Fixed scatterer momentum, $\omega_1={omega0:g}$")
    ax.legend(frameon=False, loc="upper right")
    _format_axes(ax)
    return _save(fig, filename, show)


def _plot_thermal_from_tables(
    plain: np.ndarray,
    final: np.ndarray,
    *,
    title: str,
    final_label: str,
    ylim: tuple[float, float],
    omega1_values: list[float],
    theta: float,
    mue: float,
    filename: str | None,
    show: bool,
) -> plt.Figure:
    omega = plain[:, 0]
    fig, ax = plt.subplots(figsize=(7.4, 4.8))

    omega1_labels = [fr"$\omega_1={value:g}$" for value in omega1_values]
    if _same_values(omega1_values, OMEGA1_VALUES):
        omega1_labels = OMEGA1_LABELS

    for index, label in enumerate(omega1_labels):
        color = THERMAL_COLORS[index % len(THERMAL_COLORS)]
        style = THERMAL_STYLES[index % len(THERMAL_STYLES)]
        ax.loglog(
            omega,
            _positive(plain[:, index + 1]),
            linestyle=style,
            color=color,
            linewidth=1.55,
            label=label,
        )
        ax.loglog(
            omega,
            _positive(final[:, index + 1]),
            linestyle=style,
            color=color,
            linewidth=0.9,
            alpha=0.42,
        )

    ax.set_xlim(5.0e-5, 20.0)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$\omega_3$")
    ax.set_ylabel("Kernel")
    ax.set_title(title)
    ax.text(
        0.04,
        0.94,
        fr"$\theta_e={theta:g},\ \mu_e={mue:g}$",
        transform=ax.transAxes,
        ha="left",
        va="top",
    )
    first = ax.legend(frameon=False, loc="upper right", title="incident energy")
    ax.add_artist(first)
    ax.legend(
        handles=[
            Line2D([0], [0], color="0.2", lw=1.55, label="no final-state factor"),
            Line2D([0], [0], color="0.2", lw=0.9, alpha=0.42, label=final_label),
        ],
        frameon=False,
        loc="lower left",
    )
    _format_axes(ax)
    return _save(fig, filename, show)


def plot_thermal_photon_fd(
    filename: str | None = "thermal_photon_fd_kernels.png",
    *,
    show: bool = False,
    omega1_values: Iterable[float] | None = None,
    theta: float = 0.1,
    mue: float = -1.0e2,
    n_omega: int = 120,
    omega: Iterable[float] | None = None,
    omega_range: tuple[float, float] | None = None,
    use_tables: bool = True,
) -> plt.Figure:
    omega1_values = _as_float_list(omega1_values, OMEGA1_VALUES)
    plain, stimulated = _from_tables_or_thermal(
        "kernel_ph_FD.dat",
        "kernel_ph_FD_BB.dat",
        photon=True,
        omega1_values=omega1_values,
        theta=theta,
        mue=mue,
        n_omega=n_omega,
        omega=omega,
        omega_range=omega_range,
        use_tables=use_tables,
    )
    return _plot_thermal_from_tables(
        plain,
        stimulated,
        title="Thermally averaged photon-electron kernels",
        final_label="electron blocking + photon stimulation",
        ylim=(1.0e-5, 2.0e5),
        omega1_values=omega1_values,
        theta=theta,
        mue=mue,
        filename=filename,
        show=show,
    )


def plot_thermal_neutrino_fd(
    filename: str | None = "thermal_neutrino_fd_kernels.png",
    *,
    show: bool = False,
    omega1_values: Iterable[float] | None = None,
    theta: float = 0.1,
    mue: float = -1.0e2,
    n_omega: int = 120,
    omega: Iterable[float] | None = None,
    omega_range: tuple[float, float] | None = None,
    use_tables: bool = True,
) -> plt.Figure:
    omega1_values = _as_float_list(omega1_values, OMEGA1_VALUES)
    plain, blocked = _from_tables_or_thermal(
        "kernel_nu_FD.dat",
        "kernel_nu_FD_FB.dat",
        photon=False,
        omega1_values=omega1_values,
        theta=theta,
        mue=mue,
        n_omega=n_omega,
        omega=omega,
        omega_range=omega_range,
        use_tables=use_tables,
    )
    return _plot_thermal_from_tables(
        plain,
        blocked,
        title="Thermally averaged electron-neutrino kernels",
        final_label="electron + neutrino blocking",
        ylim=(1.0e-5, 1.0e4),
        omega1_values=omega1_values,
        theta=theta,
        mue=mue,
        filename=filename,
        show=show,
    )


def _kernel_curve(
    omega0: float,
    p0: float,
    omega: np.ndarray,
    kind: str,
    support: np.ndarray | None = None,
) -> np.ndarray:
    values = np.asarray(cspack.kernel(omega0, p0, omega, kind=kind), dtype=float)
    if support is not None:
        values = np.where(support, values, np.nan)
    return _positive(values)


def plot_photon_approximations(
    filename: str | None = "photon_exact_recoil_doppler.png",
    *,
    show: bool = False,
    recoil_omega0: float = 0.1,
    recoil_p0: float = 1.0e-2,
    doppler_omega0: float = 1.0e-2,
    doppler_p0: float = 2.0e-1,
    n_omega: int = 900,
    recoil_omega: Iterable[float] | None = None,
    doppler_omega: Iterable[float] | None = None,
) -> plt.Figure:
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.1))

    omega = _kernel_grid(recoil_omega0, recoil_p0, n_omega, omega=recoil_omega)
    exact = _kernel_curve(recoil_omega0, recoil_p0, omega, "exact")
    support = np.isfinite(exact)
    recoil = _kernel_curve(recoil_omega0, recoil_p0, omega, "recoil", support)
    axes[0].semilogy(omega, exact, color="black", linewidth=1.6, label="full")
    axes[0].semilogy(
        omega,
        recoil,
        color="#c9342c",
        linestyle="--",
        linewidth=1.4,
        label="recoil",
    )
    axes[0].set_title(fr"Recoil limit: $\omega_1={recoil_omega0:g},\ p={recoil_p0:g}$")
    axes[0].set_xlabel(r"$\omega_3$")
    axes[0].set_ylabel("Kernel")
    axes[0].set_xlim(float(omega.min()), float(omega.max()))
    _auto_log_ylim(axes[0], exact, recoil)
    axes[0].legend(frameon=False)
    _format_axes(axes[0])

    omega = _kernel_grid(doppler_omega0, doppler_p0, n_omega, omega=doppler_omega)
    exact = _kernel_curve(doppler_omega0, doppler_p0, omega, "exact")
    support = np.isfinite(exact)
    doppler = _kernel_curve(doppler_omega0, doppler_p0, omega, "doppler", support)
    axes[1].semilogy(omega, exact, color="black", linewidth=1.6, label="full")
    axes[1].semilogy(
        omega,
        doppler,
        color="#2369bd",
        linestyle="--",
        linewidth=1.4,
        label="Doppler",
    )
    axes[1].set_title(fr"Doppler limit: $\omega_1={doppler_omega0:g},\ p={doppler_p0:g}$")
    axes[1].set_xlabel(r"$\omega_3$")
    axes[1].set_xlim(float(omega.min()), float(omega.max()))
    _auto_log_ylim(axes[1], exact, doppler)
    axes[1].legend(frameon=False)
    _format_axes(axes[1])

    return _save(fig, filename, show)


def plot_neutrino_channels(
    filename: str | None = "neutrino_scattering_channels.png",
    *,
    show: bool = False,
    omega0: float = 0.5,
    p0: float = 1.0,
    channels: Iterable[str] | None = None,
    n_omega: int = 1000,
    omega: Iterable[float] | None = None,
) -> plt.Figure:
    omega = _kernel_grid(omega0, p0, n_omega, omega=omega)
    default_channels = [
        ("nue_e", r"$\nu_e e^-$", "black", "-"),
        ("nue_p", r"$\nu_e e^+$", "#2369bd", "-"),
        ("nue_ep", r"$\nu_e(e^-+e^+)$", "#c9342c", "-"),
        ("numu_e", r"$\nu_{\mu,\tau} e^-$", "#7b3294", "--"),
        ("numu_p", r"$\nu_{\mu,\tau} e^+$", "#238b45", "--"),
    ]
    channel_styles = {kind: (label, color, style) for kind, label, color, style in default_channels}
    if channels is None:
        channels = [kind for kind, _, _, _ in default_channels]

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    plotted_curves = []
    for index, kind in enumerate(channels):
        label, color, style = channel_styles.get(
            kind,
            (kind, THERMAL_COLORS[index % len(THERMAL_COLORS)], THERMAL_STYLES[index % len(THERMAL_STYLES)]),
        )
        values = cspack.neutrino_kernel(omega0, p0, omega, kind=kind)
        values = _positive(values)
        plotted_curves.append(values)
        ax.plot(omega, values, color=color, linestyle=style, linewidth=1.55, label=label)

    ax.set_xlim(float(omega.min()), float(omega.max()))
    finite = np.concatenate([np.ravel(curve) for curve in plotted_curves])
    finite = finite[np.isfinite(finite)]
    if finite.size:
        ax.set_ylim(0.0, float(finite.max()) * 1.08)
    ax.set_xlabel(r"$\omega_3$")
    ax.set_ylabel("Kernel")
    ax.set_title(fr"Neutrino scattering channels, $\omega_1={omega0:g},\ p={p0:g}$")
    ax.legend(frameon=False, loc="upper left")
    _format_axes(ax)
    return _save(fig, filename, show)


def plot_scattering_matrix(
    filename: str | None = "scattering_matrix_elements.png",
    *,
    show: bool = False,
    theta: float = 0.1,
    x: Iterable[float] | None = None,
    n_x: int = 48,
    kernel: str = "exact",
    rows: Iterable[int] | None = None,
) -> plt.Figure:
    if x is None:
        x = np.logspace(-3.0, 1.0, int(n_x))
    else:
        x = np.asarray(list(x), dtype=float)
    if rows is None:
        rows = [12, 20, 28, 36]

    weights = cspack.integral_weights(x)
    matrix = cspack.scattering_matrix_explicit(
        x,
        theta=theta,
        weights=weights,
        kernel=kernel,
        epsilon=1.0e-8,
    )
    positive = _positive(matrix)
    matrix_floor = np.nanmax(positive) * 1.0e-12
    log_matrix = np.log10(np.where(positive >= matrix_floor, positive, np.nan))

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.1))
    mesh = axes[0].pcolormesh(
        x,
        x,
        log_matrix,
        shading="auto",
        cmap="viridis",
        vmin=np.nanmax(log_matrix) - 12.0,
        vmax=np.nanmax(log_matrix),
    )
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"outgoing grid $x_j$")
    axes[0].set_ylabel(r"incoming grid $x_i$")
    axes[0].set_title(fr"$\log_{{10}} M_{{ij}}$, {kernel} kernel, $\theta_e={theta:g}$")
    fig.colorbar(mesh, ax=axes[0], label=r"$\log_{10} M_{ij}$")
    _format_axes(axes[0])

    for row in rows:
        row = int(row)
        row_values = positive[row]
        row_floor = np.nanmax(row_values) * 1.0e-10
        row_values = np.where(row_values >= row_floor, row_values, np.nan)
        axes[1].loglog(x, row_values, linewidth=1.35, label=fr"$x_i={x[row]:.2g}$")
    axes[1].set_xlabel(r"outgoing grid $x_j$")
    axes[1].set_ylabel(r"$M_{ij}$")
    axes[1].set_title(fr"Selected matrix rows, $\theta_e={theta:g}$")
    _auto_log_ylim(axes[1], positive)
    axes[1].legend(frameon=False)
    _format_axes(axes[1])
    return _save(fig, filename, show)


def render_all(*, write_notebook: bool = False) -> list[Path]:
    cspack.ensure_library()
    plot_fixed_photon_neutrino()
    plot_thermal_photon_fd()
    plot_thermal_neutrino_fd()
    plot_photon_approximations()
    plot_neutrino_channels()
    plot_scattering_matrix()
    paths = sorted(FIGURE_DIR.glob("*.png"))
    if write_notebook:
        write_notebook_file()
    return paths


def _source(text: str) -> list[str]:
    return text.strip("\n").splitlines(keepends=True)


def _markdown(text: str) -> dict[str, object]:
    return {"cell_type": "markdown", "metadata": {}, "source": _source(text)}


def _code(text: str) -> dict[str, object]:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": _source(text),
    }


def write_notebook_file() -> Path:
    path = NOTEBOOK_DIR / "CSpack_examples.ipynb"
    notebook = {
        "cells": [
            _markdown(
                """
# CSpack Python Kernel Examples

The figures below are cached PNGs, so they are visible as soon as the notebook opens. Run the code cells to rebuild the local shared library and regenerate the plots from the Python wrapper.
"""
            ),
            _markdown(
                """
## References

- Chluba, Cyr, Uwabo-Niibo, and Yamaguchi, *Neutrino-electron scattering kernels in isotropic media*, arXiv:2607.25436.
- Sarkar, Chluba, and Lee, *Dissecting the Compton scattering kernel I: Isotropic media*, arXiv:1905.00868.
"""
            ),
            _markdown(
                """
## How to Run This Notebook

The images are already visible because they are saved as PNG files in `python/notebooks/figures`.

To run the examples yourself:

1. Run the `Setup` code cell once.
2. Run any example cell with `Shift-Enter`.
3. Plot cells use `show=True`, which displays the recomputed figure in the notebook.
4. Change plot arguments directly in the notebook cell, for example `plot_photon_approximations(recoil_p0=0.02, doppler_omega0=0.01, show=True)`.
5. For fixed-kernel plots, explicit `omega` grids are checked against `cspack.kernel_bounds` before plotting.
6. To refresh every cached PNG from the terminal, run `make -C python figures` from the main CSpack directory.
"""
            ),
            _markdown(
                """
## Setup

Run this cell first. It finds the repo root, imports the wrapper, imports the plotting helpers, and builds the local shared library if it is missing.
"""
            ),
            _code(
                """
from pathlib import Path
import os
import sys
import numpy as np

cwd = Path.cwd().resolve()
for base in (cwd, *cwd.parents):
    if (base / "python" / "notebooks" / "render_figures.py").exists():
        repo_root = base
        break
    if (base / "notebooks" / "render_figures.py").exists():
        repo_root = base.parent
        break
else:
    raise RuntimeError("Could not locate the CSpack repository root.")

python_dir = repo_root / "python"
notebook_dir = python_dir / "notebooks"
os.environ.setdefault("MPLCONFIGDIR", str(python_dir / ".matplotlib"))
sys.path.insert(0, str(python_dir))
sys.path.insert(0, str(notebook_dir))

import cspack
from render_figures import (
    plot_fixed_photon_neutrino,
    plot_thermal_photon_fd,
    plot_thermal_neutrino_fd,
    plot_photon_approximations,
    plot_neutrino_channels,
    plot_scattering_matrix,
    render_all,
)

cspack.ensure_library()
"""
            ),
            _markdown(
                """
## Calling the Wrapper

The wrapper accepts scalars or NumPy arrays for the kernel variables. These examples use small grids so they run quickly while showing the call signatures.
"""
            ),
            _markdown(
                """
### Photon Kernel Calls

Use `cspack.kernel` for fixed-momentum photon-electron kernels. The `kind` argument selects `exact`, `recoil`, `doppler`, or `ur`. Use `cspack.kernel_bounds` to get the physical outgoing-energy range before making a grid.
"""
            ),
            _code(
                """
omega_min, omega_max = cspack.kernel_bounds(omega0=0.1, p0=0.5)
omega = np.linspace(omega_min, omega_max, 8)

doppler_min, doppler_max = cspack.kernel_bounds(omega0=1.0e-2, p0=0.2)
doppler_omega = np.linspace(doppler_min, doppler_max, 10)[1:-1]

photon_exact = cspack.kernel(omega0=0.1, p0=0.5, omega=omega, kind="exact")
photon_doppler = cspack.kernel(
    omega0=1.0e-2,
    p0=0.2,
    omega=doppler_omega,
    kind="doppler",
)

print("omega grid:", omega)
print("exact photon kernel:", photon_exact)
print("Doppler omega grid:", doppler_omega)
print("Doppler approximation:", photon_doppler)
"""
            ),
            _markdown(
                """
### Thermal Photon Calls

Use `thermal_photon_fd_kernel` for Fermi-Dirac averaged photon kernels. `add_final_state=3` applies electron blocking and photon stimulation.
"""
            ),
            _code(
                """
thermal_omega = np.logspace(-4, -1, 12)

photon_thermal = cspack.thermal_photon_fd_kernel(
    omega0=1.0e-3,
    omega=thermal_omega,
    theta=0.1,
    mue=-100.0,
    add_final_state=0,
    kind="exact",
)
photon_stimulated = cspack.thermal_photon_fd_kernel(
    omega0=1.0e-3,
    omega=thermal_omega,
    theta=0.1,
    mue=-100.0,
    add_final_state=3,
    kind="exact",
)

print("thermal photon:", photon_thermal)
print("with photon stimulation:", photon_stimulated)
"""
            ),
            _markdown(
                """
### Neutrino Channel Calls

Use `neutrino_kernel` for fixed-momentum neutrino kernels. The channel name chooses the flavor and electron/positron target.
"""
            ),
            _code(
                """
omega = np.linspace(0.1, 0.9, 6)
channels = ["nue_e", "nue_p", "nue_ep", "numu_e", "numu_p"]

neutrino_values = {
    channel: cspack.neutrino_kernel(omega0=0.5, p0=1.0, omega=omega, kind=channel)
    for channel in channels
}

for channel, values in neutrino_values.items():
    print(channel, values)
"""
            ),
            _markdown(
                """
### Thermal Neutrino Calls

Use `thermal_neutrino_kernel` for Fermi-Dirac averaged neutrino kernels. `add_final_state=2` applies electron and neutrino blocking.
"""
            ),
            _code(
                """
nu_plain = cspack.thermal_neutrino_kernel(
    omega0=1.0e-3,
    omega=thermal_omega,
    theta=0.1,
    mue=-100.0,
    add_final_state=0,
    kind="nue_e",
)
nu_blocked = cspack.thermal_neutrino_kernel(
    omega0=1.0e-3,
    omega=thermal_omega,
    theta=0.1,
    mue=-100.0,
    add_final_state=2,
    kind="nue_e",
)

print("thermal neutrino:", nu_plain)
print("with final-state blocking:", nu_blocked)
"""
            ),
            _markdown(
                """
## Fixed-Momentum Photon and Neutrino Kernels

Paper-style comparison for `omega_1 = 0.1` and `p = 0.01, 0.05, 0.1, 0.2, 0.5`.

![Fixed-momentum photon and neutrino kernels](figures/fixed_photon_neutrino_kernels.png)
"""
            ),
            _code(
                """
# Run this cell with Shift-Enter to recompute and display this plot.
plot_fixed_photon_neutrino(
    omega0=0.1,
    p_values=[0.01, 0.05, 0.1, 0.2, 0.5],
    show=True,
)
"""
            ),
            _markdown(
                """
## Photon Kernels

Thermally averaged photon-electron redistribution for `theta_e = 0.1`, with and without final-state electron blocking plus photon stimulation.

![Thermally averaged photon kernels](figures/thermal_photon_fd_kernels.png)
"""
            ),
            _code(
                """
# Run this cell with Shift-Enter to recompute and display this plot.
plot_thermal_photon_fd(
    omega1_values=[1e-3, 1e-2, 1e-1, 1.0, 10.0],
    theta=0.1,
    show=True,
)
"""
            ),
            _markdown(
                """
### Full, Recoil, and Doppler Photon Limits

The approximation curves are shown only over the physical support of the exact kernel.
Change `recoil_omega0`, `recoil_p0`, `doppler_omega0`, or `doppler_p0` here; the plotted `omega` ranges are recomputed automatically. You can also pass `recoil_omega=` or `doppler_omega=` to provide a custom outgoing-energy grid.

![Photon full/recoil/Doppler comparison](figures/photon_exact_recoil_doppler.png)
"""
            ),
            _code(
                """
# Run this cell with Shift-Enter to recompute and display this plot.
plot_photon_approximations(
    recoil_omega0=0.1,
    recoil_p0=1.0e-2,
    doppler_omega0=1.0e-2,
    doppler_p0=0.2,
    show=True,
)
"""
            ),
            _markdown(
                """
## Neutrino Kernels

Thermally averaged electron-neutrino redistribution for `theta_e = 0.1`, with and without final-state electron and neutrino blocking.

![Thermally averaged neutrino kernels](figures/thermal_neutrino_fd_kernels.png)
"""
            ),
            _code(
                """
# Run this cell with Shift-Enter to recompute and display this plot.
plot_thermal_neutrino_fd(
    omega1_values=[1e-3, 1e-2, 1e-1, 1.0, 10.0],
    theta=0.1,
    show=True,
)
"""
            ),
            _markdown(
                """
### Neutrino Scattering Channels

Fixed-momentum channel comparison for `omega_1 = 0.5` and `p = 1`.

![Neutrino scattering channels](figures/neutrino_scattering_channels.png)
"""
            ),
            _code(
                """
# Run this cell with Shift-Enter to recompute and display this plot.
plot_neutrino_channels(
    omega0=0.5,
    p0=1.0,
    channels=["nue_e", "nue_p", "nue_ep", "numu_e", "numu_p"],
    show=True,
)
"""
            ),
            _markdown(
                """
## Scattering Matrix

This uses CSpack's quadrature weights and the explicit photon scattering-matrix wrapper at `theta_e = 0.1`. The left panel shows the full matrix; the right panel pulls out a few rows as ordinary matrix elements.

![Photon scattering matrix elements](figures/scattering_matrix_elements.png)
"""
            ),
            _code(
                """
x = np.logspace(-3, 1, 48)
weights = cspack.integral_weights(x)
matrix = cspack.scattering_matrix_explicit(
    x,
    theta=0.1,
    weights=weights,
    kernel="exact",
    epsilon=1.0e-8,
)

print("shape:", matrix.shape)
print("max element:", matrix.max())
"""
            ),
            _code(
                """
# Run this cell with Shift-Enter to recompute and display this plot.
plot_scattering_matrix(theta=0.1, show=True)
"""
            ),
            _markdown(
                """
## Regenerate All Cached Figures

This cell rewrites every PNG used by the Markdown image cells above.
"""
            ),
            _code(
                """
# Run this cell with Shift-Enter to rewrite every cached PNG in figures/.
render_all(write_notebook=False)
"""
            ),
        ],
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {
                "codemirror_mode": {"name": "ipython", "version": 3},
                "file_extension": ".py",
                "mimetype": "text/x-python",
                "name": "python",
                "nbconvert_exporter": "python",
                "pygments_lexer": "ipython3",
                "version": "3.14",
            },
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    path.write_text(json.dumps(notebook, indent=1) + "\n")
    return path


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Render cached figures for the CSpack notebook.")
    parser.add_argument("--write-notebook", action="store_true", help="Rewrite CSpack_examples.ipynb.")
    args = parser.parse_args(list(argv) if argv is not None else None)

    for path in render_all(write_notebook=args.write_notebook):
        print(path.relative_to(PYTHON_DIR))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
