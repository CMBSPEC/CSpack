from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Optional, Union


PACKAGE_DIR = Path(__file__).resolve().parent
PYTHON_DIR = PACKAGE_DIR.parent
REPO_ROOT = PYTHON_DIR.parent

SOURCES = [
    "python/cspack/_bridge.cpp",
    "Tools/Simple_routines/routines.cpp",
    "Tools/Simple_routines/parser.cpp",
    "Tools/ODE_PDE_Solver/ODE_solver_LA.cpp",
    "Tools/ODE_PDE_Solver/ODE_solver.cpp",
    "Tools/Compton_Kernel/Compton_Kernel.cpp",
    "Tools/Integration/Patterson.cpp",
    "Tools/Integration/Integration_routines.GSL.cpp",
    "Tools/Integration/Chebyshev_Int.cpp",
    "Tools/Cosmology/saikawa_shirai_gstar.cpp",
    "CSpack/CSpack_functions.cpp",
    "CSpack/CSpack_kernels.cpp",
    "CSpack/CSpack_kernel_moments.cpp",
    "CSpack/CSpack_scattering_matrix.cpp",
    "CSpack/CSpack_energy_losses.cpp",
    "CSpack/CSpack_kernel_representation.cpp",
    "CSpack/CSpack_weights.cpp",
    "CSpack/CSpack_equilibrium_solutions.cpp",
]

INCLUDE_DIRS = [
    "CSpack",
    "Tools/Definitions",
    "Tools/Integration",
    "Tools/Compton_Kernel",
    "Tools/ODE_PDE_Solver",
    "Tools/Simple_routines",
    "Tools/Cosmology",
    ".",
]

DEFAULT_PREFIXES = [
    "/opt/homebrew",
    "/usr/local",
    "/opt/local",
]


def _library_filename() -> str:
    if sys.platform == "darwin":
        return "libcspack_py.dylib"
    if os.name == "nt":
        return "cspack_py.dll"
    return "libcspack_py.so"


def default_output_path() -> Path:
    return PACKAGE_DIR / "lib" / _library_filename()


def _default_cxx() -> str:
    env_cxx = os.environ.get("CXX")
    if env_cxx:
        return env_cxx

    candidates = []
    for prefix in DEFAULT_PREFIXES:
        candidates.extend(Path(prefix, "bin").glob("g++-*"))

    def version_key(path: Path) -> tuple[int, ...]:
        match = re.fullmatch(r"g\+\+-(\d+(?:\.\d+)*)", path.name)
        if match is None:
            return ()
        return tuple(int(part) for part in match.group(1).split("."))

    for candidate in sorted(candidates, key=version_key, reverse=True):
        if version_key(candidate) and candidate.exists():
            return str(candidate)

    return shutil.which("g++") or shutil.which("c++") or "g++"


def _existing_dirs(paths: Iterable[Path]) -> list[Path]:
    return [path for path in paths if path.exists()]


def _env_or_default_dir(env_name: str, suffix: str) -> list[Path]:
    value = os.environ.get(env_name)
    if value:
        return [Path(value).expanduser()]
    return _existing_dirs(Path(prefix) / suffix for prefix in DEFAULT_PREFIXES)


def _build_command(
    output: Path,
    cxx: str,
    gsl_include: Optional[Path],
    gsl_lib: Optional[Path],
    boost_include: Optional[Path],
    openmp: bool,
) -> list[str]:
    include_dirs = [REPO_ROOT / path for path in INCLUDE_DIRS]
    include_dirs.extend(_env_or_default_dir("GSL_INC_PATH", "include"))
    include_dirs.extend(_env_or_default_dir("BOOST_INC_PATH", "include"))
    if gsl_include is not None:
        include_dirs.insert(0, gsl_include)
    if boost_include is not None:
        include_dirs.insert(0, boost_include)

    lib_dirs = _env_or_default_dir("GSL_LIB_PATH", "lib")
    if gsl_lib is not None:
        lib_dirs.insert(0, gsl_lib)

    if sys.platform == "darwin":
        shared_flags = ["-dynamiclib", "-install_name", f"@rpath/{output.name}"]
    else:
        shared_flags = ["-shared"]

    command = [
        cxx,
        "-std=c++14",
        "-Wall",
        "-pedantic",
        "-O2",
        "-fPIC",
    ]

    if openmp:
        command.extend(["-D", "OPENMP_ACTIVATED", "-fopenmp", "-pthread"])

    command.extend(f"-I{path}" for path in include_dirs)
    command.extend(str(REPO_ROOT / source) for source in SOURCES)
    command.extend(shared_flags)
    command.extend(["-o", str(output)])
    command.extend(f"-L{path}" for path in lib_dirs)
    command.extend(f"-Wl,-rpath,{path}" for path in lib_dirs)
    command.extend(["-lgsl", "-lgslcblas", "-lm"])

    return command


def _run(command: list[str], verbose: bool) -> None:
    if verbose:
        print(" ".join(command))
        subprocess.run(command, cwd=REPO_ROOT, check=True)
        return

    result = subprocess.run(command, cwd=REPO_ROOT, text=True, capture_output=True)
    if result.returncode == 0:
        return

    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
    raise subprocess.CalledProcessError(
        result.returncode,
        command,
        output=result.stdout,
        stderr=result.stderr,
    )


def build_shared_library(
    output: Optional[Union[os.PathLike[str], str]] = None,
    cxx: Optional[str] = None,
    gsl_include: Optional[Union[os.PathLike[str], str]] = None,
    gsl_lib: Optional[Union[os.PathLike[str], str]] = None,
    boost_include: Optional[Union[os.PathLike[str], str]] = None,
    openmp: Optional[bool] = None,
    verbose: bool = True,
) -> Path:
    """Compile CSpack's C ABI into a shared library for the ctypes wrapper."""
    output_path = Path(output).expanduser() if output is not None else default_output_path()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cxx_value = cxx or _default_cxx()
    gsl_include_path = Path(gsl_include).expanduser() if gsl_include is not None else None
    gsl_lib_path = Path(gsl_lib).expanduser() if gsl_lib is not None else None
    boost_include_path = Path(boost_include).expanduser() if boost_include is not None else None

    attempts = [openmp] if openmp is not None else [True, False]
    last_error: Optional[subprocess.CalledProcessError] = None
    for use_openmp in attempts:
        command = _build_command(
            output_path,
            cxx_value,
            gsl_include_path,
            gsl_lib_path,
            boost_include_path,
            bool(use_openmp),
        )
        try:
            _run(command, verbose=verbose)
            return output_path
        except subprocess.CalledProcessError as exc:
            last_error = exc
            if openmp is not None or not use_openmp:
                break
            if verbose:
                print("OpenMP build failed; retrying without OpenMP.", file=sys.stderr)

    assert last_error is not None
    raise last_error


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Build the CSpack shared library for Python.")
    parser.add_argument("--output", type=Path, default=default_output_path())
    parser.add_argument("--cxx", default=None, help="C++ compiler. Defaults to CXX, g++-15, or g++.")
    parser.add_argument("--gsl-include", type=Path, default=None)
    parser.add_argument("--gsl-lib", type=Path, default=None)
    parser.add_argument("--boost-include", type=Path, default=None)
    parser.add_argument("--no-openmp", action="store_true", help="Compile without OpenMP.")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    output = build_shared_library(
        output=args.output,
        cxx=args.cxx,
        gsl_include=args.gsl_include,
        gsl_lib=args.gsl_lib,
        boost_include=args.boost_include,
        openmp=False if args.no_openmp else None,
        verbose=not args.quiet,
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
