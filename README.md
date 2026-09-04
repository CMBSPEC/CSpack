# CSpack

CSpack is a C++ library for evaluating photon-electron Compton scattering
kernels, neutrino-electron and neutrino-positron scattering kernels, thermal
averages, moments, and scattering matrices in isotropic media.

The main documentation is in:

- [docs/cspack_kernel_library.pdf](docs/cspack_kernel_library.pdf)
- [docs/cspack_kernel_library.md](docs/cspack_kernel_library.md)
- [docs/cspack_kernel_library.tex](docs/cspack_kernel_library.tex)

Build the static library and example driver with:

```sh
make lib
make bin
```

The code examples in [CSpack_main.cpp](CSpack_main.cpp) show the basic kernel,
thermal-kernel, moment, opacity, and scattering-matrix calls.
