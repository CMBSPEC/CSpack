#include <iostream>
#include <limits>

#include "CSpack.h"

namespace {

double nan_value()
{
    return std::numeric_limits<double>::quiet_NaN();
}

kernel_ptr neutrino_kernel_ptr(int kernel_type)
{
    switch(kernel_type)
    {
        case 0: return CSpack_kernels_nu::kernel_exact_nue_e;
        case 1: return CSpack_kernels_nu::kernel_exact_nue_p;
        case 2: return CSpack_kernels_nu::kernel_exact_nue_ep;
        case 3: return CSpack_kernels_nu::kernel_exact_numu_e;
        case 4: return CSpack_kernels_nu::kernel_exact_numu_p;
        case 5: return CSpack_kernels_nu::kernel_exact_nutau_e;
        case 6: return CSpack_kernels_nu::kernel_exact_nutau_p;
        default: return nullptr;
    }
}

kernel_ptr photon_kernel_ptr(int kernel_type)
{
    switch(kernel_type)
    {
        case 0: return CSpack_kernels::kernel_exact;
        case 1: return CSpack_kernels::kernel_recoil;
        case 2: return CSpack_kernels::kernel_doppler;
        case 3: return CSpack_kernels::kernel_ur;
        default: return nullptr;
    }
}

}

extern "C" {

double cspack_omega_min(double omega0, double p0)
{
    return CSpack_functions::omegamin(omega0, p0);
}

double cspack_omega_max(double omega0, double p0)
{
    return CSpack_functions::omegamax(omega0, p0);
}

double cspack_omega_crit(double omega0, double p0)
{
    return CSpack_functions::omegacrit(omega0, p0);
}

double cspack_kernel(int kernel_type, double omega0, double p0, double omega)
{
    switch(kernel_type)
    {
        case 0: return CSpack_kernels::kernel_exact(omega0, p0, omega);
        case 1: return CSpack_kernels::kernel_recoil(omega0, p0, omega);
        case 2: return CSpack_kernels::kernel_doppler(omega0, p0, omega);
        case 3: return CSpack_kernels::kernel_ur(omega0, p0, omega);
        default: return nan_value();
    }
}

double cspack_thermal_kernel(int kernel_type, double omega0, double omega, double theta)
{
    switch(kernel_type)
    {
        case 0: return CSpack_kernels::thermal_kernel_exact(omega0, omega, theta);
        case 1: return CSpack_kernels::thermal_kernel_exact_SS_C(omega0, omega, theta);
        case 2: return CSpack_kernels::thermal_kernel_recoil(omega0, omega, theta);
        case 3: return CSpack_kernels::thermal_kernel_doppler(omega0, omega, theta);
        case 4: return CSpack_kernels::thermal_kernel_ur(omega0, omega, theta);
        case 5: return CSpack_kernels::thermal_kernel_SS_K(omega0, omega, theta);
        case 6: return CSpack_kernels::thermal_kernel_SS_C(omega0, omega, theta);
        default: return nan_value();
    }
}

double cspack_neutrino_kernel(int kernel_type, double omega0, double p0, double omega)
{
    kernel_ptr K = neutrino_kernel_ptr(kernel_type);
    return K == nullptr ? nan_value() : K(omega0, p0, omega);
}

double cspack_thermal_photon_fd_kernel(int kernel_type, double omega0, double omega,
                                       double theta, double mue, int add_FB)
{
    kernel_ptr K = photon_kernel_ptr(kernel_type);
    return K == nullptr ? nan_value() : CSpack_kernels_FD::thermal_kernel_FD(
        omega0, omega, theta, mue, K, add_FB
    );
}

double cspack_thermal_neutrino_kernel(int kernel_type, double omega0, double omega,
                                      double theta, double mue, int add_FB)
{
    kernel_ptr K = neutrino_kernel_ptr(kernel_type);
    return K == nullptr ? nan_value() : CSpack_kernels_FD::thermal_kernel_FD(
        omega0, omega, theta, mue, K, add_FB
    );
}

double cspack_neutrino_kernel_norm(int kernel_type)
{
    switch(kernel_type)
    {
        case 0:
        case 1:
        case 2:
            return 1.0;
        case 3:
        case 4:
        case 5:
        case 6:
            return CSpack_kernels_nu::Get_neutrino_kernel_norm(
                "numu_e", "cspack_neutrino_kernel_norm"
            );
        default:
            return nan_value();
    }
}

}
