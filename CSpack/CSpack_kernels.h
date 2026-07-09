//==================================================================================================
//  Created by Abir Sarkar on 15/11/2019 and modified by JC. These functions are based on
//  Sarkar, Chluba and Lee, MNRAS, 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#ifndef CSpack_kernels_h
#define CSpack_kernels_h

#include <string>

using namespace std;

typedef double (*kernel_ptr)(double, double, double);

namespace CSpack_kernels {

//==================================================================================================
// KERNELS
//==================================================================================================
double kernel_exact(double omega0, double p0, double omega);
double kernel_recoil(double omega0, double p0, double omega);
double kernel_doppler(double omega0, double p0, double omega);
double kernel_ur(double omega0, double p0, double omega);

//--------------------------------------------------------------------------------------------------
//type == 'exact', 'recoil', 'doppler', 'ur'
//--------------------------------------------------------------------------------------------------
double kernel_all(double omega0, double p0, double omega, string type);

kernel_ptr Get_kernel_pointer(string type, string calling_func);

//==================================================================================================
//  THERMALLY AVERAGED KERNELS
//==================================================================================================
double thermal_kernel_exact(double omega0, double omega, double theta);
double thermal_kernel_recoil(double omega0, double omega, double theta);
double thermal_kernel_doppler(double omega0, double omega, double theta);
double thermal_kernel_ur(double omega0, double omega, double theta);

//==================================================================================================
// kernels from Sazonov & Sunyaev 2000
//==================================================================================================
double thermal_kernel_SS_K(double omega0, double omega, double theta);
double thermal_kernel_SS_C(double omega0, double omega, double theta);

//--------------------------------------------------------------------------------------------------
// type == 'exact', 'recoil', 'doppler', 'ur', 'SS_K', 'SS_C'
//--------------------------------------------------------------------------------------------------
double thermal_kernel_all(double omega0, double omega, double theta, string type);

}

namespace CSpack_kernels_nu {

//==================================================================================================
// neutrino scattering KERNELS
//==================================================================================================
double kernel_exact(double omega0, double p0, double omega);

};

void output_kernel(string fname, int np,
                   double omega0, vector<double> p0,
                   kernel_ptr K);

void output_kernel(string fname, int np,
                   double omega0, double p0,
                   kernel_ptr K);

#endif
//==================================================================================================
//==================================================================================================
