//==================================================================================================
//  CSpack_kernel_moments.h
//
//  Created by Abir Sarkar on 15/11/2019 and modified by JC. These functions are based on
//  Sarkar, Chluba and Lee, MNRAS, 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#ifndef CSpack_kernel_moments_h
#define CSpack_kernel_moments_h

#include <string>

using namespace std;

namespace CSpack_kernel_moments {

//==================================================================================================
//  MOMENT
//==================================================================================================
double compute_exact_moments_analytical(double omega0, double p0, int l);
double compute_nr_moments_analytical(double omega0, double p0, int l);
double compute_ur_moments_analytical(double omega0, double p0, int l);
double compute_recoil_moments_analytical(double omega0, double p0, int l);
double compute_doppler_moments_analytical(double omega0, double p0, int l);

//==================================================================================================
// THERMALLY AVERAGED MOMENT
//==================================================================================================
double moment_Int_therm(double omega0, int l, double theta);
double moment_Int_app1(double omega0, int l, double theta);
double moment_Int_app2(double omega0, int l, double theta);
double moment_Int_app3(double omega0, int l, double theta);

//==================================================================================================
// THERMALLY AVERAGED MOMENTS using 2D integration
//==================================================================================================
double moment_2D_Int_therm_all(double omega0, int k, double theta, string type); // JC: not as precise
double moment_2D_Int_therm_all_II(double omega0, int k, double theta, string type, bool stim=0);

//==================================================================================================
// obtain coefficients for FP approximation dn/dy= A(x)*d^2n/dx^2 + B(x)*dn/dx + C(x)*n
//==================================================================================================
void compute_all_FP_coefficients(double omega0, double theta, string type,
                                 double &A, double &B, double &C);

// 2Sigma2-Sigma1
double G_moment_2D_Int_therm_all_II(double omega0, double theta, string type);

}

#endif

//==================================================================================================
//==================================================================================================
