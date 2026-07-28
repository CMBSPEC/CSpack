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
double kernel_exact_nue_e(double omega0, double p0, double omega);
double kernel_exact_nue_p(double omega0, double p0, double omega);
double kernel_exact_nue_ep(double omega0, double p0, double omega);

double kernel_exact_numu_e(double omega0, double p0, double omega);
double kernel_exact_numu_p(double omega0, double p0, double omega);

double kernel_exact_nutau_e(double omega0, double p0, double omega);
double kernel_exact_nutau_p(double omega0, double p0, double omega);

};

//==================================================================================================
// thermally-averaged kernels over Fermi-Dirac distribution
//==================================================================================================
namespace CSpack_kernels_FD {

double norm_FD(double theta, double mue);

double thermal_kernel_FD(double omega0, double omega, double theta, double mue,
                         kernel_ptr K, int add_FB=0);

void output_thermal_kernel(string fname, int np,
                           vector<double> omega0,
                           double theta, double mue,
                           kernel_ptr K, int add_FB=0);
};

//==================================================================================================
// Evaluation of collision terms for thermal input distributions
//==================================================================================================
namespace CSpack_Collision_Terms {

double Collision_Term(double omega0, double th_r, double th_e, double mu_e, kernel_ptr K, string sel="nu");

// here the outer integral is over p
double Moment(double omega0, int k, double th_e, double mu_e, kernel_ptr K);

// here the inner integral is over p and the outer over omega3
double Moment_II(double omega0, int k, double th_e, double mu_e, kernel_ptr K);

double Moment_Doppler(int k, double th_e, double mu_e);
double Moment_recoil(double omega1, int k);

void output_Collision_Term(string fname,
                           double xmin, double xmax, int np,
                           double th_r, double mu_e,
                           string sel="ph");

};

//==================================================================================================
// Evaluation of opacities over moments
//==================================================================================================
namespace CSpack_opacity {

// integral over Sigma_k (T/me)^2 dlnT / H from T_i/me << 1 until theta=T/me
double tau_sc_generalized(double x, double theta, int k, kernel_ptr K);

// integral over Sigma_k (T/me)^2 dlnT / H from T_i/me << 1 until theta_max. A solution vector is returned
vector<vector<double>> tau_sc_generalized_ODE(double x, double theta_max, int k, kernel_ptr K);

// get x at which tau_sc==1
double x_tau_sc_equal_unity(double theta, int k, kernel_ptr K, double x_guess=0.0);

// get theta at which tau_sc==1
double theta_tau_sc_equal_unity(double x, int k, kernel_ptr K, double t_guess=0.0);

// get x at which tau_sc==1
double theta_tau_sc_equal_unity_ODE(double x, int k, kernel_ptr K, double theta_max);
};

//==================================================================================================
// simple Kernel outputs
//==================================================================================================
void output_kernel(string fname, int np,
                   double omega0, vector<double> p0,
                   kernel_ptr K);

void output_kernel(string fname, int np,
                   double omega0, double p0,
                   kernel_ptr K);

#endif
//==================================================================================================
//==================================================================================================
