//==================================================================================================
// Created by Abir Sarkar on 15/11/2019 and modified by JC.
// Last modified: Sept 2026 [JC+Codex]
// These functions are based on Sarkar, Chluba and Lee, MNRAS, 2019
// (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#ifndef CSpack_functions_hpp
#define CSpack_functions_hpp

namespace CSpack_functions {

//==================================================================================================
// Common simple functions
//==================================================================================================
double gamma_f(double p0);
double pfunc(double g0);
double gamma_sc(double omega0, double p0, double omega);
double pfunc_sc(double omega0, double p0, double omega);

//==================================================================================================
// Critical frequencies
//==================================================================================================
double omegacrit(double omega0, double p0);
double omegamin(double omega0, double p0);
double omegatot(double omega0, double p0);
double omegamax(double omega0, double p0);

string Get_zone(double omega0, double p0, double omega);

//==================================================================================================
// Kernel functions
//==================================================================================================
double lambda_p(double p0, double omega0);
double lambda_m(double p0, double omega);
double sfunc(double x);
double ffunc(double x);
double omegabar(double omega0, double p0, double omega);
double omega0bar(double omega0, double p0, double omega);
double kappa(double omega0,  double omega, double p0, double p);

//==================================================================================================
// relativistic Maxwell-Boltzmann distribution function
//==================================================================================================
double mb_dist_func(double p, double theta);
double mb_dist_norm(double theta);
double rel_mb_dist(double p, double theta);
double pbar(double theta);

//==================================================================================================
// functions for Moments
//==================================================================================================
double alphap(double omega0, double p0);
double alpham(double omega0, double p0);
double f_moment(double omega0, double p0);
}

#endif
//==================================================================================================
//==================================================================================================
