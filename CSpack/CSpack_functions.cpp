//==================================================================================================
// Created by Abir Sarkar on 15/11/2019 and modified by JC. These routines are based on:
// Sarkar, Chluba and Lee, MNRAS, 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>

#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_bessel.h>

#include "Definitions.h"
#include "routines.h"
#include "CSpack.h"

using namespace std;

namespace CSpack_functions {

//==================================================================================================
// Common simple functions
//==================================================================================================
double gamma_f(double p0)
{ return sqrt(pow(p0, 2)+1.0); }

double pfunc(double g0) // sqrt(g^2-1)
{
    double Dg=g0-1.0;
    return sqrt(Dg*(2.0+Dg));
}

//==================================================================================================
double gamma_sc(double omega0, double p0, double omega)
{ return gamma_f(p0)+omega0-omega; }

double pfunc_sc(double omega0, double p0, double omega)
{ return sqrt( pow(p0, 2) + ((omega0-omega) + 2.0*gamma_f(p0))*(omega0-omega) ); }

//==================================================================================================
// Critical frequencies
//==================================================================================================
double omegacrit(double omega0, double p0)
{ return omega0*pow(gamma_f(p0) + p0, 2)/(1.0 + 2.0*omega0*(gamma_f(p0) + p0)); }

double omegamin(double omega0, double p0)
{ return omega0/(gamma_f(p0) + p0)/(gamma_f(p0) + p0 + 2.0*omega0); }

double omegatot(double omega0, double p0)
{ return p0*p0/(1.0+gamma_f(p0)) + omega0; }

double omegamax(double omega0, double p0)
{
    //if(omega0 > 0.5*(1.0 + p0 - gamma_f(p0))) return omegatot(omega0, p0);
    if(omega0 > 0.5*(1.0 - 1.0/(gamma_f(p0)+p0))) return omegatot(omega0, p0);

    return omegacrit(omega0, p0);
}

//==================================================================================================
string Get_zone(double omega0, double p0, double omega)
{
    double omega_min=omegamin(omega0, p0);
    double omega_cr=omegacrit(omega0, p0);
    double omega_max=omegamax(omega0, p0);

    if (omega0 > p0)
    {
        if(omega_min <= omega && omega < omega_cr) return "I";
        
        else if(omega_cr <= omega && omega <= omega0) return "II";
        
        else if(omega0 < omega && omega <= omega_max) return "III";
    }
    else
    {
        if(omega_min <= omega && omega  < omega0) return "I";

        else if(omega0 <= omega && omega <=omega_cr) return "II";

        else if(omega_cr < omega && omega <= omega_max) return "III";
    }

    return "None";
}

//==================================================================================================
// Functions for the Kernel
//==================================================================================================
double lambda_p(double p0, double omega0)
{ return p0*p0+2.0*gamma_f(p0)*omega0+omega0*omega0; }

double lambda_m(double p0, double omega)
{ return p0*p0-2.0*gamma_f(p0)*omega+omega*omega; }

//==================================================================================================
double sfunc_exact(double x)
{
    if (x>0) return asinh(sqrt(x))/sqrt(x);

    return asin(sqrt(-x))/sqrt(-x);
}

double ffunc_exact(double x)
{ return sfunc_exact(x)-sqrt(1.0+x); }

//==================================================================================================
double sfunc(double x)
{ return (abs(x) <= 1.0e-4 ? 1.0-x/6.0+3.0/40.0*pow(x, 2)-5.0/112.0*pow(x, 3) : sfunc_exact(x) ); }

double ffunc(double x)
{ return (abs(x) <= 1.0e-4 ? -2.0/3.0*x+pow(x, 2)/5.0-3.0/28.0*pow(x, 3) : ffunc_exact(x) ); }

//==================================================================================================
double omegabar(double omega0, double p0, double omega)
{
    double p=pfunc_sc(omega0, p0, omega);
    return sqrt( omega*omega0*(p+gamma_f(p))/(gamma_f(p0)+p0) );
}

double omega0bar(double omega0, double p0, double omega)
{
    double p=pfunc_sc(omega0, p0, omega);
    return sqrt( omega*omega0*(gamma_f(p0)+p0)/(p+gamma_f(p)) );
}

double kappa(double omega0, double omega, double p0, double p)
{ return (p0-p+omega0+omega)/2.0; }

//==================================================================================================
// relativistic Maxwell-Boltzmann distribution function
//==================================================================================================
double mb_dist_func(double p, double theta)
{
    double gam = gamma_f(p), func;
    func = exp(-p*p/(1.0+gam)/theta);
    return func;
}

double mb_dist_norm(double the)
{
    if(the>=0.005) return exp(1.0/the)*the*gsl_sf_bessel_Kn(2, 1.0/the);
    
    return sqrt(PI/2.0)*pow(the, 1.5)*(1.0 + the*(1.875 + the*(0.8203125
                                           + the*(-0.3076171875 + the*(0.317230224609375
                                           + the*(-0.5154991149902344 + 1.1276543140411377*the))))));
}

double rel_mb_dist(double p, double theta)
{ return mb_dist_func(p, theta)/mb_dist_norm(theta); }

double pbar(double theta)
{
    return 2.0*pow(theta, 2)*(1.0 + 3.0*theta + 3*pow(theta, 2))/mb_dist_norm(theta);
}

//==================================================================================================
// functions for Moments
//==================================================================================================
double alphap(double omega0, double p0)
{ return  1.0 + 2.0*(gamma_f(p0) + p0)*omega0; }

double alpham(double omega0, double p0)
{ return  1.0 + 2.0/(gamma_f(p0) + p0)*omega0; }

double f_moment(double omega0, double p0)
{ return gsl_sf_dilog(-2.0*(gamma_f(p0)+p0)*omega0)-gsl_sf_dilog(-2.0/(gamma_f(p0)+p0)*omega0); }

}

//==================================================================================================
extern "C" {

double nbb_func_C(double x){ return nbb_func(x); }      // == 1/[exp(x)-1]

}

//==================================================================================================
//==================================================================================================
