//
//  CSpack_functions.cpp
//  
//
//  Created by Abir Sarkar on 15/11/2019.
//

#include "CSpack_functions.h"
#include "routines.h"
#include "physical_consts.h"

using namespace std;
//======================================================================================
// General Functions
//======================================================================================
//======================================================================================
// Critical Frequencies
//======================================================================================

//double omegacrit(double omega0, double p0)
//{ return omega0*(gamma_ei(p0) + p0)/(gamma_ei(p0) - p0 + 2.0*omega0); }
//
//double omegamin(double omega0, double p0)
//{ return omega0*(gamma_ei(p0) - p0)/(gamma_ei(p0) + p0 + 2.0*omega0); }

double omegacrit(double omega0, double p0)
{ return omega0*pow(gamma_ei(p0) + p0, 2)/(1.0 + 2.0*omega0*(gamma_ei(p0) + p0)); }

double omegamin(double omega0, double p0)
{ return omega0/(gamma_ei(p0) + p0 + 2.0*omega0)/(gamma_ei(p0) + p0); }

double omegamax(double omega0, double p0)
{
    double ommax;
    if(omega0 > 0.5*(1.0 + p0 - gamma_ei(p0)))
    {
        ommax= -1.0 + gamma_ei(p0) + omega0;
    }
    else
    {
        ommax = omegacrit(omega0, p0);
    }
    
    return ommax;
    
}

//======================================================================================
// Common functions
//======================================================================================

double gamma_ei(double p0)
{
    return sqrt(pow(p0, 2)+1.0);
    
}

double beta(double p0)

{ return p0/gamma_ei(p0);}

double p_ef(double g)
{
    return sqrt(pow(g, 2)-1.0);
    
}

double gamma_ef(double omega0, double p0, double omega)
{
    return gamma_ei(p0)+omega0-omega;
}

double pfunc_low(double omega0, double pval, double omega)
{
    return sqrt(pow(pval, 2) + ((omega0-omega) + 2.0*gamma_ei(pval))*(omega0-omega));
}

//======================================================================================
// Functions for the Kernel
//======================================================================================

double lambda_p(double p0, double omega0)
{
    return pow((gamma_ei(p0)+omega0),2)-1.0;
}

double lambda_m(double p0, double omega)
{
    return pow((gamma_ei(p0)-omega),2)-1.0;
}

double sfunc(double x)
{
    double res;
    if (x>0)
        res = asinh(sqrt(x))/sqrt(x);
    else
        res = asin(sqrt(-x))/sqrt(-x);
    return res;
}

double ffunc(double x)
{
    return sfunc(x)-sqrt(1.0+x);
}


double sfunc_con(double x)
{
    return (abs(x) <= 1.0e-04 ? 1.0 - x/6.0 + 3.0/40.0*pow(x, 2) - 5.0/112.0*pow(x, 3) : sfunc(x) );
}

double ffunc_con(double x)
{
    return (abs(x) <= 1.0e-04 ? -2.0/3.0*x + pow(x, 2)/5.0 - 3.0/28.0*pow(x, 3) : ffunc(x) );
}

double omegabar(double omega0, double p0, double omega)
{
    double fact = (omega*omega0*(pfunc_low(omega0, p0, omega)+gamma_ei(pfunc_low(omega0, p0, omega))))/(gamma_ei(p0)+p0);
    return sqrt(fact);
}

double omega0bar(double omega0, double p0, double omega)
{
    double fact = (omega*omega0*(gamma_ei(p0)+p0))/(pfunc_low(omega0, p0, omega)+gamma_ei(pfunc_low(omega0, p0, omega)));
    return sqrt(fact);
}

double kappa(double omega0, double omega, double p0, double p)
{
    return  0.5*(p0-p+omega0+omega);
}

//======================================================================================
// Thermally Averaged Kernel
//======================================================================================
double mb_dist_func(double p, double theta)
{
    double g0 = gamma_ei(p), func;
    func = exp(-p*p/(1.0+g0)/theta);
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
{
    return mb_dist_func(p, theta)/mb_dist_norm(theta);
}

//======================================================================================
// Moment
//======================================================================================

double alphap(double omega0, double pval)
{
    return  1.0 + 2.0*(1.0 + beta(pval))*gamma_ei(pval)*omega0;
}

double alpham(double omega0, double pval)
{
    return  1.0 + 2.0*(1.0 - beta(pval))*gamma_ei(pval)*omega0;
}


double f_moment(double omega0, double pval)
{
    return  gsl_sf_dilog(1.0 - alphap(omega0, pval)) - gsl_sf_dilog(1.0 - alpham(omega0, pval));
}


