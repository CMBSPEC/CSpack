//==================================================================================================
//  Created by Abir Sarkar on 15/11/2019 and modified by JC. These functions are based on
//  Sarkar, Chluba and Lee, MNRAS, 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#include "routines.h"
#include "Patterson.h"

#include "CSpack.h"

using namespace std;
using namespace CSpack_functions;

//==================================================================================================
//
// COMPUTATION OF THE KERNEL
//
//==================================================================================================
double G_func(double omega, double omega0, double p0, double om0star, double omstar, double kappa)
{
    double common_fact = 3.0/(8.0*p0*gamma_f(p0)*pow(omega0,2));
    double lp=lambda_p(p0,omega0);
    double lm=lambda_m(p0,omega);

    double arg1 = pow(kappa/omstar ,2)*lp;
    double arg2 = pow(kappa/om0star,2)*lm;

    double fact1 = 2.0+pow((omstar-om0star)/(omega*omega0), 2)*(1.0+omega*omega0);
    double fact2 = sfunc(arg1)/ omstar     - sfunc(arg2)/ om0star;
    double fact3 = ffunc(arg1)/(omstar*lp) - ffunc(arg2)/(om0star*lm);

    return common_fact*kappa*(fact1+2.0*fact2+(1.0+omega*omega0)*fact3);
}

namespace CSpack_kernels {

//==================================================================================================
// EXACT KERNEL
//==================================================================================================
double kernel_exact(double omega0, double p0, double omega)
{
    double kerval=0.0;

    double omegab  = omegabar (omega0, p0, omega);
    double omegab0 = omega0bar(omega0, p0, omega);
    double p = pfunc_sc(omega0, p0, omega);
    double kappa1 = kappa(omega0, omega, p0, p);
    double kappa2 = kappa(omega, omega0, p, p0);

    double omega_min=omegamin(omega0,p0);
    double omega_cr=omegacrit(omega0,p0);
    double omega_max=omegamax(omega0,p0);

    if (omega0 > p0)
    {
        if(omega_min <= omega && omega < omega_cr)
            kerval = G_func(omega, omega0, p0, omegab0, omegab, kappa1);

        else if(omega_cr <= omega && omega <= omega0)
            kerval = G_func(omega, omega0, p0, omega, omega0, p0);

        else if(omega0 < omega && omega <= omega_max)
            kerval = G_func(omega, omega0, p0, omega0, omega, p);
    }
    else
    {
        if(omega_min <= omega && omega  < omega0)
            kerval = G_func(omega, omega0, p0, omegab0, omegab, kappa1);

        else if(omega0 <= omega && omega <=omega_cr)
            kerval = G_func(omega, omega0, p0, omegab, omegab0, kappa2);

        else if(omega_cr < omega && omega <= omega_max)
            kerval = G_func(omega, omega0, p0, omega0, omega, p);
    }

    return kerval;
}

//====================================================================================================================
// APPROXIMATE KERNELS
//====================================================================================================================
double kernel_recoil(double omega0, double p0, double omega)
{
    double delta = (omega0-omega)/omega;
    double kerval = 3.0/(8.0*pow(omega0,2))*(1.0 + pow(delta,2)/(1.0+delta) + pow((1.0-delta/omega0),2));
    return kerval;
}

double kernel_doppler(double omega0, double p0, double omega)
{
    double g0 = gamma_f(p0);
    double t = omega/omega0;
    double fact1 = (3.0+2.0*pow(p0, 2))/(2.0*p0)*(abs(log(t))-2.0*asinh(p0)) + (3.0+3.0*pow(p0, 2)+pow(p0, 4))/g0;
    double fact2 = 1.0+(10.0+8.0*pow(p0, 2)+4.0*pow(p0, 4))*t+pow(t,2);
    double kerval = 3.0/(8.0*omega0)*((1.0+t)/pow(p0, 5)*fact1 - abs(1.0-t)/(4.0*pow(p0, 6)*t)*fact2);
    return kerval;
}

double kernel_ur(double omega0, double p0, double omega)
{
    double g0 = gamma_f(p0);
    double G = 4.0*omega0*g0;
    double q = omega/(G*(g0-omega));
    double kerval = 3.0/(4.0*pow(g0,2)*omega0)*(2.0*q*log(q)+(1.0+2.0*q+pow(G*q, 2)/(2.0*(1.0 + G*q)))*(1.0-q));
    return kerval;
}

kernel_ptr Get_kernel_pointer(string type, string calling_func)
{
    if(type=="exact") return kernel_exact;
    else if(type=="recoil") return kernel_recoil;
    else if(type=="doppler") return kernel_doppler;
    else if(type=="ur") return kernel_ur;
    else throw_error(calling_func, "choose type 'exact', 'recoil', 'doppler', 'ur'", 1);

    return NULL;
}

double kernel_all(double omega0, double p0, double omega, string type)
{
    double (*kernel_ptr)(double omega0, double p0, double omega)=Get_kernel_pointer(type, "kernel_all");
    return kernel_ptr(omega0, p0, omega);
}

//==================================================================================================
//
// COMPUTATION OF THE THERMALLY AVERAGED KERNEL
//
//==================================================================================================

//==================================================================================================
// LIMITS OF INTEGRATION
//==================================================================================================
double lower_limit(double omega0, double omega)
{
    double lim=0.0;
    if(omega <= omega0)
    {
        if(omega <= omega0/(1.0+2.0*omega0))
            lim =  (omega0-omega)/2.0*sqrt((1.0+omega*omega0)/(omega*omega0)) - (omega0+omega)/2.0;

        else lim = 0.0;
    }

    else if(omega > omega0 && omega0 <= 0.5)
    {
        if(omega <= omega0/(1.0-2.0*omega0))
            lim = sqrt( (omega-omega0)*(omega-omega0+2.0) );

        else lim = (omega-omega0)/2.0*sqrt((1.0+omega*omega0)/(omega*omega0)) + (omega0+omega)/2.0;
    }

    else if(omega > omega0 && omega0 > 0.5)
    {
        lim = sqrt( (omega-omega0)*(omega-omega0+2.0) );
    }

    return lim;
}

double lower_limit_dop(double omega0, double omega)
{
    double lim=0.0;
    // Modified limits for Doppler-dominated kernel integrals
    if(omega <= omega0) lim = (omega0 - omega)/(2.0*sqrt(omega*omega0));
    else lim = (omega - omega0)/(2.0*sqrt(omega*omega0));

    return lim;
}

//==================================================================================================
// FOR ALL OF THE ABOVE IN ONE MODULE
//==================================================================================================
struct Integration_data
{
    double omega0, omega;
    double theta;
    double (*kernel_ptr)(double, double, double);
    
    Integration_data()
    {
        omega0=omega=theta=0.01;
        kernel_ptr=NULL;
    }
};

double integrand_p_all(double lp, void *q)
{
    Integration_data *d=(Integration_data *)q;
    double p=exp(lp);
    double fact= mb_dist_func(p, d->theta)*pow(p, 3) * d->kernel_ptr(d->omega0, p, d->omega);
    return fact;
}
 
double thermal_kernel_all(double omega0, double omega, double theta, string type)
{
    double epsrel=1.0e-8, epsabs=1.0e-100;
    double pmax = sqrt( (theta*log(1.0e-20) - 2.0)*theta*log(1.0e-20) );
 
    Integration_data d;
    d.omega=omega; d.omega0=omega0; d.theta=theta;
    d.kernel_ptr=Get_kernel_pointer(type, "thermal_kernel_all");

    double a, b;
    if (type=="doppler") a = lower_limit_dop(omega0, omega);
    else a = lower_limit(omega0, omega);

    a=max(sqrt(2.0*theta)*1.0e-8, a);
    b=max(pmax, a*10.0);
    
    double r=Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs, integrand_p_all, &d);
    return r/mb_dist_norm(theta);
}
 
//==================================================================================================
double thermal_kernel_exact(double omega0, double omega, double theta)
{
    return thermal_kernel_all(omega0, omega, theta, "exact");
}
 
double thermal_kernel_recoil(double omega0, double omega, double theta)
{
    return thermal_kernel_all(omega0, omega, theta, "recoil");
}

double thermal_kernel_doppler(double omega0, double omega, double theta)
{
    return thermal_kernel_all(omega0, omega, theta, "doppler");
}

double thermal_kernel_ur(double omega0, double omega, double theta)
{
    return thermal_kernel_all(omega0, omega, theta, "ur");
}

}

//==================================================================================================
//==================================================================================================
