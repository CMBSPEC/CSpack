//==================================================================================================
//  CSpack_energy_losses.h
//
//  Created by Jens Chluba March 2020.
//==================================================================================================

#include <gsl/gsl_sf_bessel.h>

#include "physical_consts.h"
#include "routines.h"
#include "Definitions.h"
#include "Patterson.h"

#include "CSpack.h"

using namespace std;
using namespace CSpack_functions;
using namespace CSpack_kernels;
using namespace CSpack_kernel_moments;

namespace CSpack_energy_losses {

//==================================================================================================
//
//  losses of photon scattering off thermal electrons with ambient blackbody radiation (Tg=Te)
//
//==================================================================================================
// d(ln rho_h)/dtau == Sigma_1^*
//==================================================================================================
double photon_energy_loss(double omega0, double the, string type, bool stim)
{ return moment_2D_Int_therm_all_II(omega0, 1, the, type, stim); }

//==================================================================================================
//
//  losses of electron scattering off ambient blackbody radiation at temperature Tg
//
//==================================================================================================

//==================================================================================================
// compute moments numerically
//==================================================================================================
struct Integration_2Ddata_losses
{
    double omega0, p0;      // omega0 and p0
    double thg;             // temperature of black body in me c^2
    double (*kernel_ptr)(double, double, double);
    bool use_stim;

    Integration_2Ddata_losses()
    {
        kernel_ptr=NULL;
        use_stim=0;
    }
};

//==================================================================================================
double integrand_losses_all_o(double lgomega, void *q)
{
    Integration_2Ddata_losses *d=(Integration_2Ddata_losses *)q;

    double omega=exp(lgomega);
    double K=d->kernel_ptr(d->omega0, d->p0, omega);
    double Dnu_nuk=omega/d->omega0-1.0;
    double stim=(d->use_stim ? 1.0/one_minus_exp_mx(omega/d->thg) : 1.0);

    return omega*K*Dnu_nuk*stim;
}

double integrand_losses_all_o0(double lgomega0, void *q)
{
    Integration_2Ddata_losses *d=(Integration_2Ddata_losses *)q;
    d->omega0=exp(lgomega0);

    double a=max(1.0e-16, omegamin(d->omega0, d->p0));
    double b=omegamax(d->omega0, d->p0);

    double epsrel=1.0e-8, epsabs=1.0e-100;
    return pow(d->omega0, 4)*nbb_func(d->omega0/d->thg)
                            *Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs,
                                                                integrand_losses_all_o, &(*d));
}

//==================================================================================================
double photon_gains(double thg, double p0, string type, bool stim)
{
    Integration_2Ddata_losses d;
    d.thg=thg;
    d.p0=p0;
    d.use_stim=stim;
    d.kernel_ptr=Get_kernel_pointer(type, "photon_gains");

    double a=1.0e-5*thg, b=50.0*thg;
    double epsrel=1.0e-8, epsabs=1.0e-50;
    double r=Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs,
                                                integrand_losses_all_o0, &d);

    double rho_CMB=pow(thg, 4)*G31_Int_pl;
    return r/rho_CMB;
}

double photon_gains_approx(double thg, double p0)
{ return 4.0/3.0*p0*p0-3.83223*thg; }

double electron_cooling(double thg, double p0, string type, bool stim)
{
    double gamma0=gamma_f(p0);
    double z=thg/const_kb_mec2/2.7255-1.0;
    double Ne=2.511e-7*(1.0-0.24)*pow(1.0+z, 3);
    double rhog_mec2=8.0*PI/pow(const_lambdac, 3)*pow(thg, 4)*G31_Int_pl;
    return photon_gains(thg, p0, type, stim) * (gamma0+1.0)/p0/p0 * rhog_mec2/Ne;
}

//==================================================================================================
double integrand_Ng_removal_o(double lgomega, void *q)
{
    Integration_2Ddata_losses *d=(Integration_2Ddata_losses *)q;

    double omega=exp(lgomega);
    double K=d->kernel_ptr(d->omega0, d->p0, omega);
    double stim=(d->use_stim ? 1.0/one_minus_exp_mx(omega/d->thg) : 1.0);

    return -omega*K*stim;
}

double integrand_Ng_removal_o0(double lgomega0, void *q)
{
    Integration_2Ddata_losses *d=(Integration_2Ddata_losses *)q;
    d->omega0=exp(lgomega0);

    double a=max(1.0e-16, omegamin(d->omega0, d->p0));
    double b=omegamax(d->omega0, d->p0);

    double epsrel=1.0e-8, epsabs=1.0e-100;
    return pow(d->omega0, 3)*nbb_func(d->omega0/d->thg)
                            *Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs,
                                                                integrand_Ng_removal_o, &(*d));
}

//==================================================================================================
double Ng_removal(double thg, double p0, string type, bool stim)
{
    Integration_2Ddata_losses d;
    d.thg=thg;
    d.p0=p0;
    d.use_stim=stim;
    d.kernel_ptr=Get_kernel_pointer(type, "Ng_removal");

    double a=1.0e-5*thg, b=50.0*thg;
    double epsrel=1.0e-8, epsabs=1.0e-50;
    double r=Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs,
                                                integrand_Ng_removal_o0, &d);

    double N_CMB=pow(thg, 3)*G21_Int_pl;
    return r/N_CMB;
}

}

//==================================================================================================
//==================================================================================================
