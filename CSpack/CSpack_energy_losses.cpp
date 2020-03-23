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
    bool removal, addition;

    Integration_2Ddata_losses()
    {
        kernel_ptr=NULL;
        use_stim=0;
        removal=addition=1;
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
double integrand_Ng_all_o(double lgomega, void *q)
{
    Integration_2Ddata_losses *d=(Integration_2Ddata_losses *)q;

    double omega=exp(lgomega);
    double Km=(d->removal  ? d->kernel_ptr(d->omega0, d->p0, omega) : 0.0);
    double Kp=(d->addition ? d->kernel_ptr(omega, d->p0, d->omega0) : 0.0);

    double x0=d->omega0/d->thg, x=omega/d->thg;
    double stimm=(d->use_stim ? 1.0/one_minus_exp_mx(x ) : 1.0);
    double stimp=(d->use_stim ? 1.0/one_minus_exp_mx(x0) : 1.0);

    return omega*( pow(x/x0, 2)*Kp*nbb_func(x)*stimp-Km*stimm*nbb_func(x0) );
}

double integrand_Ng_all_o0(double lgomega0, void *q)
{
    Integration_2Ddata_losses *d=(Integration_2Ddata_losses *)q;
    d->omega0=exp(lgomega0);

    double a=max(1.0e-16, omegamin(d->omega0, d->p0));
    double b=omegamax(d->omega0, d->p0);

    double epsrel=1.0e-8, epsabs=1.0e-100;
    return pow(d->omega0, 3)*Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs,
                                                                integrand_Ng_all_o, &(*d));
}

//==================================================================================================
double dDng_dtau(double x, double thg, double p0, string type, bool stim)
{
    Integration_2Ddata_losses d;
    d.thg=thg;
    d.p0=p0;
    d.use_stim=stim;
    d.kernel_ptr=Get_kernel_pointer(type, "Ng_removal");
    d.removal=1;
    d.addition=1;

    double omega0=x*thg;

    return integrand_Ng_all_o0(log(omega0), &d)/pow(omega0, 3);
}

//==================================================================================================
double dng_dtau_removal(double x, double thg, double p0, string type, bool stim)
{
    Integration_2Ddata_losses d;
    d.thg=thg;
    d.p0=p0;
    d.use_stim=stim;
    d.kernel_ptr=Get_kernel_pointer(type, "Ng_removal");
    d.removal=1;
    d.addition=0;

    double omega0=x*thg;

    return integrand_Ng_all_o0(log(omega0), &d)/pow(omega0, 3);
}

//==================================================================================================
double DNg_dtau(double thg, double p0, string type, bool stim)
{
    Integration_2Ddata_losses d;
    d.thg=thg;
    d.p0=p0;
    d.use_stim=stim;
    d.kernel_ptr=Get_kernel_pointer(type, "Ng_removal");
    d.removal=1;
    d.addition=1;

    double a=1.0e-5*thg, b=200.0*thg;
    double epsrel=1.0e-8, epsabs=1.0e-50;
    double r=Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs,
                                                integrand_Ng_all_o0, &d);

    double N_CMB=pow(thg, 3)*G21_Int_pl;
    return r/N_CMB;
}

//==================================================================================================
double Ng_dtau_removal(double thg, double p0, string type, bool stim)
{
    Integration_2Ddata_losses d;
    d.thg=thg;
    d.p0=p0;
    d.use_stim=stim;
    d.kernel_ptr=Get_kernel_pointer(type, "Ng_removal");
    d.removal=1;
    d.addition=0;

    double a=1.0e-5*thg, b=200.0*thg;
    double epsrel=1.0e-8, epsabs=1.0e-50;
    double r=Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs,
                                                integrand_Ng_all_o0, &d);

    double N_CMB=pow(thg, 3)*G21_Int_pl;
    return r/N_CMB;
}

//==================================================================================================
struct Integration_data_DC
{
    double omega0, p, g;
    double lam0, lam1, phi1;
    double Thg;

    Integration_data_DC() {}
};

double integrand_GDC_phi1lam0lam1(double lam1, void *q)
{
    Integration_data_DC *d=(Integration_data_DC *)q;
    d->lam1=lam1;

    // angular functions
    double mu0 =(d->g-d->lam0)/d->p;
    double mu1 =(d->g-d->lam1)/d->p;
    double mu01=mu0*mu1 + cos(d->phi1) * sin(1.0-mu0*mu0) * sin(1.0-mu1*mu1);
    double a01=1.0-mu01;

    double xi=d->omega0*a01/d->lam1;
    double zeta=d->lam0*d->omega0*xi/(1.0+xi);
    double Del=1.0+zeta;
    double p0=sqrt(zeta*(2.0+zeta));
    double F=Del/p0 * log(Del+p0)-1.0;
    double w=a01/(d->lam1*d->lam0);

    return d->lam0 * F/pow(d->lam1*(1.0+xi), 2) * ( 1.0/(1.0+xi)+(1.0+xi) + w*(w-2.0) );
}

double integrand_GDC_phi1lam0(double lam0, void *q)
{
    Integration_data_DC *d=(Integration_data_DC *)q;
    d->lam0=lam0;

    double epsrel=1.0e-6, epsabs=1.0e-100;
    return Integrate_using_Patterson_adaptive(d->g-d->p, d->g+d->p, epsrel, epsabs,
                                              integrand_GDC_phi1lam0lam1, &(*d));
}

double integrand_GDC_phi1(double phi1, void *q)
{
    Integration_data_DC *d=(Integration_data_DC *)q;
    d->phi1=phi1;

    double epsrel=3.0e-6, epsabs=1.0e-100;
    return Integrate_using_Patterson_adaptive(d->g-d->p, d->g+d->p, epsrel, epsabs,
                                              integrand_GDC_phi1lam0, &(*d));
}

double integrand_GDC(double omega0, double p0)
{
    Integration_data_DC d;
    d.omega0=omega0;
    d.p=p0;
    d.g=gamma_f(p0);

    double normf=2.0*9.0/(64.0*PI)/pow(omega0*p0, 2)/d.g;

    double epsrel=1.0e-5, epsabs=1.0e-100;
    return normf*Integrate_using_Patterson_adaptive(0.0, PI, epsrel, epsabs,
                                                    integrand_GDC_phi1, &d);
}

double integrand_GDC_nplphi1(double lgx, void *q)
{
    Integration_data_DC *d=(Integration_data_DC *)q;
    double x=exp(lgx), omega=x*d->Thg;

    double G=integrand_GDC(omega, d->p);

    return G*pow(x, 5)*nbb_func(x)/24.8863;
}

double integrand_GDC_npl(double Thg, double p0)
{
    Integration_data_DC d;
    d.p=p0;
    d.g=gamma_f(p0);
    d.Thg=Thg;

    double epsrel=1.0e-4, epsabs=1.0e-100;
    return Integrate_using_Patterson_adaptive(log(1.0e-4), log(1.0e+2), epsrel, epsabs,
                                              integrand_GDC_nplphi1, &d);
}

//==================================================================================================
double dng_dtau_DC_add(double x, double thg, double p0, string type, bool stim)
{
    double fstim=(stim ? 1.0/one_minus_exp_mx(x) : 1.0);
    return 4.0*const_alpha/(3.0*PI)*(1.0+2.0*p0*p0)*24.8863*thg*thg*fstim/pow(x, 3);
}

double Ng_DC_add(double thg, double p0, string type, bool stim)
{
    double xmin=0.001, xmax=200.0;
    double fstim=(stim ? log(xmax/xmin)+1.0/xmin-1.0/xmax : log(xmax/xmin));

    double Gtherm=integrand_GDC_npl(thg, p0);
    double r=4.0*const_alpha/(3.0*PI)*24.8863*thg*thg*Gtherm*fstim;
    return r/G21_Int_pl;
}

double Ng_DC_add_approx(double thg, double p0, string type, bool stim)
{
    double xmin=0.001, xmax=200.0;
    double fstim=(stim ? log(xmax/xmin)+1.0/xmin-1.0/xmax : log(xmax/xmin));

    double Gtherm=(1.0+2.0*p0*p0);
    double r=4.0*const_alpha/(3.0*PI)*24.8863*thg*thg*Gtherm*fstim;
    return r/G21_Int_pl;
}


}

//==================================================================================================
//==================================================================================================
