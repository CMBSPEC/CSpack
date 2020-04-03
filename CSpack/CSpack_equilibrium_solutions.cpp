//==================================================================================================
//  Created by Jens Chluba March 2020.
//==================================================================================================

#include "physical_consts.h"
#include "routines.h"
#include "Definitions.h"
#include "Patterson.h"

#include "CSpack.h"

using namespace std;
using namespace CSpack_functions;

namespace CSpack_equilibrium_solutions {

//==================================================================================================
// compute moments numerically
//==================================================================================================
struct Integration_equilibrium_sols
{
    double muc;      // constant chemical potential for equilibrium solution
    int alpha;

    Integration_equilibrium_sols()
    {
        alpha=2;
        muc=0.0;
    }
};

//==================================================================================================
double integrand_equilibrium(double lgx, void *q)
{
    Integration_equilibrium_sols *d=(Integration_equilibrium_sols *)q;

    double x=exp(lgx);
    // x^alpha * f= x^alpha * ( 1/[exp(x+muc)-1] - 1/[exp(x)-1] )
    return -pow(x, d->alpha+1)*nbb_func(x)*one_minus_exp_mx(d->muc)/one_minus_exp_mx(x+d->muc);
}

//==================================================================================================
double DNumber_integral(double muc)
{
    Integration_equilibrium_sols d;
    d.muc=muc;
    d.alpha=2;

    double a=log(1.0e-16), b=log(200.0);
    double epsrel=1.0e-8, epsabs=1.0e-100;
    double r=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, integrand_equilibrium, &d);

    return r/G21_Int_pl;
}

double Number_integral(double muc){ return 1.0+DNumber_integral(muc); }

double DEnergy_integral(double muc)
{
    Integration_equilibrium_sols d;
    d.muc=muc;
    d.alpha=3;

    double a=log(1.0e-16), b=log(200.0);
    double epsrel=1.0e-8, epsabs=1.0e-100;
    double r=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, integrand_equilibrium, &d);

    return r/G31_Int_pl;
}

double Energy_integral(double muc){ return 1.0+DEnergy_integral(muc); }

//==================================================================================================
struct Root_equilibrium_sols
{
    double DTg_Te;
    double DN_N, Drho_rho;
    bool number;

    Root_equilibrium_sols()
    {
        DTg_Te=1.0;
        DN_N=Drho_rho=0.0;
        number=1;
    }
};

double root_func_equilibrium(double *muc, void *q)
{
    Root_equilibrium_sols *d=(Root_equilibrium_sols *)q;

    double Dphi=d->DTg_Te;

    double r=( d->number ? DNumber_integral(*muc) : DEnergy_integral(*muc) );

    if(d->number) return Dphi*(3.0+Dphi*(3.0+Dphi))-(r-d->DN_N    )/(1.0+d->DN_N    );
    return    (2.0+Dphi)*Dphi*(2.0+Dphi*(2.0+Dphi))-(r-d->Drho_rho)/(1.0+d->Drho_rho);
}

double Compute_muc_N(double DTg_Te, double DN_N)
{
    if(DTg_Te<=-1.0) throw_error("Compute_muc_N", "inconsistent values", 1);
    Root_equilibrium_sols d;
    d.DTg_Te=DTg_Te;
    d.DN_N=DN_N;
    d.number=1;

    return find_root_brent(root_func_equilibrium, &d, 0.0, 200.0, 1.0e-16);
}

double Compute_muc_rho(double DTg_Te, double Drho_rho)
{
    if(DTg_Te<=-1.0) throw_error("Compute_muc_rho", "inconsistent values", 1);
    Root_equilibrium_sols d;
    d.DTg_Te=DTg_Te;
    d.Drho_rho=Drho_rho;
    d.number=0;

    return find_root_brent(root_func_equilibrium, &d, 0.0, 200.0, 1.0e-16);
}

}

//==================================================================================================
//==================================================================================================
