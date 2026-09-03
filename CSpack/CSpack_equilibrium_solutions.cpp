//==================================================================================================
// Created by Jens Chluba March 2020.
// Last modified: Sept 2026 [JC+Codex]
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

    return find_root_brent(root_func_equilibrium, &d, -0.1, 200.0, 1.0e-16,
                           "Compute_muc_N");
}

double Compute_muc_rho(double DTg_Te, double Drho_rho)
{
    if(DTg_Te<=-1.0) throw_error("Compute_muc_rho", "inconsistent values", 1);
    Root_equilibrium_sols d;
    d.DTg_Te=DTg_Te;
    d.Drho_rho=Drho_rho;
    d.number=0;

    return find_root_brent(root_func_equilibrium, &d, -0.1, 200.0, 1.0e-16,
                           "Compute_muc_rho");
}

double Compute_Drho_rho(double DTg_Te, double DN_N)
{
    double muc=Compute_muc_N(DTg_Te, DN_N);
    double DE_E_int=DEnergy_integral(muc);
    double DTT4=(2.0+DTg_Te)*DTg_Te*(2.0+DTg_Te*(2.0+DTg_Te));

    return (DE_E_int-DTT4)/(1.0+DTT4);
}

//==================================================================================================
double root_func_equilibrium_rhoN(double *muc, void *q)
{
    Root_equilibrium_sols *d=(Root_equilibrium_sols *)q;

    double rN=DNumber_integral(*muc), rrho=DEnergy_integral(*muc);
    double er=d->Drho_rho*(3.0+d->Drho_rho*(3.0+d->Drho_rho));
    double eN=d->DN_N*(2.0+d->DN_N)*(2.0+d->DN_N*(2.0+d->DN_N));
    double DN=rN*(2.0+rN)*(2.0+rN*(2.0+rN));
    double Dr=rrho*(3.0+rrho*(3.0+rrho));

    double r=er-eN;
    return r+DN*(1.0+er)-Dr*(1.0+eN);
}

double Compute_muc_rhoN(double Drho_rho, double DN_N)
{
    Root_equilibrium_sols d;
    d.Drho_rho=Drho_rho;
    d.DN_N=DN_N;

    return find_root_brent(root_func_equilibrium_rhoN, &d, 0.0, 10.0, 1.0e-12,
                           "Compute_muc_rhoN");
}

double Compute_DT_T(double Drho_rho, double DN_N, double muc)
{
    double rN=DNumber_integral(muc), rrho=DEnergy_integral(muc);

    return ( Drho_rho+rN*(1.0+rrho)-(DN_N+rrho*(1.0+DN_N)) )/(1.0+DN_N)/(1.0+rrho);
}

//==================================================================================================
double root_func_equilibrium_DT_T(double *DT_T, void *q)
{
    Root_equilibrium_sols *d=(Root_equilibrium_sols *)q;
    double del2=*DT_T;

    return del2*(3.0+del2*(3.0+del2))-d->DN_N;
}

double Compute_DT_T(double Drho_rho, double DN_N)
{
    Root_equilibrium_sols d;

    double muc=Compute_muc_rhoN(Drho_rho, DN_N);
    double rN=DNumber_integral(muc);
    d.DN_N=(rN-DN_N)/(1.0+DN_N);

    return find_root_brent(root_func_equilibrium_DT_T, &d, -200.0, 200.0, 1.0e-16,
                           "Compute_DT_T");
}

double Compute_Drho_rho_inj(double DTg_Te, double muc)
{
    double rrho=DEnergy_integral(muc);

    return ( rrho-DTg_Te*(2.0+DTg_Te)*(2.0+DTg_Te*(2.0+DTg_Te)) )/pow(1.0+DTg_Te, 4);
}

double Compute_DN_N_inj(double DTg_Te, double muc)
{
    double rN=DNumber_integral(muc);

    return ( rN-DTg_Te*(3.0+DTg_Te*(3.0+DTg_Te)) )/pow(1.0+DTg_Te, 3);
}

}

//==================================================================================================
//==================================================================================================
