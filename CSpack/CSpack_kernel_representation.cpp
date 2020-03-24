//==================================================================================================
//  Created by JC in March 2020.
//==================================================================================================

#include <string>
#include <vector>

#include "routines.h"
#include "Patterson.h"

#include "CSpack.h"

using namespace std;
using namespace CSpack_functions;
using namespace CSpack_kernels;
using namespace CSpack_kernel_moments;

//==================================================================================================
Kernel_representation::~Kernel_representation(){ spline_up=spline_down=-1; }

Kernel_representation::Kernel_representation(double omin, double om0, double omax, int np,
                                             double The,
                                             double eps_thresh, double eps_interpol,
                                             int maxMom)
{
    omega0=om0; Theta=The;
    spline_up=spline_down=-1;

    string type="exact";
    P0=thermal_kernel_all(omega0, omega0, Theta, type);

    create_kernel_splines(omin, eps_thresh, np, type);
    create_kernel_splines(omax, eps_thresh, np, type);

    // compute moments
    for(int k=0; k<=maxMom; k++)
        Moments.push_back(compute_moment(k));

    for(int k=0; k<=maxMom; k++) cout << k << " " << Moments[k] << " "
                                      << moment_2D_Int_therm_all_II(omega0, k, Theta, type) << endl;
}

//==================================================================================================
// function for root finding process
//==================================================================================================
struct rootData
{
    double P0eps, omega0, Theta;
    string type;
};

double root_func(double *lw, void *p)
{
    rootData *d=(rootData *) p;
    double w=exp(*lw);
    return thermal_kernel_all(d->omega0, d->omega0*w, d->Theta, d->type)/d->P0eps-1.0;
}

//==================================================================================================
void Kernel_representation::create_kernel_splines(double omega_lim, double eps_thresh,
                                                  int np, string type)
{
    rootData d;
    d.P0eps=P0*eps_thresh;
    d.omega0=omega0; d.Theta=Theta;
    d.type=type;

    //==============================================================================================
    // estimate of kernel width
    //==============================================================================================
    double wfac=( omega_lim<omega0 ? omegamin(omega0, pbar(Theta))/omega0
                                   : omegamax(omega0, pbar(Theta))/omega0 );

    double lwstart=0.0, lwlim=log(omega_lim/omega0), lwsig=log(wfac), lwc;
    double P=thermal_kernel_all(omega0, omega0*exp(lwsig), Theta, type);

    if(lwlim<0.0) // omega =< omega0
    {
        if(P<=P0*eps_thresh) lwc=find_root_brent(root_func, &d, lwsig, lwstart, 1.0e-3);
        else
        {
            lwstart=lwsig;
            while(P>P0*eps_thresh && lwsig>=lwlim)
            {
                lwsig*=2.0;
                P=thermal_kernel_all(omega0, omega0*exp(lwsig), Theta, type);
            }

            lwsig=max(lwlim, lwsig);
            lwc=find_root_brent(root_func, &d, lwsig, lwstart, 1.0e-3);
        }

        wmin=exp(lwc);
        //cout << " min " << omegamin(omega0, pbar(Theta))/omega0 << " " << wmin << endl;
    }
    else  // omega >= omega0
    {
        if(P<=P0*eps_thresh) lwc=find_root_brent(root_func, &d, lwstart, lwsig, 1.0e-3);
        else
        {
            lwstart=lwsig;
            while(P>P0*eps_thresh && lwsig<=lwlim)
            {
                lwsig*=2.0;
                P=thermal_kernel_all(omega0, omega0*exp(lwsig), Theta, type);
            }

            lwsig=min(lwlim, lwsig);
            lwc=find_root_brent(root_func, &d, lwstart, lwsig, 1.0e-3);
        }

        wmax=exp(lwc);
        //cout << " max " << omegamax(omega0, pbar(Theta))/omega0 << " " << wmax << endl;
    }

    //==============================================================================================
    // compute kernel in log-log
    //==============================================================================================
    vector<double> lw(np), lP(np);

    if(lwc>0.0)
    {
        init_xarr(0.0, lwc, &lw[0], np, 0, 0);
        lP[0]=log(P0);
        for(int i=1; i<np; i++)
            lP[i]=log(thermal_kernel_all(omega0, omega0*exp(lw[i]), Theta, type));

        // setup splines
        if(spline_up==-1) spline_up=calc_spline_coeffies_JC(np, &lw[0], &lP[0], "P+");
        else update_spline_coeffies_JC(spline_up, np, &lw[0], &lP[0], "P+");
    }
    else
    {
        init_xarr(lwc, 0.0, &lw[0], np, 0, 0);
        lP.back()=log(P0);
        for(int i=0; i<np-1; i++)
            lP[i]=log(thermal_kernel_all(omega0, omega0*exp(lw[i]), Theta, type));

        // setup splines
        if(spline_down==-1) spline_down=calc_spline_coeffies_JC(np, &lw[0], &lP[0], "P-");
        else update_spline_coeffies_JC(spline_down, np, &lw[0], &lP[0], "P-");
    }

    return;
}

//==================================================================================================
double Kernel_representation::Kernel(double om)
{
    double w=om/omega0;
    if(w<=wmin || w>=wmax) return 0.0;
    if(w==1.0) return P0;
    if(w>1.0) return exp(calc_spline_JC(log(w), spline_up, "P+ interpol"));
    return exp(calc_spline_JC(log(w), spline_down, "P- interpol"));
}

//==================================================================================================
struct momentData
{
    int k;
    int spline_up, spline_down;
    bool stim;
    double omega0, The;
    double lwmin, lwmax;

    momentData() { stim =0; }
};

double moment_func(double lw, void *p)
{
    momentData *d=(momentData *) p;
    double w=exp(lw);
    double lP=( w>1.0 ? calc_spline_JC(lw, d->spline_up  , "P+ interpol")
                      : calc_spline_JC(lw, d->spline_down, "P- interpol") );

    double Dnuk=pow(w-1.0, d->k);

    return w*Dnuk*exp(lP);
}

double moment_func_flipped(double lw, void *p)
{
    momentData *d=(momentData *) p;
    double w=exp(lw);
    double lPp=( lw>=d->lwmax ? 0.0 : calc_spline_JC( lw, d->spline_up  , "P+ interpol") );
    double lPm=( lw>=d->lwmin ? 0.0 : calc_spline_JC(-lw, d->spline_down, "P- interpol") );

    double Dnuk=pow(w-1.0, d->k);

    return w*Dnuk*( exp(lPp) + exp(lPm)/pow(-w, d->k+2) );
}

double Kernel_representation::compute_moment(int k)
{
    double epsrel=1.0e-8, epsabs=1.0e-100;

    momentData d;
    d.k=k;
    d.spline_up  =spline_up;
    d.spline_down=spline_down;
    d.lwmin=-log(wmin);
    d.lwmax= log(wmax);

    double r=Integrate_using_Patterson_adaptive(-d.lwmin, 0.0, epsrel, epsabs, moment_func, &d);
    r+=Integrate_using_Patterson_adaptive(0.0, d.lwmax, epsrel, epsabs, moment_func, &d);

//    double r=Integrate_using_Patterson_adaptive(0.0, min(d.lwmin, d.lwmax),
//                                                epsrel, epsabs, moment_func_flipped, &d);
//
//    r+=Integrate_using_Patterson_adaptive(min(d.lwmin, d.lwmax), max(d.lwmin, d.lwmax),
//                                          epsrel, epsabs, moment_func_flipped, &d);

    return omega0*r;
}


//==================================================================================================
//==================================================================================================
