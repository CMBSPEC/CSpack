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
    double sigma=sqrt(2.0*Theta); // <--> < (Dnu/nu)^2 >

    double lwstart=0.0, lwlim=log(omega_lim/omega0), lwsig=10.0*log(1.0+sigma);
    double lwc;

    if(lwlim<0.0) lwsig*=-1.0;
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
Kernel_representation::~Kernel_representation(){ spline_up=spline_down=-1; }

Kernel_representation::Kernel_representation(double omin, double om0, double omax,
                                             double The,
                                             double eps_thresh, double eps_interpol,
                                             int maxMom)
{
    int np=500;
    omega0=om0; Theta=The;
    spline_up=spline_down=-1;

    string type="exact";
    P0=thermal_kernel_all(omega0, omega0, Theta, type);

    create_kernel_splines(omin, eps_thresh, np, type);
    create_kernel_splines(omax, eps_thresh, np, type);

    // compute moments
    // Moments.push_back(m);
}

double Kernel_representation::Kernel(double om)
{
    double w=om/omega0;
    if(w==1.0) return P0;
    if(w>1.0) return exp(calc_spline_JC(log(w), spline_up, "P+ interpol"));
    return exp(calc_spline_JC(log(w), spline_down, "P- interpol"));;
}

//==================================================================================================
//==================================================================================================
