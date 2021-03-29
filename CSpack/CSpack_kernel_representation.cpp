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
Kernel_representation::~Kernel_representation()
{
    if(spline_up  !=-1) free_spline_JC(spline_up  , "P+");
    if(spline_down!=-1) free_spline_JC(spline_down, "P-");
}

Kernel_representation::Kernel_representation()
{
    spline_up=spline_down=np=-1;
}

Kernel_representation::Kernel_representation(double omin, double om0, double omax, int npv,
                                             double The,
                                             double eps_thresh, double eps_interpol,
                                             int maxMom, bool stim)
{
    spline_up=spline_down=np=-1;
    this->init(omin, om0, omax, npv, The, eps_thresh, eps_interpol, maxMom, stim);
}

//==================================================================================================
// for openmp runs this should be ran serial before the init call...
//==================================================================================================
void Kernel_representation::allocate_splines(int npv)
{
    if(spline_up  !=-1) free_spline_JC(spline_up  , "P+");
    if(spline_down!=-1) free_spline_JC(spline_down, "P-");

    np=npv;
    vector<double> lw(np), lP(np, 0.0);
    init_xarr(0.0, log(10.0), &lw[0], np, 0, 0);

    spline_up  =calc_spline_coeffies_JC(np, &lw[0], &lP[0], "P+");
    spline_down=calc_spline_coeffies_JC(np, &lw[0], &lP[0], "P-");

    return;
}

void Kernel_representation::init(double omin, double om0, double omax, int npv, double The,
                                 double eps_thresh, double eps_interpol, int maxMom, bool stim)
{
    omega0=om0; Theta=The;
    wmin=omin/omega0; wmax=omax/omega0;

    string type="exact";
    P0=thermal_kernel_all(omega0, omega0, Theta, type);

    if(np==-1) allocate_splines(npv);
    create_kernel_splines(omin, eps_thresh, npv, type);
    create_kernel_splines(omax, eps_thresh, npv, type);

    // compute moments
    for(int k=0; k<=maxMom; k++) Moments.push_back(compute_moment(k, stim));
    G=compute_G(stim);
    H=compute_H(stim);

//    for(int k=0; k<=maxMom; k++)
//        cout << k << " " << Moments[k]
//                  << " " << moment_2D_Int_therm_all_II(omega0, k, Theta, type)
//             << endl;
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

double find_root_CS(double (* func)(double *, void *p), void *p,
                    double x1, double x2, double xacc)
{
    if(x1==x2) return x1;
    
    double x;
    if(x1>x2) x=find_root_brent(func, p, x2, x1, xacc);
    else x=find_root_brent(func, p, x2, x1, xacc);

    return x;
}

//==================================================================================================
void Kernel_representation::create_kernel_splines(double omega_lim, double eps_thresh,
                                                  int npv, string type)
{
    if(np!=-1 && np!=npv) throw_error("create_kernel_splines", "memory not correctly allocated", 1);

    //==============================================================================================
    // no need to set up kernel part in this case
    //==============================================================================================
    if(omega_lim==omega0) return;

    //==============================================================================================
    // data for generation of splines
    //==============================================================================================
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

    if(P<=P0*eps_thresh) lwc=find_root_CS(root_func, &d, lwsig, lwstart, 1.0e-3);
    else
    {
        lwstart=lwsig;
        while(P>P0*eps_thresh && fabs(lwsig)<=fabs(lwlim))
        {
            lwsig*=2.0;
            P=thermal_kernel_all(omega0, omega0*exp(lwsig), Theta, type);
        }

        // if w found within range that brackets null --> solve
        if(P<P0*eps_thresh) lwc=find_root_CS(root_func, &d, lwsig, lwstart, 1.0e-3);
        else lwc=lwsig;
    }

    if(omega_lim<omega0) wmin=max(wmin, exp(lwc))/1.01;
    else wmax=min(wmax, exp(lwc))*1.01;

    //==============================================================================================
    // compute kernel in log-log
    //==============================================================================================
    vector<double> lw(np), lP(np);

    if(lwc>0.0)
    {
        init_xarr(0.0, log(wmax), &lw[0], np, 0, 0);
        lP[0]=log(P0);
        for(int i=1; i<np; i++)
            lP[i]=log(thermal_kernel_all(omega0, omega0*exp(lw[i]), Theta, type));

        // setup splines
        update_spline_coeffies_JC(spline_up, np, &lw[0], &lP[0], "P+");
    }
    else
    {
        init_xarr(log(wmin), 0.0, &lw[0], np, 0, 0);
        lP.back()=log(P0);
        for(int i=0; i<np-1; i++)
            lP[i]=log(thermal_kernel_all(omega0, omega0*exp(lw[i]), Theta, type));

        // setup splines
        update_spline_coeffies_JC(spline_down, np, &lw[0], &lP[0], "P-");
    }

    return;
}

//==================================================================================================
double Kernel_representation::Kernel(double om)
{
    double w=om/omega0;

    if(w<wmin || w>wmax) return 0.0;
    
    if(w>1.0 && spline_up!=-1) return exp(calc_spline_JC(log(w), spline_up, "P+ interpol"));
    if(w<1.0 && spline_down!=-1) return exp(calc_spline_JC(log(w), spline_down, "P- interpol"));

    return P0;
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

    double f=(d->stim ? one_minus_exp_mx(d->omega0/d->The)/one_minus_exp_mx(w*d->omega0/d->The) : 1.0);
    double Dnuk=pow(w-1.0, d->k);

    return w*Dnuk*exp(lP)*f;
}

double moment_func_flipped(double lw, void *p)
{
    momentData *d=(momentData *) p;
    double w=exp(lw);
    double lPp=( lw>=d->lwmax ? 0.0 : calc_spline_JC( lw, d->spline_up  , "P+ interpol") );
    double lPm=( lw>=d->lwmin ? 0.0 : calc_spline_JC(-lw, d->spline_down, "P- interpol") );

    double f=(d->stim ? one_minus_exp_mx(d->omega0/d->The)/one_minus_exp_mx(w*d->omega0/d->The) : 1.0);
    double Dnuk=pow(w-1.0, d->k);

    return w*Dnuk*( exp(lPp) + exp(lPm)/pow(-w, d->k+2) )*f;
}

double Kernel_representation::compute_moment(int k, bool stim)
{
    double epsrel=1.0e-8, epsabs=1.0e-100;

    momentData d;
    d.k=k;
    d.spline_up  =spline_up;
    d.spline_down=spline_down;
    d.lwmin=-log(wmin);
    d.lwmax= log(wmax);
    d.omega0=omega0; d.The=Theta;
    d.stim=stim;

    double r=Integrate_using_Patterson_adaptive(-d.lwmin, d.lwmax, epsrel, epsabs, moment_func, &d);
    
//    double r=Integrate_using_Patterson_adaptive(-d.lwmin, 0.0, epsrel, epsabs, moment_func, &d);
//    r+=Integrate_using_Patterson_adaptive(0.0, d.lwmax, epsrel, epsabs, moment_func, &d);

//    double r=Integrate_using_Patterson_adaptive(0.0, min(d.lwmin, d.lwmax),
//                                                epsrel, epsabs, moment_func_flipped, &d);
//
//    r+=Integrate_using_Patterson_adaptive(min(d.lwmin, d.lwmax), max(d.lwmin, d.lwmax),
//                                          epsrel, epsabs, moment_func_flipped, &d);

    return omega0*r;
}

double moment_G_func(double lw, void *p)
{
    // G=2 Sigma_2 - Sigma_1
    momentData *d=(momentData *) p;
    double w=exp(lw);
    double lP=( w>1.0 ? calc_spline_JC(lw, d->spline_up  , "P+ interpol")
                      : calc_spline_JC(lw, d->spline_down, "P- interpol") );

    double f=(d->stim ? one_minus_exp_mx(d->omega0/d->The)/one_minus_exp_mx(w*d->omega0/d->The) : 1.0);
    double Dnuk=(w-1.0)*(2.0*w-3.0);

    return w*Dnuk*exp(lP)*f;
}

double Kernel_representation::compute_G(bool stim)
{
    double epsrel=1.0e-9, epsabs=1.0e-100;

    momentData d;
    d.spline_up  =spline_up;
    d.spline_down=spline_down;
    d.lwmin=-log(wmin);
    d.lwmax= log(wmax);
    d.omega0=omega0; d.The=Theta;
    d.stim=stim;

    double r=Integrate_using_Patterson_adaptive(-d.lwmin, 0.0, epsrel, epsabs, moment_G_func, &d);
    r+=Integrate_using_Patterson_adaptive(0.0, d.lwmax, epsrel, epsabs, moment_G_func, &d);

    return omega0*r;
}

double moment_H_func(double lw, void *p)
{
    // H=4 Sigma_2 - Sigma_1
    momentData *d=(momentData *) p;
    double w=exp(lw);
    double lP=( w>1.0 ? calc_spline_JC(lw, d->spline_up  , "P+ interpol")
                      : calc_spline_JC(lw, d->spline_down, "P- interpol") );

    double f=(d->stim ? one_minus_exp_mx(d->omega0/d->The)/one_minus_exp_mx(w*d->omega0/d->The) : 1.0);
    double Dnuk=(w-1.0)*(4.0*w-5.0);

    return w*Dnuk*exp(lP)*f;
}

double Kernel_representation::compute_H(bool stim)
{
    double epsrel=1.0e-9, epsabs=1.0e-100;

    momentData d;
    d.spline_up  =spline_up;
    d.spline_down=spline_down;
    d.lwmin=-log(wmin);
    d.lwmax= log(wmax);
    d.omega0=omega0; d.The=Theta;
    d.stim=stim;

    double r=Integrate_using_Patterson_adaptive(-d.lwmin, 0.0, epsrel, epsabs, moment_H_func, &d);
    r+=Integrate_using_Patterson_adaptive(0.0, d.lwmax, epsrel, epsabs, moment_H_func, &d);

    return omega0*r;
}

//==================================================================================================
//
// Kernel representation class for multiple values of Theta
//
//==================================================================================================
void Kernel_representation_Te::init_Kernel_Table(double omin, double om0, double omax, int npom,
                                                 double Theta_min, double Theta_max, int logdens_The,
                                                 double eps_thresh, double eps_interpol,
                                                 int maxMom, bool stim)
{
    this->Theta_min =Theta_min;
    this->Theta_max =Theta_max;
    this->logdens_The=logdens_The;

    npThe=init_xarr_dens(Theta_min, Theta_max, The_arr, logdens_The, 1);

    Kernels_The.resize(npThe);

    // allocate memory for splines (done in serial)
    for(int k=0; k<npThe; k++) Kernels_The[k].allocate_splines(npom);

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int k=0; k<npThe; k++) Kernels_The[k].init(omin, om0, omax, npom, The_arr[k],
                                                   eps_thresh, eps_interpol, maxMom, stim);
}


Kernel_representation_Te::Kernel_representation_Te(double omin, double om0, double omax, int npom,
                                                   double Theta_min, double Theta_max, int logdens_The,
                                                   double eps_thresh, double eps_interpol,
                                                   int maxMom, bool stim)
{
    init_Kernel_Table(omin, om0, omax, npom,
                      Theta_min, Theta_max, logdens_The,
                      eps_thresh, eps_interpol, maxMom, stim);
}

//==================================================================================================
// workhorse for interpolation accross The
//==================================================================================================
double Kernel_representation_Te::do_interpol(double lgx, const double *lgxa, const double *ya)
{
    //===========================================================================
    // output with 4 point interpolation
    //===========================================================================
    double DD3=pow(lgxa[1]-lgxa[0], 3), r=0.0;

    for(int m=0; m<4; m++) D[m]=lgx-lgxa[m];

    a[0]= D[1]*D[2]*D[3]/(-6.0*DD3);
    a[1]= D[0]*D[2]*D[3]/( 2.0*DD3);
    a[2]= D[0]*D[1]*D[3]/(-2.0*DD3);
    a[3]= D[0]*D[1]*D[2]/( 6.0*DD3);

    for(int m=0; m<4; m++) r+=a[m]*ya[m];
    return r;
}

//==================================================================================================
// return kernel value for om0 -> om and The
//==================================================================================================
double Kernel_representation_Te::Kernel(double om, double The)
{
    if(npThe==0) throw_error("Kernel_representation_Te::Kernel", "kernel data not set", 1);

    bool show=0;
    bool rangeok=check_xrange(The, The_arr, show, "The");
    if(!rangeok) throw_error("Kernel_representation_Te::Kernel", "The outside of range", 2);

    // find index around The
    unsigned int intj=get_start_index_interpol(The, The_arr, 4);

    // store interpolation data
    for(int ix=0; ix<4; ix++)
    {
        lgThe[ix]=log(The_arr[intj+ix]);
        lgK  [ix]=log(Kernels_The[intj+ix].Kernel(om));
    }

    return exp(do_interpol(log(The), lgThe, lgK));
}

//==================================================================================================
//==================================================================================================
