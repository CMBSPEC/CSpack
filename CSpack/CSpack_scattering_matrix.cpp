//==================================================================================================
//  Created by Abir Sarkar on 15/11/2019 and modified by JC. These functions are based on
//  Sarkar, Chluba and Lee, MNRAS, 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#include "routines.h"

#include "CSpack.h"

using namespace std;
using namespace CSpack_functions;
using namespace CSpack_kernels;

//==================================================================================================
// main functions in namespace
//==================================================================================================
namespace CSpack_scattering_matrix {

int verbosity_scat_matrix=0;

void set_verbosity(int verb){ verbosity_scat_matrix=verb; return; }

//==================================================================================================
// weights using f(x) = L(x) [5 point Lagrange polynomial] and integrating over
// x in [xa, xb] around xj
//==================================================================================================
double Int_li_xj_II(double x1, double x2, double x3, double x4, double xj,
                    double xa, double xb)
{
    double r=(pow(xb, 5)-pow(xa, 5))/5.0;
    r+=-(x1+x2+x3+x4)/4.0 * (pow(xb, 4)-pow(xa, 4));
    r+=(x1*(x2+x3+x4)+x2*(x3+x4)+x3*x4)/3.0 * (pow(xb, 3)-pow(xa, 3));
    r+=-(x1*(x2*x3+x2*x4+x3*x4)+x2*x3*x4)/2.0 * (xb-xa) * (xb+xa);
    r+= x1*x2*x3*x4 * (xb-xa);

    return r /(x1-xj)/(x2-xj)/(x3-xj)/(x4-xj);
}

double Int_li_xj_II(int j, int k, vector<double> &xi)
{
    //    double xa=(xi[j]+xi[j-1])/2.0, xb=(xi[j+1]+xi[j])/2.0;
    double xa=xi[j], xb=xi[j+1];
    double x1, x2, x3, x4;

    if(k==-2){ x1=xi[j-1]; x2=xi[j]; x3=xi[j+1]; x4=xi[j+2]; }
    else if(k==-1){ x1=xi[j-2]; x2=xi[j]; x3=xi[j+1]; x4=xi[j+2]; }
    else if(k== 0){ x1=xi[j-2]; x2=xi[j-1]; x3=xi[j+1]; x4=xi[j+2]; }
    else if(k== 1){ x1=xi[j-2]; x2=xi[j-1]; x3=xi[j]; x4=xi[j+2]; }
    else if(k== 2){ x1=xi[j-2]; x2=xi[j-1]; x3=xi[j]; x4=xi[j+1]; }
    else { return 0.0; }

    double r=Int_li_xj_II(x1, x2, x3, x4, xi[j+k], xa, xb);

    return r;
}

//==================================================================================================
void Integral_weights_trapz(vector<double> &xarr, vector<double> &Int_wi)
{
    int np=xarr.size();
    Int_wi.resize(np, 0.0);

    //================================================================
    // initial points around lower boundary
    //================================================================
    int ix=0;
    Int_wi[ix]+=(xarr[ix+1]-xarr[ix])/2.0; ix++;
    for(; ix<np-1; ix++) Int_wi[ix]+=(xarr[ix+1]-xarr[ix-1])/2.0;
    Int_wi[ix]+=(xarr[ix]-xarr[ix-1])/2.0;

    return;
}

//==================================================================================================
void Integral_weights(vector<double> &xarr, vector<double> &Int_wi)
{
    int np=xarr.size();
    Int_wi.resize(np, 0.0);

    //================================================================
    // initial points around lower boundary
    //================================================================
    int ix=0;
    Int_wi[ix]+=(xarr[ix+1]-xarr[ix])/2.0; ix++;
    for(; ix<2; ix++) Int_wi[ix]+=(xarr[ix+1]-xarr[ix-1])/2.0;
    Int_wi[ix]+=(xarr[ix]-xarr[ix-1])/2.0;

    //================================================================
    // internal points with 5-point formula
    //================================================================
    for(ix=2; ix<np-3; ix++)
        for(int k=-2; k<3; k++) Int_wi[ix+k]+=Int_li_xj_II(ix, k, xarr);

    //================================================================
    // finish off using simple trapeziodal rule
    //================================================================
    Int_wi[ix]+=(xarr[ix+1]-xarr[ix])/2.0; ix++;
    for(; ix<np-1; ix++) Int_wi[ix]+=(xarr[ix+1]-xarr[ix-1])/2.0;
    Int_wi[ix]+=(xarr[ix]-xarr[ix-1])/2.0;

    return;
}

//==================================================================================================
// Routines for scattering matrix setups
//--------------------------------------------------------------------------------------------------
// inputs :
// xarr   : contains frequency grid points x=h nu/kTe = omega/theta
// theta  : kTe/mc^2
// Int_wi : Integral weight factors to turn int f(x) dx == sum Int_wi f(xi)
// nK     : defines number of points per kernel wing for Kernel-representation method
//
// outputs: Msc = wj Pij theta
//          KR  = setup Kernel_representation vector on given grid and temperature
//
// epsilon: optional parameter to compress matrix density [eps<1.0e-4 recommended]
// stim   : include stimulated factors from blackbody in moments
//==================================================================================================
void compute_scattering_matrix(const vector<double> &xarr, double theta,
                               const vector<double> &Int_wi, int nK,
                               vector<vector<double> > &Msc,
                               vector<Kernel_representation> &KR,
                               double epsilon, bool stim)
{
    if(verbosity_scat_matrix>0)
        cout << " compute_scattering_matrix :: setting up scattering matrix for The= " << theta << endl;

    int npx=xarr.size();

    // create matrix
    if((int)Msc.size()!=npx)
    {
        Msc.clear();
        vector<double> zeros(npx, 0.0);
        for(int i=0; i<npx; i++) Msc.push_back(zeros);
    }

    // make vector of Kernel representations
    KR.resize(npx);
    double omin=xarr[0]*theta/300.0, omax=xarr.back()*theta*300.0;
    //double omin=xarr[0]*theta, omax=xarr.back()*theta;
    for(int i=0; i<npx; i++) KR[i].allocate_splines(nK);

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int i=0; i<npx; i++)
        KR[i].init(omin, xarr[i]*theta, omax, nK, theta, epsilon, epsilon, 2, stim);

    for(int i=0; i<npx; i++) //x
    {
        // reset matrix
        for(int j=0; j<npx; j++) Msc[i][j]=0.0;

        double omega_fac=Int_wi[i] * theta; // dnu' weight
        double om0=xarr[i]*theta;
        // diagonal element for reference
        Msc[i][i]= KR[i].Kernel(om0) * omega_fac;

        for(int j=i+1; j<npx; j++) //xp>x
        {
            double omega_fac=Int_wi[j] * theta; // dnu' weight
            double omp=xarr[j]*theta;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc[i][j]= KR[i].Kernel(omp) * omega_fac;

            if(abs(Msc[i][j]/Msc[i][i])<epsilon) break;
        }

        for(int j=i-1; j>=0; j--) //xp<x
        {
            double omega_fac=Int_wi[j] * theta; // dnu' weight
            double omp=xarr[j]*theta;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc[i][j]= KR[i].Kernel(omp) * omega_fac;

            if(abs(Msc[i][j]/Msc[i][i])<epsilon) break;
        }
    }

    if(verbosity_scat_matrix>0)
        cout << " compute_scattering_matrix :: done." << endl << endl;
}

//==================================================================================================
// Routines for scattering matrix setups
//--------------------------------------------------------------------------------------------------
// inputs:
// xarr  : contains frequency grid points x=h nu/kTe = omega/theta
// theta : kTe/mc^2
// Int_wi: Integral weight factors to turn int f(x) dx == sum Int_wi f(xi)
//
// outputs: Msc = wj Pij theta
//
// type   : type of kernel to be used explicitly ['exact', 'SS_K', 'SS_C']
// epsilon: optional parameter to compress matrix density [eps<1.0e-4 recommended]
//==================================================================================================
void compute_scattering_matrix(const vector<double> &xarr, double theta,
                               const vector<double> &Int_wi,
                               vector<vector<double> > &Msc,
                               string type,
                               double epsilon)
{
    if(verbosity_scat_matrix>0)
        cout << " compute_scattering_matrix :: setting up scattering matrix for The= " << theta << endl;

    int npx=xarr.size();

    double (*kernel)(double omega0, double omega, double theta);

    if(type=="exact") kernel=thermal_kernel_exact;
    else if(type=="SS_K") kernel=thermal_kernel_SS_K;
    else if(type=="SS_C") kernel=thermal_kernel_SS_C;
    else throw_error("compute_scattering_matrix", "kernel type not available", 1);

    // create matrix
    if((int)Msc.size()!=npx)
    {
        Msc.clear();
        vector<double> zeros(npx, 0.0);
        for(int i=0; i<npx; i++) Msc.push_back(zeros);
    }

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int i=0; i<npx; i++) //x
    {
        // reset matrix
        for(int j=0; j<npx; j++) Msc[i][j]=0.0;

        double omega_fac=Int_wi[i] * theta; // dnu' weight
        double om0=xarr[i]*theta;
        // diagonal element for reference
        Msc[i][i]= kernel(om0, om0, theta) * omega_fac;

        for(int j=i+1; j<npx; j++) //xp>x
        {
            double omega_fac=Int_wi[j] * theta; // dnu' weight
            double omp=xarr[j]*theta;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc[i][j]= kernel(om0, omp, theta) * omega_fac;

            if(abs(Msc[i][j]/Msc[i][i])<epsilon) break;
        }

        for(int j=i-1; j>=0; j--) //xp<x
        {
            double omega_fac=Int_wi[j] * theta; // dnu' weight
            double omp=xarr[j]*theta;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc[i][j]= kernel(om0, omp, theta) * omega_fac;

            if(abs(Msc[i][j]/Msc[i][i])<epsilon) break;
        }
    }
#ifdef OPENMP_ACTIVATED
#pragma omp barrier
#endif

    if(verbosity_scat_matrix>0)
        cout << " compute_scattering_matrix :: done." << endl << endl;
}

//==================================================================================================
// total cross section
//--------------------------------------------------------------------------------------------------
// inputs:
// xarr  : contains frequency grid points x = hnu/kTe = omega/theta for which Msc is defined
// Msc   : corresponding scattering matrix
//
// outputs: sigarr_i = sum_j Msc_ij * stim
//
// add_stim: optional parameter to add stimulated scattering effect in blackbody radiation field
// Te_Tg   : allow for difference between blackbody and electron temperature
//==================================================================================================
void compute_sigma_tot(const vector<double> &xarr,
                       const vector<vector<double> > &Msc,
                       vector<double> &sigarr,
                       bool add_stim,
                       double Te_Tg)
{
    int npx=xarr.size();

    if((int)sigarr.size()!=npx) sigarr.resize(npx);

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int i=0; i<npx; i++)
    {
        sigarr[i]=0.0;

        for(int j=0; j<npx; j++)
        {
            double stim=(add_stim ? one_minus_exp_mx(xarr[i]*Te_Tg)/one_minus_exp_mx(xarr[j]*Te_Tg) : 1.0);
            sigarr[i]+= Msc[i][j] * stim;
        }
    }

    return;
}

}

//==================================================================================================
// global memory that is required for C part of the routines
//==================================================================================================
vector<vector<double> > _global_Msc_temp;
vector<Kernel_representation> _global_KR_temp;

//==================================================================================================
// C versions of the scattering matrix setup
//==================================================================================================
extern "C" {

//==================================================================================================
// computes weights to turn int f(x) dx == sum Int_wi f(xi) on the grid given by xarr
//--------------------------------------------------------------------------------------------------
// inputs  :
// *xarr   : pointer to frequency grid points x=h nu/kTe = omega/theta
// npointsx: points in x
//--------------------------------------------------------------------------------------------------
// output  :
// *Int_wi : pointer to Integral weight factors to turn int f(x) dx == sum Int_wi f(xi).
//--------------------------------------------------------------------------------------------------
// comment : Memory has to be allocated before calling the function
//==================================================================================================
void Integral_weights(const double *xarr, double *Int_wi, int npointsx)
{
    vector<double> vecx(npointsx), vecw;
    for(int i=0; i<npointsx; i++) vecx[i]=xarr[i];

    //CSpack_scattering_matrix::Integral_weights_trapz(vecx, vecw); // simplest rule
    CSpack_scattering_matrix::Integral_weights(vecx, vecw); // 5 points rule [more accurate]

    for(int i=0; i<npointsx; i++) Int_wi[i]=vecw[i];

    return;
}

//==================================================================================================
// Routine for scattering matrix setup using explicit computation but with threshold
//--------------------------------------------------------------------------------------------------
// inputs  :
// *xarr   : pointer to frequency grid points x=h nu/kTe = omega/theta
// *Int_wi : pointer to Integral weight factors to turn int f(x) dx == sum Int_wi f(xi)
// npointsx: points in x
// theta   : kTe/mc^2
// type    : type of kernel to be used explicitly [0: 'exact', 1: 'SS_K', 2: 'SS_C']
// epsilon : optional parameter to compress matrix density [eps<1.0e-4 recommended]
//--------------------------------------------------------------------------------------------------
// outputs : Msc = wj Pij theta
//--------------------------------------------------------------------------------------------------
// comment : Memory has to be allocated before calling the function
//==================================================================================================
void compute_scattering_matrix_explicit(const double *xarr, const double *Int_wi, int npointsx,
                                        double theta, int kernel_type, double epsilon,
                                        double **Msc)
{
    vector<double> vecx(npointsx), vecw(npointsx);
    for(int i=0; i<npointsx; i++){ vecx[i]=xarr[i]; vecw[i]=Int_wi[i]; }

    string type;
    if(kernel_type==0) type="exact";
    else if(kernel_type==1) type="SS_K";
    else if(kernel_type==2) type="SS_C";

    CSpack_scattering_matrix::compute_scattering_matrix(vecx, theta, vecw,
                                                        _global_Msc_temp,
                                                        type, epsilon);

    for(int i=0; i<npointsx; i++)
        for(int j=0; j<npointsx; j++) Msc[i][j]=_global_Msc_temp[i][j];

    return;
}

//==================================================================================================
// Routine for scattering matrix setup using kernel representation routines
//--------------------------------------------------------------------------------------------------
// inputs  :
// *xarr   : pointer to frequency grid points x=h nu/kTe = omega/theta
// *Int_wi : pointer to Integral weight factors to turn int f(x) dx == sum Int_wi f(xi)
// npointsx: points in x
// nK      : defines number of points per kernel wing for Kernel-representation method
// theta   : kTe/mc^2
// epsilon : optional parameter to compress matrix density [eps<1.0e-4 recommended]
//--------------------------------------------------------------------------------------------------
// outputs : Msc = wj Pij theta
//--------------------------------------------------------------------------------------------------
// comment : Memory has to be allocated before calling the function
//==================================================================================================
void compute_scattering_matrix_KR(const double *xarr, const double *Int_wi, int npointsx,
                                  int nK, double theta, double epsilon,
                                  double **Msc)
{
    vector<double> vecx(npointsx), vecw(npointsx);
    for(int i=0; i<npointsx; i++){ vecx[i]=xarr[i]; vecw[i]=Int_wi[i]; }

    CSpack_scattering_matrix::compute_scattering_matrix(vecx, theta, vecw, nK,
                                                        _global_Msc_temp, _global_KR_temp,
                                                        epsilon, 0);

    for(int i=0; i<npointsx; i++)
        for(int j=0; j<npointsx; j++) Msc[i][j]=_global_Msc_temp[i][j];

    return;
}

}
//==================================================================================================
//==================================================================================================
