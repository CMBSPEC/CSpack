//==================================================================================================
//  Created by Abir Sarkar on 15/11/2019 and modified by JC. These functions are based on
//  Sarkar, Chluba and Lee, MNRAS, 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#include "routines.h"

#include "CSpack.h"

using namespace std;
using namespace CSpack_functions;
using namespace CSpack_kernels;

namespace CSpack_scattering_matrix {

int verbosity_scat_matrix=0;

void set_verbosity(int verb){ verbosity_scat_matrix=verb; return; }

//==================================================================================================
// Routines for scattering matrix setups
//--------------------------------------------------------------------------------------------------
// inputs:
// xarr  : contains frequency grid points x=h nu/kTe = omega/theta
// theta : kTe/mc^2
// Int_wi: Integral weight factors to turn int f(x) dx == sum Int_wi f(xi)
// nK     : defines number of points per kernel wing for Kernel-representation method
//
// outputs: Msc = wj Pij theta
//          KR  = setup Kernel_representation vector on given grid and temperature
//
// epsilon: optional parameter to compress matrix density [eps<1.0e-4 recommended]
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
// epsilon: optional parameter to compress matrix density [eps<1.0e-4 recommended]
//==================================================================================================
void compute_scattering_matrix(const vector<double> &xarr, double theta,
                               const vector<double> &Int_wi,
                               vector<vector<double> > &Msc,
                               double epsilon)
{
    if(verbosity_scat_matrix>0)
        cout << " compute_scattering_matrix :: setting up scattering matrix for The= " << theta << endl;

    int npx=xarr.size();

    //double (*kernel)(double omega0, double omega, double theta)=thermal_kernel_exact;
    double (*kernel)(double omega0, double omega, double theta)=thermal_kernel_SS_K;
    //double (*kernel)(double omega0, double omega, double theta)=thermal_kernel_SS_C;

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
//==================================================================================================
