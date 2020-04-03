//==================================================================================================
//  Created by Abir Sarkar on 15/11/2019 and modified by JC. These functions are based on
//  Sarkar, Chluba and Lee, MNRAS, 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#ifndef CSpack_scattering_matrix_h
#define CSpack_scattering_matrix_h

#include <string>
#include <vector>

#include "CSpack_kernel_representation.h"

using namespace std;

namespace CSpack_scattering_matrix {

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
                               double epsilon=1.0e-40, bool stim=0);

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
                               double epsilon=1.0e-40);

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
                       bool add_stim=0,
                       double Te_Tg=1.0);

}

#endif
//==================================================================================================
//==================================================================================================
