//==================================================================================================
//
// Integral weights using Lagrange interpolation polynomial
//
//==================================================================================================
// Integral int g(x) dx is written as sum_i wi gi. The function g(x) is expressed in terms of
// Lagrange polynomials. An extra weigt function f(x) can be introduced to write
// int f(x) g(x) dx = sum_i wfi gi.
//==================================================================================================
//
// Author: Jens Chluba
// first implementation: June 2021
// last modification   : June 2021
//
//==================================================================================================

#ifndef CSpack_weights_h
#define CSpack_weights_h

#include <string>
#include <vector>

#include "CSpack.h"

using namespace std;

namespace CSpack_weights {

void set_verbosity(int verb=0);

//==================================================================================================
// integral weights for scattering matrix
//==================================================================================================
void Integral_weights_trapz(vector<double> &xarr, vector<double> &Int_wi);
void Integral_weights(vector<double> &xarr, vector<double> &Int_wi);
void Integral_weights_logx(vector<double> &xarr, vector<double> &Int_wi);

//==================================================================================================
// explicit computation of weights
//==================================================================================================
void compute_all_Lagrange_Polynomial_weights(int order, const vector<double> &xa, vector<double> &w);

//==================================================================================================
// compute weights using extra weight function
//==================================================================================================
void Lagrange_Polynomial_weights(int i, int order,
                                 const vector<double> &xa, vector<double> &w,
                                 double (*fw)(double x, void *pfw), void *pfw,
                                 bool do_norm=1);

void compute_all_Lagrange_Polynomial_weights(int order, const vector<double> &xa, vector<double> &w,
                                             double (*fw)(double x, void *pfw), void *pfw=NULL);
}

#endif
//==================================================================================================
//==================================================================================================
