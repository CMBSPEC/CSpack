//==================================================================================================
//  Created by Abir Sarkar on 15/11/2019 and modified by JC. These functions are based on
//  Sarkar, Chluba and Lee, MNRAS, 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#ifndef cspack_h
#define cspack_h

//==================================================================================================
// available C functions
//==================================================================================================

//==================================================================================================
// Planck function, written in such a way that is always gives very accurate results
//==================================================================================================
double nbb_func_C(double x);                                                      // == 1/(exp(x)-1)

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
void Integral_weights(const double *xarr, double *Int_wi, int npointsx);

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
                                        double **Msc);

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
                                  double **Msc);
#endif

//==================================================================================================
//==================================================================================================
