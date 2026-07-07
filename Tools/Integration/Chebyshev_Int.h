//==================================================================================================
// Author: Jens Chluba
// first implementation: Jan 2010
// last modification   : Aug 2021
//==================================================================================================
// Aug  5 th, 2021: made integrator backward compatible
// Mar  30th, 2020: generalized integrator to allow for parameters.
// Oct  20th, 2011: One can now give pmax (refinement level) to the integrator.
// June 14th, 2011: Nested table will only be updated if more abscissae and weights are needed.

#ifndef CHEBYSHEV_INT_H
#define CHEBYSHEV_INT_H

#include <iostream>
#include <string>

using namespace std;

void no_Chebyshev_alert();
void Chebyshev_alert();

//==================================================================================================
//
// some simple test integrands
//
//==================================================================================================
double Chebyshev_test_fcn(double x, void *pars);
double Chebyshev_test_fcn_log(double lgx, void *pars);

//==================================================================================================
// setup nested integral weight table before actual computations
//==================================================================================================
void setup_nested_Chebyshev_table(int pmax, int mess=0);

//==================================================================================================
// Driver that loops over the different Chebyshev rules
//
// a : lower bound of integral
// b : upper bound of integral (constraint: a<b)
//
// epsrel: requested relative accuracy [epsrel<0 --> will not show error for failed convergence]
// epsabs: requested absolute accuracy
//
// *fptr : pointer to function f(x, pars) [pars can be set optionally]
// *neval: on exit it will contain the number of function evals
// *r: on exit it will be equal to the best estimate of the integral
//
// *pars: user parameters to be passed on in function call. Inside the function a recast is assumed
//  pmax: explicitly state maximal level of 2^pmax for function evaluations
//        If pmax < 0 --> will not use error checking 
//
// mess: optional message label for debugging
//==================================================================================================
int Integrate_using_Chebyshev(double a, double b, int pmax, double epsrel, double epsabs,
                              double (*fptr)(double, void *pars), void *pars,
                              int *neval, double *r, string mess="");

int Integrate_using_Chebyshev(double a, double b, double epsrel, double epsabs,
                              double (*fptr)(double, void *pars), void *pars,
                              int *neval, double *r, string mess="");

//==================================================================================================
int Integrate_using_Chebyshev(double a, double b, int pmax, double epsrel, double epsabs,
                              double (*fptr)(double, void *pars),
                              int *neval, double *r, string mess="");

int Integrate_using_Chebyshev(double a, double b, double epsrel, double epsabs,
                              double (*fptr)(double, void *pars),
                              int *neval, double *r, string mess="");

//==================================================================================================
int Integrate_using_Chebyshev(double a, double b, int pmax, double epsrel, double epsabs,
                              double (*fptr)(double), int *neval, double *r, string mess="");

int Integrate_using_Chebyshev(double a, double b, double epsrel, double epsabs,
                              double (*fptr)(double), int *neval, double *r, string mess="");

//==================================================================================================
// Driver that loops over the different Chebyshev rules
//
// a : lower bound of integral
// x0: point inside the integration region (a< x0 < (a+b)/2); it allows
//     deforming the integration domain. If set outside this range x0=(3*a+b)/4.
// b : upper bound of integral (constraint: a<b)
//
// epsrel: requested relative accuracy [epsrel<0 --> will not show error for failed convergence]
// epsabs: requested absolute accuracy
//
// *fptr : pointer to function f(x, pars) [pars can be set optionally]
// *neval: on exit it will contain the number of function evals
// *r: on exit it will be equal to the best estimate of the integral
//
// *pars: user parameters to be passed on in function call. Inside the function a recast is assumed
//
// mess: optional message label for debugging
//==================================================================================================
int Integrate_using_Chebyshev(double a, double x0, double b, double epsrel, double epsabs,
                              double (*fptr)(double, void *pars), void *pars,
                              int *neval, double *r, string mess="");

int Integrate_using_Chebyshev(double a, double x0, double b, double epsrel, double epsabs,
                              double (*fptr)(double, void *pars),
                              int *neval, double *r, string mess="");

int Integrate_using_Chebyshev(double a, double x0, double b, double epsrel, double epsabs,
                              double (*fptr)(double), int *neval, double *r, string mess="");

#endif

//==================================================================================================
//==================================================================================================
