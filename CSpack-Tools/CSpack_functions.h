//
//  CSpack_functions.hpp
//  
//
//  Created by Abir Sarkar on 15/11/2019.
//

#ifndef CSpack_functions_hpp
#define CSpack_functions_hpp

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <limits.h>
#include <vector>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_dilog.h>
#include <limits>
#include <complex>
//======================================================================================
// General Functions
//======================================================================================
//======================================================================================
// Critical Frequencies
//======================================================================================
double omegacrit(double omega0, double p0);
double omegamin(double omega0, double p0);
double omegamax(double omega0, double p0);
//======================================================================================
// Common functions
//======================================================================================
double gamma_ei(double p0);
double p_ef(double g);
double beta(double p0);
double gamma_ef(double omega0, double p0, double omega);
double pfunc_low(double omega0, double p0, double omega);


//======================================================================================
// Kernel
//======================================================================================
double lambda_p(double p0, double omega0);
double lambda_m(double p0, double omega);
double sfunc(double x);
double ffunc(double x);
double sfunc_con(double x);
double ffunc_con(double x);
double omegabar(double omega0, double p0, double omega);
double omega0bar(double omega0, double p0, double omega);
double kappa(double omega0,  double omega, double p0, double p);

double mb_dist_func(double p, double theta);
double mb_dist_norm(double the);
double rel_mb_dist(double p, double theta);

//======================================================================================
// Moment
//======================================================================================
double alphap(double omega0, double p0);
double alpham(double omega0, double p0);
double f_moment(double omega0, double p0);


//

#endif /* CSpack_functions_hpp */
