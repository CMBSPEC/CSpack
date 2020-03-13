//==================================================================================================
// CSpack.h
//
//
// Created by Abir Sarkar on 12/11/2019.
// Modifications by Jens Chluba, Feb 2020
//==================================================================================================

#ifndef cspack_h
#define cspack_h

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
#include <CSpack_functions.h>
#include "physical_consts.h"

using namespace std;
//==================================================================================================
// KERNEL
//==================================================================================================
double kernel_exact(double omega0, double p0, double omega);
double kernel_recoil(double omega0, double p0, double omega);
double kernel_doppler(double omega0, double p0, double omega);
double kernel_ur(double omega0, double p0, double omega);

//==================================================================================================
//  MOMENT
//==================================================================================================
double compute_exact_moments_analytical(double omega0, double p0, int l);
double compute_nr_moments_analytical(double omega0, double p0, int l);
double compute_ur_moments_analytical(double omega0, double p0, int l);
double kernel_recoil_moment(double omega0, double p0, int l);
double kernel_dop_moment(double omega0, double p0, int l);

//==================================================================================================
//  THERMALLY AVERAGED KERNEL
//==================================================================================================
double thermal_kernel_exact(double omega0, double omega, double theta);
double thermal_kernel_recoil(double omega0, double omega, double theta);
double thermal_kernel_ur(double omega0, double omega, double theta);
double thermal_kernel_doppler(double omega0, double omega, double theta);
double thermal_kernel_all(double omega0, double omega, double theta, int type);

//double thermal_kernel_intpol(double omega0, double omega, double theta);
double thermal_kernel_intpol(vector<double> &om0arr, vector<double> &omarr, vector<double> &the0arr, vector<vector<vector<double>>> &int_mat, double omega0, double omega, double theta);
double thermal_kernel_intpol_spline(vector<double> &om0arr, vector<double> &omarr, vector<double> &the0arr, vector<vector<vector<double>>> &int_mat, double omega0, double omega, double theta);

//==========================================================================================
//INFO FOR DATABASE
//==========================================================================================
struct ker_info {
    const int nomom0 = 100, nomthe = 100;
    double om0min = 0.01, om0max = 150.1;
    double the0min = 5.0/510.9989461, the0max = 160.0/510.9989461;
    double Te = 1.1509e+07, Tg = Te;
    string ifname = "./data_for_interpolation/data_for_interpolation_4.dat";
};

//==================================================================================================
// THERMALLY AVERAGED MOMENT
//==================================================================================================
double moment_Int_therm(double omega0, int l, double theta);
double moment_Int_app1(double omega0, int l, double theta);
double moment_Int_app2(double omega0, int l, double theta);
double moment_Int_app3(double omega0, int l, double theta);

//==================================================================================================
// THERMALLY AVERAGED MOMENTS using 2D integration
//==================================================================================================
double moment_2D_Int_therm_all(double omega0, int k, double theta, int type); // JC: not as precise
double moment_2D_Int_therm_all_II(double omega0, int k, double theta, int type);

// 2Sigma2-Sigma1
double G_moment_2D_Int_therm_all_II(double omega0, double theta, int type);

//==================================================================================================
// obtain coefficients for FP approximation dn/dy= A(x)*d^2n/dx^2 + B(x)*dn/dx + C(x)*n
//==================================================================================================
void compute_all_FP_coefficients(double omega0, double theta, int type,
                                 double &A, double &B, double &C);

#endif /* cspack_h */
