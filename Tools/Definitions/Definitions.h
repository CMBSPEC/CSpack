//-------------------------------------------------------------------------------------------------------
// Author Jens Chluba Oct 2010
//-------------------------------------------------------------------------------------------------------
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <cstdlib>

#ifdef COSMORECPATH
#define COSMORECDIR ((string)COSMORECPATH)
#else
#define COSMORECDIR ((string)"")
#endif

//-------------------------------------------------------------------------------------------------------
#define SQRT_PI 1.7724538509055160273 
#define PI 3.1415926535897932385
#define PIS 2.5066282746310005024
#define PI2 9.8696044010893586188
#define PI3 31.006276680299820175
#define PI4 97.409091034002437236
#define TWOPI 6.2831853071795864769
#define FOURPI 12.566370614359172954
#define EULERGAMMA 0.57721566490153286061
//-------------------------------------------------------------------------------------------------------
#define SQRT_2 1.4142135623730950488

const double G11_Int_pl=1.644934066848226;
const double G21_Int_pl=2.404113806319189;
const double G31_Int_pl=6.493939402266829;
const double G41_Int_pl=24.88626612344088;
const double G51_Int_pl=122.0811674381339;

const double alpha_rho_pl=G21_Int_pl/G31_Int_pl;
const double alpha_mu=0.4561442592067353;
const double I4_Int_pl =25.975757609067315;
const double kappa_Int_pl=2.141851504502378;

#endif
