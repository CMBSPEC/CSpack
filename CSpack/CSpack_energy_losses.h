//==================================================================================================
//  CSpack_energy_losses.h
//
//  Created by Jens Chluba March 2020.
//==================================================================================================

#ifndef CSpack_energy_losses_h
#define CSpack_energy_losses_h

#include <string>

using namespace std;

namespace CSpack_energy_losses {

//==================================================================================================
//  losses of photon scattering off thermal electrons with ambient blackbody radiation (Tg=Te)
//==================================================================================================
// d(ln rho_h)/dtau
//==================================================================================================
double photon_energy_loss(double omega0, double the, string type, bool stim=0);

//==================================================================================================
//  losses of electron scattering off ambient blackbody radiation at temperature Tg
//==================================================================================================
double electron_cooling(double thg, double p0, string type, bool stim=0);
double photon_gains(double thg, double p0, string type, bool stim=0);
double photon_gains_approx(double thg, double p0);

double integrand_GDC(double omega0, double p0);
double integrand_GDC_npl(double Thg, double p0);

double dDng_dtau(double x, double thg, double p0, string type, bool stim);
double dng_dtau_removal(double x, double thg, double p0, string type, bool stim);

double DNg_dtau(double thg, double p0, string type, bool stim);
double Ng_dtau_removal(double thg, double p0, string type, bool stim);

double dng_dtau_DC_add(double x, double thg, double p0, string type, bool stim);
double Ng_DC_add(double thg, double p0, string type, bool stim);
double Ng_DC_add_approx(double thg, double p0, string type, bool stim);

}

#endif

//==================================================================================================
//==================================================================================================
