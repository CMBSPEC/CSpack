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

}

#endif

//==================================================================================================
//==================================================================================================
