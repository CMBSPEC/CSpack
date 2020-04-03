//==================================================================================================
//  Created by Jens Chluba March 2020.
//==================================================================================================

#ifndef CSpack_equilibrium_solutions_h
#define CSpack_equilibrium_solutions_h

#include <string>

using namespace std;

namespace CSpack_equilibrium_solutions {

double Compute_muc_N   (double DTg_Te, double DN_N    =0.0);
double Compute_muc_rho (double DTg_Te, double Drho_rho=0.0);
double Compute_Drho_rho(double DTg_Te, double DN_N    =0.0);

}

#endif

//==================================================================================================
//==================================================================================================
