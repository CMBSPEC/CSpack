//==================================================================================================
//  Created by Jens Chluba March 2020.
//==================================================================================================

#ifndef CSpack_equilibrium_solutions_h
#define CSpack_equilibrium_solutions_h

#include <string>

using namespace std;

namespace CSpack_equilibrium_solutions {

double Number_integral(double muc);
double DNumber_integral(double muc);

double Energy_integral(double muc);
double DEnergy_integral(double muc);

double Compute_muc_N   (double DTg_Te, double DN_N    =0.0);
double Compute_muc_rho (double DTg_Te, double Drho_rho=0.0);

double Compute_Drho_rho(double DTg_Te  , double DN_N  =0.0);
double Compute_muc_rhoN(double Drho_rho, double DN_N  =0.0);
double Compute_DT_T    (double Drho_rho, double DN_N  =0.0);
double Compute_DT_T    (double Drho_rho, double DN_N, double muc);

double Compute_Drho_rho_inj(double DTg_Te, double muc);
double Compute_DN_N_inj(double DTg_Te, double muc);

}
#endif

//==================================================================================================
//==================================================================================================
