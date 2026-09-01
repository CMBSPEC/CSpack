//==================================================================================================
// Created by Abir Sarkar on 15/11/2019 and modified by JC.
// Last modified: Sept 2026 [JC+Codex]
// These functions are based on Sarkar, Chluba and Lee, MNRAS, 2019
// (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#ifndef CSpack_scattering_matrix_h
#define CSpack_scattering_matrix_h

#include <string>
#include <vector>

#include "ODE_solver_LA.h"
#include "CSpack_kernel_representation.h"

using namespace std;
using namespace ODE_solver_LA;

namespace CSpack_scattering_matrix {

void set_verbosity(int verb=0);

//==================================================================================================
// Routines for scattering matrix setups
//--------------------------------------------------------------------------------------------------
// inputs:
// xarr  : contains frequency grid points x=h nu/kTe = omega/theta
// theta : kTe/mc^2
// Int_wi: Integral weight factors to turn int f(x) dx == sum Int_wi f(xi)
// nK     : defines number of points per kernel wing for Kernel-representation method
//
// outputs: Msc = wj Pij theta
//          KR  = setup Kernel_representation vector on given grid and temperature
//
// epsilon: optional parameter to compress matrix density [eps<1.0e-4 recommended]
// stim   : include stimulated factors from blackbody in moments
//==================================================================================================
// Generic Kernel_representation setup. The thermal kernel can be a photon wrapper, neutrino wrapper
// or any other function matching thermal_kernel_ptr.
void compute_scattering_matrix(const vector<double> &xarr, double theta,
                               const vector<double> &Int_wi, int nK,
                               vector<vector<double> > &Msc,
                               vector<Kernel_representation> &KR,
                               thermal_kernel_ptr K, void *p,
                               double epsilon=1.0e-40, bool stim=0);

// Photon-compatible wrapper using the default exact kernel representation.
void compute_scattering_matrix(const vector<double> &xarr, double theta,
                               const vector<double> &Int_wi, int nK,
                               vector<vector<double> > &Msc,
                               vector<Kernel_representation> &KR,
                               double epsilon=1.0e-40, bool stim=0);

//==================================================================================================
// Routines for scattering matrix setups
//--------------------------------------------------------------------------------------------------
// inputs:
// xarr  : contains frequency grid points x=h nu/kTe = omega/theta
// theta : kTe/mc^2
// Int_wi: Integral weight factors to turn int f(x) dx == sum Int_wi f(xi)
//
// outputs: Msc = wj Pij theta
//
// type   : type of kernel to be used explicitly ['exact', 'exact+SS_C', 'recoil', 'doppler', 'ur', 'SS_K', 'SS_C']
// epsilon: optional parameter to compress matrix density [eps<1.0e-4 recommended]
//==================================================================================================
void compute_scattering_matrix(const vector<double> &xarr, double theta,
                               const vector<double> &Int_wi,
                               vector<vector<double> > &Msc,
                               string type="exact",
                               double epsilon=1.0e-40);

//==================================================================================================
// Routines for scattering matrix setups
//--------------------------------------------------------------------------------------------------
// inputs:
// xarr  : contains frequency grid points x=h nu/kTe = omega/theta
// theta : kTe/mc^2
//
// outputs: Msc = int Pij dxj_bin * theta
//
// type   : type of kernel to be used explicitly ['exact', 'exact+SS_C', 'recoil', 'doppler', 'ur', 'SS_K', 'SS_C']
// epsilon: optional parameter to compress matrix density [eps<1.0e-4 recommended]
// add_stim: optional parameter to add stimulated scattering effect in blackbody radiation field
//==================================================================================================
void compute_scattering_matrix_bin_averaged(const vector<double> &xarr, double theta,
                                            vector<vector<double> > &Msc,
                                            string type,
                                            double epsilon,
                                            bool add_stim=0);

//==================================================================================================
// Routines for scattering matrix setups
//--------------------------------------------------------------------------------------------------
// inputs:
// xarr    : contains frequency grid points x=h nu/kTg = omega/theta_g
// theta   : kTe/mc^2
// theta_g : kTg/mc^2
//
// outputs: Msc = int Pij dxj_bin * theta_g
//
// type   : type of photon kernel to be used explicitly ['exact', 'exact+SS_C', 'recoil', 'doppler', 'ur', 'SS_K', 'SS_C']
// K, p   : generic thermal kernel evaluator and parameters
// epsilon: optional parameter to compress matrix density [eps<1.0e-4 recommended]
// add_stim: optional parameter to add stimulated scattering effect in blackbody radiation field
//==================================================================================================
// Generic bin-averaged setup. Photon/neutrino details enter only through the supplied thermal kernel.
void compute_scattering_matrix_bin_averaged_II(const vector<double> &xarr,
                                               double theta, double theta_g,
                                               vector<vector<double> > &Msc,
                                               vector<Kernel_representation> &KR,
                                               thermal_kernel_ptr K, void *p,
                                               double epsilon);

// Photon wrapper. Use type="exact+SS_C" when the low-energy SS_C replacement is desired.
void compute_scattering_matrix_bin_averaged_II(const vector<double> &xarr,
                                               double theta, double theta_g,
                                               vector<vector<double> > &Msc,
                                               vector<Kernel_representation> &KR,
                                               string type,
                                               double epsilon, bool add_stim=0);

//==================================================================================================
// Routines for neutrino scattering matrix setups using FD-averaged kernels
//--------------------------------------------------------------------------------------------------
// inputs:
// xarr    : contains frequency grid points x=omega/theta_g
// theta   : scatterer temperature in m_e c^2 units
// theta_g : grid temperature scale in m_e c^2 units
// K/type  : neutrino scattering kernel ['nue_e', 'nue_p', 'nue_ep', 'numu_e', 'numu_p',
//           'nutau_e', 'nutau_p']
// mue     : chemical potential for the FD scatterer distribution
// add_FB  : final-state factors in thermal_kernel_FD [0: none, 1: e blocking,
//           2: e and nu blocking, 3: e blocking and photon stimulation]
//==================================================================================================
void compute_scattering_matrix_neutrino_bin_averaged_II(const vector<double> &xarr,
                                                        double theta, double theta_g,
                                                        vector<vector<double> > &Msc,
                                                        vector<Kernel_representation> &KR,
                                                        kernel_ptr K,
                                                        double mue=0.0, int add_FB=0,
                                                        double epsilon=1.0e-40);

void compute_scattering_matrix_neutrino_bin_averaged_II(const vector<double> &xarr,
                                                        double theta, double theta_g,
                                                        vector<vector<double> > &Msc,
                                                        vector<Kernel_representation> &KR,
                                                        string type="nue_ep",
                                                        double mue=0.0, int add_FB=0,
                                                        double epsilon=1.0e-40);

//==================================================================================================
// total cross section
//--------------------------------------------------------------------------------------------------
// inputs:
// xarr  : contains frequency grid points x = hnu/kTe = omega/theta for which Msc is defined
// Msc   : corresponding scattering matrix
//
// outputs: sigarr_i = sum_j Msc_ij * stim
//
// add_stim: optional parameter to add stimulated scattering effect in blackbody radiation field
// Te_Tg   : allow for difference between blackbody and electron temperature
//==================================================================================================
void compute_sigma_tot(const vector<double> &xarr,
                       const vector<vector<double> > &Msc,
                       vector<double> &sigarr,
                       bool add_stim=0,
                       double Te_Tg=1.0);

}

//==================================================================================================
//
// Scattering matrix representation class
//
//==================================================================================================
class Msc_representation
{
private:

    int verbosity_scat_matrix;
    double The, eps_thresh;

    vector<vector<double> > Sigmas;  // relevant moments

    ODE_solver_matrix Msc_sparse;   // scattering matrix in sparse matrix form

public:

    //==============================================================================================
    ~Msc_representation(){};
    Msc_representation(){ Sigmas.clear(); The=0; eps_thresh=1.0e-6; verbosity_scat_matrix=0; }
    
    Msc_representation(const vector<double> &xearr,
                       const vector<double> &Int_wi,
                       double The, double eps_thresh,
                       int maxMom=0, bool stimMom=0);

    //==============================================================================================
    // This version uses parallel setup by filling a full matrix first and then copying
    //==============================================================================================
    void init(const vector<double> &xearr,
              const vector<double> &Int_wi,
              double The, double eps_thresh,
              int maxMom=0, bool stimMom=0);

    //==============================================================================================
    // This version is useful for setting up multiple scattering matrices,
    // where each core can handle one from the calling program
    //==============================================================================================
    void init_serial(const vector<double> &xearr,
                     const vector<double> &Int_wi,
                     double The, double eps_thresh,
                     int maxMom=0, bool stimMom=0);

    //==============================================================================================
    double Msc(int i, int j){ return Msc_sparse.Get_element(i, j); }
    double Get_Theta(){ return The; }
    double Get_eps_thresh(){ return eps_thresh; }

    const ODE_solver_matrix& Get_Msc() { return Msc_sparse; }
    void Get_Msc(vector<vector<double> > &Msc);   // copy into standard full matrix [usually slower]

    void update_moments(const vector<double> &xearr, int maxMom=0, bool stimMom=0, double Te_Tg=0.0);
    double Get_Sigmak(int k, int i);
    double Get_Sigma0(int i){ return Get_Sigmak(0, i); }
    double Get_Sigma1(int i){ return Get_Sigmak(1, i); }
    double Get_Sigma2(int i){ return Get_Sigmak(2, i); }
};

//==================================================================================================
//
// Scattering matrix representation class for multiple values of Theta
//
//==================================================================================================
class Msc_representation_Te : Msc_representation
{
private:

    int npx;
    int logdens_The, npThe;
    double The_min, The_max, The_curr;
    double eps_interpol;                           // The threshold that decides about interpolation

    vector<double> The_arr;
    vector<vector<double> > Msc_full;
    ODE_solver_LA::ODE_solver_matrix Msc_sparse;   // scattering matrix in sparse matrix form

    vector<Msc_representation> Msc_The; // scattering matrix data for various values of The

    double do_interpol(double lgx, const double *lgxa, const double *ya);
    double Msc(int i, int j, unsigned int intj, double The);

public:

    //==============================================================================================
    ~Msc_representation_Te(){};
    Msc_representation_Te(){ npThe=0; Msc_The.resize(0); eps_interpol=1.0e-4; };
    Msc_representation_Te(const vector<double> &xearr,
                          const vector<double> &Int_wi,
                          double The_min, double The_max, int logdens_The,
                          double eps_thresh, double eps_interpol,
                          int maxMom=0, bool stimMom=0);

    void init(const vector<double> &xearr,
              const vector<double> &Int_wi,
              double The_min, double The_max, int logdens_The,
              double eps_thresh, double eps_interpol,
              int maxMom=0, bool stimMom=0);

    int Get_npThe(){ return npThe; }

    double Msc(int i, int j, double The);
    const ODE_solver_LA::ODE_solver_matrix& Get_Msc(double The);
    void Get_Msc(double The, vector<vector<double> > &Msc);   // copy into standard full matrix
};

#endif
//==================================================================================================
//==================================================================================================
