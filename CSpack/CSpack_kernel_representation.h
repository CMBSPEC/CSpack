//==================================================================================================
// Created by JC in March 2020.
// Last modified: Sept 2026 [JC+Codex]
//==================================================================================================

#ifndef CSpack_kernel_representation_h
#define CSpack_kernel_representation_h

#include <string>
#include <vector>

#include "CSpack_kernels.h"

using namespace std;

//==================================================================================================
// Generic thermal kernel interface used by Kernel_representation. The last argument can point to
// kernel-specific parameters, which lets the same spline machinery handle photons, FD-averaged
// neutrinos, or another thermal kernel without changing the representation class.
//==================================================================================================
typedef double (*thermal_kernel_ptr)(double, double, double, void *);

//==================================================================================================
// Parameter block for the FD-averaged neutrino wrapper below. The kernel itself stores the
// redistribution shape; sigma_norm*x0^2 applies the neutrino cross-section normalization relative to
// Thomson scattering, with x0=omega0/theta_g.
//==================================================================================================
struct Kernel_representation_nu_params
{
    double mue, theta_g, sigma_norm;
    kernel_ptr K;
    int add_FB;

    Kernel_representation_nu_params(){ mue=0.0; theta_g=0.0; sigma_norm=1.0; K=NULL; add_FB=0; }
};

double thermal_kernel_photon_KR(double omega0, double omega, double theta, void *p);
double thermal_kernel_neutrino_KR(double omega0, double omega, double theta, void *p);

//==================================================================================================
class Kernel_representation
{
private:

    double omega0, Theta;
    double P0, wmin, wmax;
    int spline_up, spline_down, np;

    vector<double> Moments;
    double G, H;

    void create_kernel_splines(double omega_lim, double eps_thresh, int np,
                               thermal_kernel_ptr K, void *p);
    double compute_moment(int k, bool stim);
    double compute_G(bool stim);
    double compute_H(bool stim);

    
public:

    //==============================================================================================
    ~Kernel_representation();
    Kernel_representation();

    // Generic constructor/init use a caller-supplied thermal kernel. The following photon-only
    // overloads are kept for backward compatibility and call the generic version internally.
    Kernel_representation(double omin, double om0, double omax, int np,
                          double The,
                          double eps_thresh, double eps_interpol,
                          thermal_kernel_ptr K, void *p,
                          int maxMom=0, bool stim=0);
    Kernel_representation(double omin, double om0, double omax, int np,
                          double The,
                          double eps_thresh, double eps_interpol,
                          int maxMom=0, bool stim=0);

    //==============================================================================================
    // for openmp runs this function should be ran serial before the init call...
    //==============================================================================================
    void allocate_splines(int np);

    void init(double omin, double om0, double omax, int np, double The,
              double eps_thresh, double eps_interpol,
              thermal_kernel_ptr K, void *p,
              int maxMom=0, bool stim=0);
    void init(double omin, double om0, double omax, int np, double The,
              double eps_thresh, double eps_interpol, int maxMom=0, bool stim=0);

    //==============================================================================================
    double Kernel(double om);
    double Get_Theta(){ return Theta; }
    double Get_omega0(){ return omega0; }
    double Get_P0(){ return P0; }
    double Get_omega_min() { return omega0*wmin; }
    double Get_omega_max() { return omega0*wmax; }
    double Get_Mom(int k){ return (k<(int)Moments.size() ? Moments[k] : 0.0); }
    double Get_G(){ return G; }
    double Get_H(){ return H; }
};

//==================================================================================================
//
// Kernel representation class for multiple values of Theta
//
//==================================================================================================
class Kernel_representation_Te : Kernel_representation
{
private:

    int logdens_The, npThe;
    double Theta_min, Theta_max;

    vector<double> The_arr;
    vector<Kernel_representation> Kernels_The; // Kernel data for various values of The

    void init_Kernel_Table(double omin, double om0, double omax, int npom,
                           double Theta_min, double Theta_max, int logdens_The,
                           double eps_thresh, double eps_interpol,
                           int maxMom, bool stim);

    double lgThe[4], lgK[4], D[4], a[4];
    double do_interpol(double lgx, const double *lgxa, const double *ya);

public:

    //==============================================================================================
    ~Kernel_representation_Te(){};
    Kernel_representation_Te(){ npThe=0; Kernels_The.resize(0); };
    Kernel_representation_Te(double omin, double om0, double omax, int npom,
                             double Theta_min, double Theta_max, int logdens_The,
                             double eps_thresh, double eps_interpol,
                             int maxMom=0, bool stim=0);

    double Kernel(double om, double The);
};

#endif
//==================================================================================================
//==================================================================================================
