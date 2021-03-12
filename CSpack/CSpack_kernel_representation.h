//==================================================================================================
//  Created by JC in March 2020.
//==================================================================================================

#ifndef CSpack_kernel_representation_h
#define CSpack_kernel_representation_h

#include <string>
#include <vector>

using namespace std;

//==================================================================================================
class Kernel_representation
{
private:

    double omega0, Theta;
    double P0, wmin, wmax;
    int spline_up, spline_down, np;

    vector<double> Moments;
    double G, H;

    void create_kernel_splines(double omega_lim, double eps_thresh, int np, string type);
    double compute_moment(int k, bool stim);
    double compute_G(bool stim);
    double compute_H(bool stim);

    
public:

    //==============================================================================================
    ~Kernel_representation();
    Kernel_representation();
    Kernel_representation(double omin, double om0, double omax, int np,
                          double The,
                          double eps_thresh, double eps_interpol,
                          int maxMom=0, bool stim=0);

    //==================================================================================================
    // for openmp runs this should be ran serial before the init call...
    //==================================================================================================
    void allocate_splines(int np);

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

#endif
//==================================================================================================
//==================================================================================================
