//====================================================================================================================
//
//  CSpack_main.cpp
//
//  Created by Abir Sarkar on 21/01/2020 with modifications from Jens Chluba.
//  This code illustrates some of the computations with CSpack
//
//====================================================================================================================
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "physical_consts.h"
#include "routines.h"

#include "CSpack.h"

using namespace std;
using namespace CSpack_functions;
using namespace CSpack_kernels;
using namespace CSpack_kernels_FD;
using namespace CSpack_Collision_Terms;
using namespace CSpack_opacity;
using namespace CSpack_kernel_moments;
using namespace CSpack_energy_losses;
using namespace CSpack_equilibrium_solutions;

//====================================================================================================================
int main(int narg, char *args[])
{
    create_directory_if_it_does_not_exist("./outputs", 1);

    vector<double> xa;
    init_xarr(1.0e-4, 1.0e+4, xa, 30, 1, 0);
    
    //ofstream ofiles("./outputs/theta_sc_0_1_exact.dat");
    //ofstream ofiles("./outputs/theta_sc_1_exact_minus.dat");
    ofstream ofiles("./outputs/theta_sc_0_1_recoil.dat");
    ofiles.precision(8);

    double tsc_0=0.0, tsc_1=0.0;
    for(int k=0; k<(int)xa.size(); k++)
    {
        cout << " x = " << xa[k] << " " << flush;
        ofiles << xa[k] << " ";
        
        tsc_0=theta_tau_sc_equal_unity(xa[k], 0, CSpack_kernels_nu::kernel_exact_nue_ep, tsc_0);
        tsc_1=theta_tau_sc_equal_unity(xa[k], 1, CSpack_kernels_nu::kernel_exact_nue_ep, tsc_1);

        //tsc_0=theta_tau_sc_equal_unity_ODE(xa[k], 0, CSpack_kernels_nu::kernel_exact_nue_ep, 500.0);
        //tsc_1=theta_tau_sc_equal_unity_ODE(xa[k], 1, CSpack_kernels_nu::kernel_exact_nue_ep, 500.0);
        
        cout << tsc_0 << " " << tsc_1 << endl;
        ofiles << tsc_0 << " " << tsc_1 << endl;
    }
    ofiles.close();
/*
    vector<double> thea;
    init_xarr(0.1, 300.0, thea, 50, 1, 0);
    
    ofstream ofiles("./outputs/x_sc_0_1_exact.dat");
    ofiles.precision(8);

    double xsc_0=0.0, xsc_1=0.0;
    for(int k=0; k<(int)thea.size(); k++)
    {
        cout << " theta = " << thea[k] << " " << flush;
        ofiles << thea[k] << " ";
        
        xsc_0=x_tau_sc_equal_unity(thea[k], 0, CSpack_kernels_nu::kernel_exact_nue_ep, xsc_0);
        xsc_1=x_tau_sc_equal_unity(thea[k], 1, CSpack_kernels_nu::kernel_exact_nue_ep, xsc_1);
        
        cout << xsc_0 << " " << xsc_1 << endl;
        ofiles << xsc_0 << " " << xsc_1 << endl;
    }
    ofiles.close();
*/
    exit(0);
    
    double Th_e=1.0;
    cout << Moment(0.1, 0, Th_e, 0.0, CSpack_kernels_nu::kernel_exact_nue_e) << endl;
    cout << Moment(0.1, 1, Th_e, 0.0, CSpack_kernels_nu::kernel_exact_nue_e) << endl;
    cout << tau_sc_generalized(1.0, 2.0, 0, CSpack_kernels_nu::kernel_exact_nue_ep) << endl;
    exit(0);

    double Th_r=1.0;
    string add="_1.0";
    output_Collision_Term("./outputs/Coll_ph_no_FB"+add+".dat", 0.01, 50.0, 400, Th_r, -1.0e+2, "ph");
    output_Collision_Term("./outputs/Coll_ph_FB"+add+".dat", 0.01, 50.0, 400, Th_r, 0.0, "ph");
    output_Collision_Term("./outputs/Coll_nu_no_FB"+add+".dat", 0.01, 50.0, 400, Th_r, -1.0e+2, "nu");
    output_Collision_Term("./outputs/Coll_nu_FB"+add+".dat", 0.01, 50.0, 400, Th_r, 0.0, "nu");
    exit(0);

    output_kernel("./outputs/kernel_nue_e_sc.dat", 2000, 0.5, 1.0, CSpack_kernels_nu::kernel_exact_nue_e);
    output_kernel("./outputs/kernel_nue_p_sc.dat", 2000, 0.5, 1.0, CSpack_kernels_nu::kernel_exact_nue_p);
    output_kernel("./outputs/kernel_numutau_e_sc.dat", 2000, 0.5, 1.0, CSpack_kernels_nu::kernel_exact_nutau_e);
    output_kernel("./outputs/kernel_numutau_p_sc.dat", 2000, 0.5, 1.0, CSpack_kernels_nu::kernel_exact_nutau_e);

    exit(0);
    
    double omega0=0.1;
    vector<double> p0a={omega0/10.0, omega0/2.0, omega0, omega0*2.0, omega0*5.0};
    
    output_kernel("./outputs/kernel_ph.dat", 2000, omega0, p0a, CSpack_kernels::kernel_exact);
    output_kernel("./outputs/kernel_nu.dat", 2000, omega0, p0a, CSpack_kernels_nu::kernel_exact_nue_e);

    vector<double> o0a={0.001, 0.01, 0.1, 1.0, 10.0};

    output_thermal_kernel("./outputs/kernel_ph_FD.dat", 500, o0a, 0.1, -1.0e+2, CSpack_kernels::kernel_exact);
    output_thermal_kernel("./outputs/kernel_nu_FD.dat", 500, o0a, 0.1, -1.0e+2, CSpack_kernels_nu::kernel_exact_nue_e);

    output_thermal_kernel("./outputs/kernel_ph_FD_BB.dat", 500, o0a, 0.1, -1.0e+2, CSpack_kernels::kernel_exact, 3);
    output_thermal_kernel("./outputs/kernel_nu_FD_FB.dat", 500, o0a, 0.1, -1.0e+2, CSpack_kernels_nu::kernel_exact_nue_e, 2);

    vector<double> o0b={0.01, 0.1, 1.0, 10.0, 100.0};

    output_thermal_kernel("./outputs/kernel_ph_FD_2.0.dat", 500, o0b, 2.0, -1.0e+2, CSpack_kernels::kernel_exact);
    output_thermal_kernel("./outputs/kernel_nu_FD_2.0.dat", 500, o0b, 2.0, -1.0e+2, CSpack_kernels_nu::kernel_exact_nue_e);

    output_thermal_kernel("./outputs/kernel_ph_FD_2.0_FB_mue_0.0.dat", 500, o0b, 2.0, 0.0e+0, CSpack_kernels::kernel_exact, 1);
    output_thermal_kernel("./outputs/kernel_nu_FD_2.0_FB_mue_0.0.dat", 500, o0b, 2.0, 0.0e+0, CSpack_kernels_nu::kernel_exact_nue_e, 1);

    output_thermal_kernel("./outputs/kernel_ph_FD_2.0_FB_mue_0.0_BB.dat", 500, o0b, 2.0, 0.0e+0, CSpack_kernels::kernel_exact, 3);
    output_thermal_kernel("./outputs/kernel_nu_FD_2.0_FB_mue_0.0_FB.dat", 500, o0b, 2.0, 0.0e+0, CSpack_kernels_nu::kernel_exact_nue_e, 2);

    exit(1);

    output_kernel_moments("./outputs/kernel_moment_0_ph.dat", 500, o0a, 0, CSpack_kernels::kernel_exact);
    output_kernel_moments("./outputs/kernel_moment_1_ph.dat", 500, o0a, 1, CSpack_kernels::kernel_exact);
    output_kernel_moments("./outputs/kernel_moment_2_ph.dat", 500, o0a, 2, CSpack_kernels::kernel_exact);

    output_kernel_moments("./outputs/kernel_moment_0_nu.dat", 500, o0a, 0, CSpack_kernels_nu::kernel_exact_nue_e);
    output_kernel_moments("./outputs/kernel_moment_1_nu.dat", 500, o0a, 1, CSpack_kernels_nu::kernel_exact_nue_e);
    output_kernel_moments("./outputs/kernel_moment_2_nu.dat", 500, o0a, 2, CSpack_kernels_nu::kernel_exact_nue_e);

    exit(1);

    Kernel_representation Kth(1.0e-4, 0.01, 0.01, 100, 0.1, 1.0e-10, 1.0e-6, 2);

    cout << Kth.Kernel(0.01*0.9999) << " " << Kth.Kernel(0.01*0.99) << " "
         << Kth.Kernel(0.01*0.95) << " " << Kth.Kernel(0.01*1.1) << " " << Kth.Kernel(0.01) << endl;

    Kernel_representation Kth2(1.0e-6, 1.0, 2.0, 100, 0.1, 1.0e-10, 1.0e-6, 2);
    cout << Kth2.Kernel(1.0*0.9999) << " " << Kth2.Kernel(1.0*0.99) << " "
         << Kth2.Kernel(1.0*0.95) << " " << Kth2.Kernel(1.0*1.1) << endl;

    cout << Compute_muc_N(-0.001) << " " << Compute_muc_rho(-0.001) << " "
         << Compute_Drho_rho(-0.0001) << " " << Compute_muc_N(-0.0001)/1.401 << endl;

//    double omega0=0.000001, p0=20.0;

//    cout << compute_exact_moments_analytical(omega0, p0, 0) << endl;
//    cout << compute_exact_moments_analytical(omega0, p0, 1) << endl;
//    cout << compute_exact_moments_analytical(omega0, p0, 2) << endl;

//    //================================================================================================================
//    int np=500;
//    bool stim=1;
//
//    vector<double> oarr(np);
//    init_xarr(1.0e-3, 1.0e+9, &oarr[0], np, 1, 0); // photon or electron kinetic energy in keV
//
////    ofstream ofile("./outputs/photon_losses.dat");
//    ofstream ofile("./outputs/electron_losses.dat");
//    ofile.precision(8);
//
//    int nz=5;
//    //double za[7]={1.0e+3, 1.0e+4, 5.0e+4, 1.0e+5, 1.0e+6, 5.0e+6, 1.0e+7};
//    double za[5]={1.0e+3, 1.0e+4, 1.0e+5, 1.0e+6, 1.0e+7};
//
//    for(int k=0; k<np; k++)
//    {
//        // MeV units
//        ofile << oarr[k]/1.0e+3 << " ";
//
//        for(int l=0; l<nz; l++)
//        {
//            double theta=2.7255*(1.0+za[l])*const_kb_mec2;
//            //ofile << -photon_energy_loss(oarr[k]/const_me, theta, "exact", stim) << " ";
//            double p0=pfunc(1.0+oarr[k]/const_me);
//            ofile << electron_cooling(theta, p0, "exact", stim) << " ";
//        }
//
//        ofile << endl;
//    }

    //================================================================================================================
    int np=500;
    bool stim=1;
    double Ecr=1.0e+6;

    vector<double> xarr(np);
    init_xarr(1.0e-4, 1.0e+6, &xarr[0], np, 1, 0); // CMB photon energy x=hnu/kTg

    ofstream ofile("./outputs/CS_removal_DC_addition.dat");
    ofile.precision(8);

    int nz=5;
    double za[5]={1.0e+3, 1.0e+4, 1.0e+5, 1.0e+6, 1.0e+7};

    for(int k=0; k<np; k++)
    {
        ofile << xarr[k] << " ";

        double x3=pow(xarr[k], 3);
        for(int l=0; l<nz; l++)
        {
            double thg=2.7255*(1.0+za[l])*const_kb_mec2;
            double p0=pfunc(1.0+Ecr/const_me);
            //double dI_dtau_CS=x3*dng_dtau_removal(xarr[k], thg, p0, "exact", stim);
            double dI_dtau_DC=x3*dng_dtau_DC_add(xarr[k], thg, p0, "exact", stim);
            double dI_dtau_CS=x3*dDng_dtau(xarr[k], thg, p0, "exact", stim);
            ofile << dI_dtau_CS << " " << dI_dtau_DC << " " << dI_dtau_CS+dI_dtau_DC << " ";
        }

        ofile << endl;
    }

//    //================================================================================================================
//    int np=500;
//    bool stim=1;
//
////    p0=0.01;
////    cout << integrand_GDC(omega0, p0) << " " << 1.0+2.0*p0*p0 << endl;
////    cout << Ng_DC_add(0.01, p0, "exact", stim) << " " << Ng_DC_add_approx(0.01, p0, "exact", stim) << endl;
////    //exit(1);
//
//    vector<double> oarr(np);
//    init_xarr(1.0e-3, 1.0e+9, &oarr[0], np, 1, 0); // photon or electron kinetic energy in keV
//
//    ofstream ofile("./outputs/N_CS_removal_DC_addition_approx.dat");
//    ofile.precision(8);
//
//    int nz=5;
//    double za[5]={1.0e+3, 1.0e+4, 1.0e+5, 1.0e+6, 1.0e+7};
//
//    for(int k=0; k<np; k++)
//    {
//        // MeV units
//        ofile << oarr[k]/1.0e+3 << " ";
//
//        for(int l=0; l<nz; l++)
//        {
//            double theta=2.7255*(1.0+za[l])*const_kb_mec2;
//            double p0=pfunc(1.0+oarr[k]/const_me);
//            double Ng_CS=DNg_dtau(theta, p0, "exact", stim);
//            double Ng_DC=Ng_DC_add_approx(theta, p0, "exact", stim);
//
//            ofile << -Ng_CS << " " << Ng_DC << " " << Ng_CS+Ng_DC << " ";
//        }
//
//        ofile << endl;
//    }

//    vector<double> zarr(np);
//    init_xarr(1.0e+3, 1.0e+7, &zarr[0], np, 1, 0);
//
//    create_directory_if_it_does_not_exist("./outputs", 1);
//    ofstream ofile("./outputs/photon_losses.dat");
//    ofile.precision(8);
//
//    int no=7;
//    double oa[7]={1.0e+0, 1.0e+1, 1.0e+2, 1.0e+3, 1.0e+4, 1.0e+5, 1.0e+6};
//
//    for(int k=0; k<np; k++)
//    {
//        double theta=2.7255*(1.0+zarr[k])*const_kb_mec2;
//        ofile << zarr[k] << " " << theta << " "
//
//        for(int l=0; l<no; l++)
//            ofile << photon_energy_loss(omega0[l]/const_me, theta, "exact", stim) << " "
//
//        ofile << endl;
//    }

    ofile.close();
    //================================================================================================================

    return 0;
}

//====================================================================================================================
//====================================================================================================================
