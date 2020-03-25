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
using namespace CSpack_kernel_moments;
using namespace CSpack_energy_losses;

//====================================================================================================================
int main(int narg, char *args[])
{
    Kernel_representation Kth(1.0e-4, 0.01, 0.01, 100, 0.1, 1.0e-10, 1.0e-6, 2);

    cout << Kth.Kernel(0.01*0.9999) << " " << Kth.Kernel(0.01*0.99) << " "
         << Kth.Kernel(0.01*0.95) << " " << Kth.Kernel(0.01*1.1) << " " << Kth.Kernel(0.01) << endl;

    Kernel_representation Kth2(1.0e-6, 1.0, 2.0, 100, 0.1, 1.0e-10, 1.0e-6, 2);
    cout << Kth2.Kernel(1.0*0.9999) << " " << Kth2.Kernel(1.0*0.99) << " "
         << Kth2.Kernel(1.0*0.95) << " " << Kth2.Kernel(1.0*1.1) << endl;

    double omega0=0.000001, p0=20.0;

//    cout << compute_exact_moments_analytical(omega0, p0, 0) << endl;
//    cout << compute_exact_moments_analytical(omega0, p0, 1) << endl;
//    cout << compute_exact_moments_analytical(omega0, p0, 2) << endl;

    create_directory_if_it_does_not_exist("./outputs", 1);

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
