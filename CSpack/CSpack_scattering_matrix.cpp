//==================================================================================================
//  Created by Abir Sarkar on 15/11/2019 and modified by JC. These functions are based on
//  Sarkar, Chluba and Lee, MNRAS, 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#include "routines.h"
#include "Patterson.h"

#include "CSpack.h"

using namespace std;
using namespace CSpack_functions;
using namespace CSpack_kernels;
using namespace CSpack_weights;

//==================================================================================================
// main functions in namespace
//==================================================================================================
namespace CSpack_scattering_matrix {

int verbosity_scat_matrix=0;

void set_verbosity(int verb){ verbosity_scat_matrix=verb; return; }

//==================================================================================================
// Routines for scattering matrix setups
//--------------------------------------------------------------------------------------------------
// inputs :
// xarr   : contains frequency grid points x=h nu/kTe = omega/theta
// theta  : kTe/mc^2
// Int_wi : Integral weight factors to turn int f(x) dx == sum Int_wi f(xi)
// nK     : defines number of points per kernel wing for Kernel-representation method
//
// outputs: Msc = wj Pij theta
//          KR  = setup Kernel_representation vector on given grid and temperature
//
// epsilon: optional parameter to compress matrix density [eps<1.0e-4 recommended]
// stim   : include stimulated factors from blackbody in moments
//==================================================================================================
void compute_scattering_matrix(const vector<double> &xarr, double theta,
                               const vector<double> &Int_wi, int nK,
                               vector<vector<double> > &Msc,
                               vector<Kernel_representation> &KR,
                               double epsilon, bool stim)
{
    if(verbosity_scat_matrix>0)
        cout << " compute_scattering_matrix :: setting up scattering matrix for The= " << theta << endl;

    int npx=xarr.size();

    // create matrix
    if((int)Msc.size()!=npx)
    {
        Msc.clear();
        vector<double> zeros(npx, 0.0);
        for(int i=0; i<npx; i++) Msc.push_back(zeros);
    }

    // make vector of Kernel representations
    KR.resize(npx);
    double omin=xarr[0]*theta/3.0, omax=xarr.back()*theta*3.0;
    for(int i=0; i<npx; i++) KR[i].allocate_splines(nK);

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int i=0; i<npx; i++)
        KR[i].init(omin, xarr[i]*theta, omax, nK, theta, epsilon, epsilon, -1, stim);
#ifdef OPENMP_ACTIVATED
#pragma omp barrier
#endif

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int i=0; i<npx; i++) //x
    {
        // reset matrix
        for(int j=0; j<npx; j++) Msc[i][j]=0.0;

        double omp=xarr[i]*theta;
        // diagonal element for reference
        Msc[i][i]= Int_wi[i] * theta * KR[i].Kernel(omp);

        for(int j=i+1; j<npx; j++) //xp>x
        {
            omp=xarr[j]*theta;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc[i][j]= Int_wi[j] * theta * KR[i].Kernel(omp);

            if(abs(Msc[i][j]/Msc[i][i])<epsilon) break;
        }

        for(int j=i-1; j>=0; j--) //xp<x
        {
            omp=xarr[j]*theta;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc[i][j]= Int_wi[j] * theta * KR[i].Kernel(omp);

            if(abs(Msc[i][j]/Msc[i][i])<epsilon) break;
        }
    }
#ifdef OPENMP_ACTIVATED
#pragma omp barrier
#endif

    if(verbosity_scat_matrix>0)
        cout << " compute_scattering_matrix :: done." << endl << endl;
}

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
// type   : type of kernel to be used explicitly ['exact', 'SS_K', 'SS_C']
// epsilon: optional parameter to compress matrix density [eps<1.0e-4 recommended]
//==================================================================================================
void compute_scattering_matrix(const vector<double> &xarr, double theta,
                               const vector<double> &Int_wi,
                               vector<vector<double> > &Msc,
                               string type,
                               double epsilon)
{
    if(verbosity_scat_matrix>0)
        cout << " compute_scattering_matrix :: setting up scattering matrix for The= " << theta << endl;

    int npx=xarr.size();

    double (*kernel)(double omega0, double omega, double theta)=NULL;

    if(type=="exact") kernel=thermal_kernel_exact;
    else if(type=="SS_K") kernel=thermal_kernel_SS_K;
    else if(type=="SS_C") kernel=thermal_kernel_SS_C;
    else throw_error("compute_scattering_matrix", "kernel type not available", 1);

    // create matrix
    if((int)Msc.size()!=npx)
    {
        Msc.clear();
        vector<double> zeros(npx, 0.0);
        for(int i=0; i<npx; i++) Msc.push_back(zeros);
    }

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int i=0; i<npx; i++) //x
    {
        // reset matrix
        for(int j=0; j<npx; j++) Msc[i][j]=0.0;

        double om0=xarr[i]*theta, omp=xarr[i]*theta;

        // avoid errors at low energies using exact expressions
        if(type=="exact")
        {
            if(om0<1.0e-4 && theta<1.0e-4) kernel=thermal_kernel_SS_C;
            else kernel=thermal_kernel_exact;
        }

        // diagonal element for reference
        Msc[i][i]= Int_wi[i] * theta * kernel(om0, omp, theta);

        for(int j=i+1; j<npx; j++) //xp>x
        {
            omp=xarr[j]*theta;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc[i][j]= Int_wi[j] * theta * kernel(om0, omp, theta);

            if(abs(Msc[i][j]/Msc[i][i])<epsilon && j-(i+1)>=4) break;
        }

        for(int j=i-1; j>=0; j--) //xp<x
        {
            omp=xarr[j]*theta;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc[i][j]= Int_wi[j] * theta * kernel(om0, omp, theta);

            if(abs(Msc[i][j]/Msc[i][i])<epsilon && i-1-j>=4) break;
        }
    }
#ifdef OPENMP_ACTIVATED
#pragma omp barrier
#endif

    if(verbosity_scat_matrix>0)
        cout << " compute_scattering_matrix :: done." << endl << endl;
}

//==================================================================================================
// routines to do bin averaging
//==================================================================================================
struct Integral_Mij_bin_average
{
    Kernel_representation *KR;
    bool add_stim;
    double x0, om0, theta, theta_g;

    Integral_Mij_bin_average()
    {
        KR=NULL;
        add_stim=0;
    }
};

double Kernel_weight_func(double x, void *p) // x == hnu/kTg == om/theta_g
{
    Integral_Mij_bin_average &IM=*(Integral_Mij_bin_average *)p;

    if(x==-1.0e+300) return IM.KR->Get_omega_min()/IM.theta_g;
    if(x== 1.0e+300) return IM.KR->Get_omega_max()/IM.theta_g;

    double om=x*IM.theta_g;
    double stim=(IM.add_stim ? one_minus_exp_mx(IM.x0)/one_minus_exp_mx(x) : 1.0);

    double K=0.0;
    if(IM.om0<1.0e-4 && IM.theta<1.0e-4) K=thermal_kernel_SS_C(IM.om0, om, IM.theta);
    else K=IM.KR->Kernel(om);
    //else K=thermal_kernel_exact(IM.om0, om, IM.theta);
    //K=IM.KR->Kernel(om);

    return K * stim * IM.theta_g; // domp = theta_g dx --> factor of theta_g
}

void fill_Msc_bin_averaged_II(int i, int j, const vector<double> &xarr,
                              double theta, double theta_g,
                              vector<vector<double> > &Msc,
                              vector<Kernel_representation> &KR,
                              bool add_stim=0)
{
    Integral_Mij_bin_average IM;
    IM.KR=&KR[i];
    IM.add_stim=add_stim;
    IM.x0 =xarr[i];
    IM.om0=IM.x0*theta_g;
    IM.theta  =theta;
    IM.theta_g=theta_g;

    //--------------------------------------------------------------------------
    // order >=2 is recommended to get good energy conservation
    // comment JC: 2-5 all seem to be giving similar results really...
    //--------------------------------------------------------------------------
    Lagrange_Polynomial_weights(j, 3, xarr, Msc[i], Kernel_weight_func, &IM, 0);

    return;
}

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
// type   : type of kernel to be used explicitly ['exact', 'SS_K', 'SS_C']
// epsilon: optional parameter to compress matrix density [eps<1.0e-4 recommended]
// add_stim: optional parameter to add stimulated scattering effect in blackbody radiation field
//==================================================================================================
void compute_scattering_matrix_bin_averaged_II(const vector<double> &xarr,
                                               double theta, double theta_g,
                                               vector<vector<double> > &Msc,
                                               vector<Kernel_representation> &KR,
                                               string type,
                                               double epsilon, bool add_stim)
{
    string funcname="compute_scattering_matrix_bin_averaged";
    if(verbosity_scat_matrix>0)
        cout << " " + funcname + " :: setting up scattering matrix for The= " << theta << endl;

    int npx=xarr.size();

    // create matrix
    if((int)Msc.size()!=npx)
    {
        Msc.clear();
        vector<double> zeros(npx, 0.0);
        for(int i=0; i<npx; i++) Msc.push_back(zeros);
    }

    if((int)KR.size()!=npx) KR.resize(npx);

    //--------------------------------------------------------------------------
    // make vector of Kernel representations
    //--------------------------------------------------------------------------
    int nK=80; // JC: for higher precision, this parameter should be increased
    double omin=xarr[0]*theta_g/3.0, omax=xarr.back()*theta_g*3.0;
    for(int i=0; i<npx; i++) KR[i].allocate_splines(nK);

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int i=0; i<npx; i++)
        KR[i].init(omin, xarr[i]*theta_g, omax, nK, theta, epsilon/2.0, epsilon/2.0, -1, 0);
#ifdef OPENMP_ACTIVATED
#pragma omp barrier
#endif

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int i=0; i<npx; i++) //x
    {
        // reset matrix
        for(int j=0; j<npx; j++) Msc[i][j]=0.0;

        // diagonal element for reference
        fill_Msc_bin_averaged_II(i, i, xarr, theta, theta_g, Msc, KR, add_stim);

        for(int j=i+1; j<npx; j++) //xp>x
        {
            fill_Msc_bin_averaged_II(i, j, xarr, theta, theta_g, Msc, KR, add_stim);

            if(abs(Msc[i][j]/Msc[i][i])<epsilon) break;
        }

        for(int j=i-1; j>=0; j--) //xp<x
        {
            fill_Msc_bin_averaged_II(i, j, xarr, theta, theta_g, Msc, KR, add_stim);

            if(abs(Msc[i][j]/Msc[i][i])<epsilon) break;
        }
    }
#ifdef OPENMP_ACTIVATED
#pragma omp barrier
#endif

    if(verbosity_scat_matrix>0)
        cout << " " + funcname + " :: done." << endl << endl;
}

//==================================================================================================
//
// old versions
//
//==================================================================================================
struct Integral_Info_Mij
{
    int i, k; // k== 0, 1
    const vector<double> *xa;
    double theta, om0, theta_g;
    double (*kernel)(double omega0, double omega, double theta);
    Kernel_representation *KR;
    bool add_stim;

    Integral_Info_Mij()
    {
        xa=NULL;
        kernel=NULL;
        KR=NULL;
        add_stim=0;
    }
};

void Get_oml_omu(int i, const vector<double> &xa, double theta, double &oml, double &omu)
{
    if(i==0)
    {
        oml=xa[0]*theta;
        omu=xa[1]*theta;
    }
    else if(i==(int)xa.size()-1)
    {
        oml=xa[(int)xa.size()-2]*theta;
        omu=xa[(int)xa.size()-1]*theta;
    }
    else
    {
        oml=(xa[i]+xa[i-1])/2.0*theta;
        omu=(xa[i]+xa[i+1])/2.0*theta;
    }

    return;
}

double Weight_Polynomial(int i, int k, double x, const vector<double> &xa)
{
    //return 1.0;
    if(k==0) return 1.0;
    if(k==1) return (x-xa[i])/(xa[i+1]-xa[i-1]);
    return 0.0;
}

//==================================================================================================
double dMij_bin_averaged(double om, void *p)
{
    Integral_Info_Mij &d=*(Integral_Info_Mij *)p;

    double stim=(d.add_stim ? one_minus_exp_mx(d.om0/d.theta_g)/one_minus_exp_mx(om/d.theta_g) : 1.0);

    return d.kernel(d.om0, om, d.theta) * stim * Weight_Polynomial(d.i, d.k, om/d.theta_g, *d.xa);
}

double Mij_bin_averaged(int i, int k, const vector<double> &xa,
                        double om0, double oml, double omp, double omu, double theta,
                        double (*kernel)(double omega0, double omega, double theta),
                        bool add_stim=0)
{
    Integral_Info_Mij d;
    d.i=i;
    d.k=k;
    d.xa=&xa;
    d.theta=theta;
    d.om0=om0;
    d.theta_g=theta;
    d.kernel=kernel;
    d.add_stim=add_stim;

    //    // use analytic formula for non-relativistic limit
    //    if(d.om0<1.0e-4 && d.theta<1.0e-4) d.kernel=thermal_kernel_SS_C;
    //    else d.kernel=thermal_kernel_exact;
    //    d.kernel=thermal_kernel_SS_C;

    double a=oml, b=omu;
    double epsrel=1.0e-6, epsabs=1.0e-50;
    double r=0.0;

    if(a<om0 && om0<b) // split integral across cusp...
    {
        r=Integrate_using_Patterson_adaptive(a, om0, epsrel, epsabs, dMij_bin_averaged, &d);
        r+=Integrate_using_Patterson_adaptive(om0, b, epsrel, epsabs, dMij_bin_averaged, &d);
    }
    else r=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dMij_bin_averaged, &d);

    return r;
}

//==================================================================================================
double dMij_bin_averaged_KR(double om, void *p)
{
    Integral_Info_Mij &d=*(Integral_Info_Mij *)p;

    double stim=(d.add_stim ? one_minus_exp_mx(d.om0/d.theta_g)/one_minus_exp_mx(om/d.theta_g) : 1.0);

    return d.KR->Kernel(om) * stim * Weight_Polynomial(d.i, d.k, om/d.theta_g, *d.xa);
}

double Mij_bin_averaged(int i, int k, const vector<double> &xa,
                        double om0, double oml, double omp, double omu, double theta,
                        Kernel_representation &KR,
                        bool add_stim=0)
{
    Integral_Info_Mij d;
    d.i=i;
    d.k=k;
    d.xa=&xa;
    d.theta=theta;
    d.om0=om0;
    d.theta_g=theta;
    d.KR=&KR;
    d.add_stim=add_stim;

    double a=max(oml, KR.Get_omega_min()), b=min(omu, KR.Get_omega_max());
    double epsrel=1.0e-6, epsabs=1.0e-50;
    double r=0.0;

    if(a>=b) return 0.0;

    if(a<om0 && om0<b) // split integral across cusp...
    {
        r=Integrate_using_Patterson_adaptive(a, om0, epsrel, epsabs, dMij_bin_averaged_KR, &d);
        r+=Integrate_using_Patterson_adaptive(om0, b, epsrel, epsabs, dMij_bin_averaged_KR, &d);
    }
    else r=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dMij_bin_averaged_KR, &d);

    //    if(isinf(r) || isnan(r))
    //    {
    //        cout << i << " " << k << " " << theta << endl;
    //        throw_error("matrix element bad", "", 1);
    //    }

    return r;
}

void fill_Msc_bin_averaged(int i, int j, const vector<double> &xarr, double theta,
                           vector<vector<double> > &Msc,
                           double (*kernel)(double omega0, double omega, double theta),
                           bool add_stim=0)
{
    double om0=xarr[i]*theta, oml, omp=xarr[j]*theta, omu;

    // get limits of frequency bin
    Get_oml_omu(j, xarr, theta, oml, omu);

    Msc[i][j]+=Mij_bin_averaged(j, 0, xarr, om0, oml, omp, omu, theta, kernel, add_stim);

//    if(j>0 && j<(int)xarr.size()-1)
//    {
//        double Mtemp=Mij_bin_averaged(j, 1, xarr, om0, oml, omp, omu, theta, kernel, add_stim);
//        Msc[i][j+1]+=Mtemp;
//        Msc[i][j-1]-=Mtemp;
//    }

    return;
}

void fill_Msc_bin_averaged(int i, int j, const vector<double> &xarr, double theta,
                           vector<vector<double> > &Msc,
                           vector<Kernel_representation> &KR,
                           bool add_stim=0)
{
    double om0=xarr[i]*theta, oml, omp=xarr[j]*theta, omu;

    // get limits of frequency bin
    Get_oml_omu(j, xarr, theta, oml, omu);

    Msc[i][j]+=Mij_bin_averaged(j, 0, xarr, om0, oml, omp, omu, theta, KR[i], add_stim);

    if(j>0 && j<(int)xarr.size()-1)
    {
        double Mtemp=Mij_bin_averaged(j, 1, xarr, om0, oml, omp, omu, theta, KR[i], add_stim);
        Msc[i][j+1]+=Mtemp;
        Msc[i][j-1]-=Mtemp;
    }

    return;
}

//==================================================================================================
// Routines for scattering matrix setups
//--------------------------------------------------------------------------------------------------
// inputs:
// xarr  : contains frequency grid points x=h nu/kTe = omega/theta
// theta : kTe/mc^2
//
// outputs: Msc = int Pij dxj_bin * theta
//
// type   : type of kernel to be used explicitly ['exact', 'SS_K', 'SS_C']
// epsilon: optional parameter to compress matrix density [eps<1.0e-4 recommended]
// add_stim: optional parameter to add stimulated scattering effect in blackbody radiation field
//==================================================================================================
void compute_scattering_matrix_bin_averaged(const vector<double> &xarr, double theta,
                                            vector<vector<double> > &Msc,
                                            string type,
                                            double epsilon,
                                            bool add_stim)
{
    string funcname="compute_scattering_matrix_bin_averaged";
    if(verbosity_scat_matrix>0)
        cout << " " + funcname + " :: setting up scattering matrix for The= " << theta << endl;

    int npx=xarr.size();

//    double (*kernel)(double omega0, double omega, double theta)=NULL;
//    if(type=="exact") kernel=thermal_kernel_exact;
//    else if(type=="SS_K") kernel=thermal_kernel_SS_K;
//    else if(type=="SS_C") kernel=thermal_kernel_SS_C;
//    else throw_error(funcname, "kernel type not available", 1);

    // create matrix
    if((int)Msc.size()!=npx)
    {
        Msc.clear();
        vector<double> zeros(npx, 0.0);
        for(int i=0; i<npx; i++) Msc.push_back(zeros);
    }

    // make vector of Kernel representations
    int nK=30;
    vector<Kernel_representation> KR(npx);
    double omin=xarr[0]*theta/3.0, omax=xarr.back()*theta*3.0;
    for(int i=0; i<npx; i++) KR[i].allocate_splines(nK);

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int i=0; i<npx; i++)
        KR[i].init(omin, xarr[i]*theta, omax, nK, theta, epsilon, epsilon, -1, 0);
#ifdef OPENMP_ACTIVATED
#pragma omp barrier
#endif

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int i=0; i<npx; i++) //x
    {
        // reset matrix
        for(int j=0; j<npx; j++) Msc[i][j]=0.0;

        // diagonal element for reference
        //fill_Msc_bin_averaged(i, i, xarr, theta, Msc, kernel, add_stim);
        fill_Msc_bin_averaged(i, i, xarr, theta, Msc, KR, add_stim);

        for(int j=i+1; j<npx; j++) //xp>x
        {
            //fill_Msc_bin_averaged(i, j, xarr, theta, Msc, kernel, add_stim);
            fill_Msc_bin_averaged(i, j, xarr, theta, Msc, KR, add_stim);

            if(abs(Msc[i][j]/Msc[i][i])<epsilon) break;
        }

        for(int j=i-1; j>=0; j--) //xp<x
        {
            //fill_Msc_bin_averaged(i, j, xarr, theta, Msc, kernel, add_stim);
            fill_Msc_bin_averaged(i, j, xarr, theta, Msc, KR, add_stim);

            if(abs(Msc[i][j]/Msc[i][i])<epsilon) break;
        }
    }
#ifdef OPENMP_ACTIVATED
#pragma omp barrier
#endif

    if(verbosity_scat_matrix>0)
        cout << " " + funcname + " :: done." << endl << endl;
}

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
                       bool add_stim,
                       double Te_Tg)
{
    int npx=xarr.size();

    if((int)sigarr.size()!=npx) sigarr.resize(npx);

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int i=0; i<npx; i++)
    {
        sigarr[i]=0.0;

        for(int j=0; j<npx; j++)
        {
            double stim=(add_stim ? one_minus_exp_mx(xarr[i]*Te_Tg)/one_minus_exp_mx(xarr[j]*Te_Tg) : 1.0);
            sigarr[i]+= Msc[i][j] * stim;
        }
    }

    return;
}

}

//==================================================================================================
//
// Scattering matrix representation class
//
//==================================================================================================
Msc_representation::Msc_representation(const vector<double> &xearr,
                                       const vector<double> &Int_wi,
                                       double The, double eps_thresh,
                                       int maxMom, bool stimMom)
{
    this->init(xearr, Int_wi, The, eps_thresh, maxMom, stimMom);
}

//==================================================================================================
// main setup routine directly producing sparse matrix. This version is useful for setting up
// multiple scattering matrices, where each core can handle one from the calling program
//==================================================================================================
void Msc_representation::init_serial(const vector<double> &xearr,
                                     const vector<double> &Int_wi,
                                     double The, double eps_thresh,
                                     int maxMom, bool stimMom)
{
    string func=" Msc_representation::init_serial : ";

    if(verbosity_scat_matrix>0)
        cout << func + " setting up scattering matrix for The= " << The << endl;

    this->The=The;
    this->eps_thresh=eps_thresh;

    double (*kernel)(double omega0, double omega, double theta)=thermal_kernel_exact;
    //double (*kernel)(double omega0, double omega, double theta)=thermal_kernel_SS_K;
    //double (*kernel)(double omega0, double omega, double theta)=thermal_kernel_SS_C;

    // create matrix
    int npx=xearr.size();
    vector<double> Msc_rows(npx); // rows of matrix;
    Msc_sparse.clear();

    for(int i=0; i<npx; i++) //x
    {
        // reset rows of matrix
        for(int j=0; j<npx; j++) Msc_rows[j]=0.0;

        double omega_fac=Int_wi[i] * The; // dnu' weight
        double om0=xearr[i]*The;

        // diagonal element for reference
        Msc_rows[i]= kernel(om0, om0, The) * omega_fac;

        // up-scattering wing
        for(int j=i+1; j<npx; j++) //xp>x
        {
            double omega_fac=Int_wi[j] * The; // dnu' weight
            double omp=xearr[j]*The;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc_rows[j]= kernel(om0, omp, The) * omega_fac;

            if(abs(Msc_rows[j]/Msc_rows[i])<eps_thresh) break;
        }

        // down-scattering wing
        for(int j=i-1; j>=0; j--) //xp<x
        {
            double omega_fac=Int_wi[j] * The; // dnu' weight
            double omp=xearr[j]*The;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc_rows[j]= kernel(om0, omp, The) * omega_fac;

            if(abs(Msc_rows[j]/Msc_rows[i])<eps_thresh) break;
        }

        // saving step
        for(int j=0; j<npx; j++)
            if(Msc_rows[j]!=0.0) Msc_sparse.save_info(npx, i, j, Msc_rows[j]);
    }

    // initialize moments up to some maximal order [Te=Tg]
    update_moments(xearr, maxMom, stimMom, 1.0);

    if(verbosity_scat_matrix>0) cout << func + " done." << endl << endl;

    return;
}

//==============================================================================================
// This version uses parallel setup by filling a full matrix first and then copying
//==============================================================================================
void Msc_representation::init(const vector<double> &xearr,
                              const vector<double> &Int_wi,
                              double The, double eps_thresh,
                              int maxMom, bool stimMom)
{
    string func=" Msc_representation::init :";
    if(verbosity_scat_matrix>0)
        cout << func + " setting up scattering matrix for The= " << The << endl;

    this->The=The;
    this->eps_thresh=eps_thresh;

    // create matrix
    int npx=xearr.size();
    vector<vector<double> > Msc_full(npx, vector<double>(npx, 0.0));

    // although this is first setting up the full matrix, because this is done in parallel,
    // it is faster for single matrices. Direct sparse matrix version for multiple The.
    CSpack_scattering_matrix::compute_scattering_matrix(xearr, The, Int_wi, Msc_full, "exact", eps_thresh);

    // saving step
    Msc_sparse.clear();
    for(int i=0; i<npx; i++)
        for(int j=0; j<npx; j++)
            if(Msc_full[i][j]!=0.0) Msc_sparse.save_info(npx, i, j, Msc_full[i][j]);

    // initialize moments up to some maximal order [Te=Tg]
    update_moments(xearr, maxMom, stimMom, 1.0);

    if(verbosity_scat_matrix>0) cout << func + " done." << endl << endl;

    return;
}

//==================================================================================================
void Msc_representation::update_moments(const vector<double> &xearr,
                                        int maxMom, bool stimMom, double Te_Tg)
{
    if(maxMom<0) return;

    int npx=xearr.size();

    Sigmas.clear();
    Sigmas.resize(maxMom+1, vector<double>(npx, 0.0));

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int i=0; i<npx; i++)
    {
        for(int j=0; j<npx; j++)
        {
            double Mij=Msc(i, j);

            if(Mij!=0.0)
            {
                double xgi=xearr[i]*Te_Tg, xgj=xearr[j]*Te_Tg;
                double stim=(stimMom ? one_minus_exp_mx(xgi)/one_minus_exp_mx(xgj) : 1.0);

                for(int k=0; k<=maxMom; k++) Sigmas[k][i]+=Msc(i, j) * stim * pow(xgj/xgi-1.0, k);
            }
        }
    }

    return;
}

//==================================================================================================
void Msc_representation::Get_Msc(vector<vector<double> > &Msc_full)
{
    int npx=Msc_sparse.dim;

    for(int i=0; i<npx; i++)
        for(int j=0; j<npx; j++) Msc_full[i][j]=Msc(i, j);

    return;
}

double Msc_representation::Get_Sigmak(int k, int i)
{
    if(k<(int)Sigmas.size()) return Sigmas[k][i];
    else throw_error("Msc_representation::Get_Sigmak", "requested moment not setup", 1);

    return 0.0;
}

//==================================================================================================
//
// Scattering matrix representation class for multiple values of Theta
//
//==================================================================================================
Msc_representation_Te :: Msc_representation_Te(const vector<double> &xearr,
                                               const vector<double> &Int_wi,
                                               double The_min, double The_max, int logdens_The,
                                               double eps_thresh, double eps_interpol,
                                               int maxMom, bool stimMom)
{
    init(xearr, Int_wi,
         The_min, The_max, logdens_The,
         eps_thresh, eps_interpol,
         maxMom, stimMom);
}

//==================================================================================================
// workhorse for interpolation accross The
//==================================================================================================
double Msc_representation_Te::do_interpol(double lgx, const double *lgxa, const double *ya)
{
    //===========================================================================
    // output with 4 point interpolation
    //===========================================================================
    double a[4], D[4], DD3=pow(lgxa[1]-lgxa[0], 3), r=0.0;

    for(int m=0; m<4; m++) D[m]=lgx-lgxa[m];

    a[0]= D[1]*D[2]*D[3]/(-6.0*DD3);
    a[1]= D[0]*D[2]*D[3]/( 2.0*DD3);
    a[2]= D[0]*D[1]*D[3]/(-2.0*DD3);
    a[3]= D[0]*D[1]*D[2]/( 6.0*DD3);

    for(int m=0; m<4; m++) r+=a[m]*ya[m];
    return r;
}

//==================================================================================================
void Msc_representation_Te :: init(const vector<double> &xearr,
                                   const vector<double> &Int_wi,
                                   double The_min, double The_max, int logdens_The,
                                   double eps_thresh, double eps_interpol,
                                   int maxMom, bool stimMom)
{
    this->The_min=The_min;
    this->The_max=The_max;
    this->The_curr=0.0;
    this->logdens_The=logdens_The;
    this->eps_interpol=eps_interpol;
    this->npx=xearr.size();

    npThe=init_xarr_dens(The_min, The_max, The_arr, logdens_The, 0);

    Msc_The.resize(npThe);
    Msc_full.clear();
    Msc_full.resize(xearr.size(), vector<double>(xearr.size(), 0.0));
    Msc_sparse.clear();

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int iT=0; iT<npThe; iT++)
        Msc_The[iT].init_serial(xearr, Int_wi, The_arr[iT], eps_thresh, maxMom, stimMom);

    return;    
}

//==================================================================================================
double Msc_representation_Te :: Msc(int i, int j, unsigned int iT, double The)
{
    // store interpolation data
    double lgThe[4], lgM[4];
    for(int ix=0; ix<4; ix++)
    {
        lgThe[ix]=log(The_arr[iT+ix]);

        double Mij=Msc_The[iT+ix].Msc(i, j);
        if(Mij==0.0) return 0.0; // if any of the matrix elements = 0 don't interpolate

        lgM[ix]=log(Mij);

    }

    return exp(do_interpol(log(The), lgThe, lgM));
}

//==================================================================================================
double Msc_representation_Te :: Msc(int i, int j, double The)
{
    if(npThe==0) throw_error("Msc_representation_Te :: Msc", "data not set", 1);

    bool show=1;
    bool rangeok=check_xrange(The, The_arr, show, "The");
    if(!rangeok) throw_error("Msc_representation_Te :: Msc", "The outside of range", 2);

    // find index around The
    unsigned int iT=get_start_index_interpol(The, The_arr, 4);

    return Msc(i, j, iT, The);
}

//==================================================================================================
void Msc_representation_Te :: Get_Msc(double The, vector<vector<double> > &Msc)
{
    if((int)Msc.size()!=npx) Msc.resize(npx, vector<double>(npx));

    if(fabs(The_curr/The-1.0)>eps_interpol)
    {
        The_curr=The;

        if(npThe==0) throw_error("Msc_representation_Te :: Get_Msc", "data not set", 1);

        bool show=1;
        bool rangeok=check_xrange(The, The_arr, show, "The");
        if(!rangeok) throw_error("Msc_representation_Te :: Get_Msc", "The outside of range", 2);

        // find index around The
        unsigned int iT=get_start_index_interpol(The, The_arr, 4);

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
        // interpolation across Te
        for(int i=0; i<npx; i++)
            for(int j=0; j<npx; j++)
                Msc[i][j]=this->Msc(i, j, iT, The);
#ifdef OPENMP_ACTIVATED
#pragma omp barrier
#endif
    }

    return;
}

//==================================================================================================
const ODE_solver_LA::ODE_solver_matrix& Msc_representation_Te :: Get_Msc(double The)
{
    if(fabs(The_curr/The-1.0)>eps_interpol || Msc_sparse.dim==0)
    {
        Msc_sparse.clear();

        Get_Msc(The, Msc_full);

        // save into sparse matrix
        for(int i=0; i<npx; i++)
            for(int j=0; j<npx; j++)
                if(Msc_full[i][j]!=0.0) Msc_sparse.save_info(npx, i, j, Msc_full[i][j]);
    }

    return Msc_sparse;
}

//==================================================================================================
//
// C routines [only partially setup]
//
//==================================================================================================

//==================================================================================================
// global memory that is required for C part of the routines
//==================================================================================================
vector<vector<double> > _global_Msc_temp;
vector<Kernel_representation> _global_KR_temp;

//==================================================================================================
// C versions of the scattering matrix setup
//==================================================================================================
extern "C" {

//==================================================================================================
// computes weights to turn int f(x) dx == sum Int_wi f(xi) on the grid given by xarr
//--------------------------------------------------------------------------------------------------
// inputs  :
// *xarr   : pointer to frequency grid points x=h nu/kTe = omega/theta
// npointsx: points in x
//--------------------------------------------------------------------------------------------------
// output  :
// *Int_wi : pointer to Integral weight factors to turn int f(x) dx == sum Int_wi f(xi).
//--------------------------------------------------------------------------------------------------
// comment : Memory has to be allocated before calling the function
//==================================================================================================
void Integral_weights(const double *xarr, double *Int_wi, int npointsx)
{
    vector<double> vecx(npointsx), vecw;
    for(int i=0; i<npointsx; i++) vecx[i]=xarr[i];

    //CSpack_weights::Integral_weights_trapz(vecx, vecw); // simplest rule
    CSpack_weights::Integral_weights(vecx, vecw); // 5 points rule [more accurate]

    for(int i=0; i<npointsx; i++) Int_wi[i]=vecw[i];

    return;
}

//==================================================================================================
// Routine for scattering matrix setup using explicit computation but with threshold
//--------------------------------------------------------------------------------------------------
// inputs  :
// *xarr   : pointer to frequency grid points x=h nu/kTe = omega/theta
// *Int_wi : pointer to Integral weight factors to turn int f(x) dx == sum Int_wi f(xi)
// npointsx: points in x
// theta   : kTe/mc^2
// type    : type of kernel to be used explicitly [0: 'exact', 1: 'SS_K', 2: 'SS_C']
// epsilon : optional parameter to compress matrix density [eps<1.0e-4 recommended]
//--------------------------------------------------------------------------------------------------
// outputs : Msc = wj Pij theta
//--------------------------------------------------------------------------------------------------
// comment : Memory has to be allocated before calling the function
//==================================================================================================
void compute_scattering_matrix_explicit(const double *xarr, const double *Int_wi, int npointsx,
                                        double theta, int kernel_type, double epsilon,
                                        double **Msc)
{
    vector<double> vecx(npointsx), vecw(npointsx);
    for(int i=0; i<npointsx; i++){ vecx[i]=xarr[i]; vecw[i]=Int_wi[i]; }

    string type;
    if(kernel_type==0) type="exact";
    else if(kernel_type==1) type="SS_K";
    else if(kernel_type==2) type="SS_C";

    CSpack_scattering_matrix::compute_scattering_matrix(vecx, theta, vecw,
                                                        _global_Msc_temp,
                                                        type, epsilon);

    for(int i=0; i<npointsx; i++)
        for(int j=0; j<npointsx; j++) Msc[i][j]=_global_Msc_temp[i][j];

    return;
}

//==================================================================================================
// Routine for scattering matrix setup using kernel representation routines
//--------------------------------------------------------------------------------------------------
// inputs  :
// *xarr   : pointer to frequency grid points x=h nu/kTe = omega/theta
// *Int_wi : pointer to Integral weight factors to turn int f(x) dx == sum Int_wi f(xi)
// npointsx: points in x
// nK      : defines number of points per kernel wing for Kernel-representation method
// theta   : kTe/mc^2
// epsilon : optional parameter to compress matrix density [eps<1.0e-4 recommended]
//--------------------------------------------------------------------------------------------------
// outputs : Msc = wj Pij theta
//--------------------------------------------------------------------------------------------------
// comment : Memory has to be allocated before calling the function
//==================================================================================================
void compute_scattering_matrix_KR(const double *xarr, const double *Int_wi, int npointsx,
                                  int nK, double theta, double epsilon,
                                  double **Msc)
{
    vector<double> vecx(npointsx), vecw(npointsx);
    for(int i=0; i<npointsx; i++){ vecx[i]=xarr[i]; vecw[i]=Int_wi[i]; }

    CSpack_scattering_matrix::compute_scattering_matrix(vecx, theta, vecw, nK,
                                                        _global_Msc_temp, _global_KR_temp,
                                                        epsilon, 0);

    for(int i=0; i<npointsx; i++)
        for(int j=0; j<npointsx; j++) Msc[i][j]=_global_Msc_temp[i][j];

    return;
}

}
//==================================================================================================
//==================================================================================================
