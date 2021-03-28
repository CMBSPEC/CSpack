//==================================================================================================
//  Created by Abir Sarkar on 15/11/2019 and modified by JC. These functions are based on
//  Sarkar, Chluba and Lee, MNRAS, 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#include "routines.h"

#include "CSpack.h"
#include "ODE_solver_LA.h"

using namespace std;
using namespace CSpack_functions;
using namespace CSpack_kernels;

//==================================================================================================
// main functions in namespace
//==================================================================================================
namespace CSpack_scattering_matrix {

int verbosity_scat_matrix=0;

void set_verbosity(int verb){ verbosity_scat_matrix=verb; return; }

//==================================================================================================
// weights using f(x) = L(x) [5 point Lagrange polynomial] and integrating over
// x in [xa, xb] around xj
//==================================================================================================
double Int_li_xj_II(double x1, double x2, double x3, double x4, double xj,
                    double xa, double xb)
{
    double r=(pow(xb, 5)-pow(xa, 5))/5.0;
    r+=-(x1+x2+x3+x4)/4.0 * (pow(xb, 4)-pow(xa, 4));
    r+=(x1*(x2+x3+x4)+x2*(x3+x4)+x3*x4)/3.0 * (pow(xb, 3)-pow(xa, 3));
    r+=-(x1*(x2*x3+x2*x4+x3*x4)+x2*x3*x4)/2.0 * (xb-xa) * (xb+xa);
    r+= x1*x2*x3*x4 * (xb-xa);

    return r /(x1-xj)/(x2-xj)/(x3-xj)/(x4-xj);
}

double Int_li_xj_II(int j, int k, vector<double> &xi)
{
    //    double xa=(xi[j]+xi[j-1])/2.0, xb=(xi[j+1]+xi[j])/2.0;
    double xa=xi[j], xb=xi[j+1];
    double x1, x2, x3, x4;

    if(k==-2){ x1=xi[j-1]; x2=xi[j]; x3=xi[j+1]; x4=xi[j+2]; }
    else if(k==-1){ x1=xi[j-2]; x2=xi[j]; x3=xi[j+1]; x4=xi[j+2]; }
    else if(k== 0){ x1=xi[j-2]; x2=xi[j-1]; x3=xi[j+1]; x4=xi[j+2]; }
    else if(k== 1){ x1=xi[j-2]; x2=xi[j-1]; x3=xi[j]; x4=xi[j+2]; }
    else if(k== 2){ x1=xi[j-2]; x2=xi[j-1]; x3=xi[j]; x4=xi[j+1]; }
    else { return 0.0; }

    double r=Int_li_xj_II(x1, x2, x3, x4, xi[j+k], xa, xb);

    return r;
}

//==================================================================================================
void Integral_weights_trapz(vector<double> &xarr, vector<double> &Int_wi)
{
    int np=xarr.size();
    Int_wi.resize(np, 0.0);

    //================================================================
    // initial points around lower boundary
    //================================================================
    int ix=0;
    Int_wi[ix]+=(xarr[ix+1]-xarr[ix])/2.0; ix++;
    for(; ix<np-1; ix++) Int_wi[ix]+=(xarr[ix+1]-xarr[ix-1])/2.0;
    Int_wi[ix]+=(xarr[ix]-xarr[ix-1])/2.0;

    return;
}

//==================================================================================================
void Integral_weights(vector<double> &xarr, vector<double> &Int_wi)
{
    int np=xarr.size();
    Int_wi.resize(np, 0.0);

    //================================================================
    // initial points around lower boundary
    //================================================================
    int ix=0;
    Int_wi[ix]+=(xarr[ix+1]-xarr[ix])/2.0; ix++;
    for(; ix<2; ix++) Int_wi[ix]+=(xarr[ix+1]-xarr[ix-1])/2.0;
    Int_wi[ix]+=(xarr[ix]-xarr[ix-1])/2.0;

    //================================================================
    // internal points with 5-point formula
    //================================================================
    for(ix=2; ix<np-3; ix++)
        for(int k=-2; k<3; k++) Int_wi[ix+k]+=Int_li_xj_II(ix, k, xarr);

    //================================================================
    // finish off using simple trapeziodal rule
    //================================================================
    Int_wi[ix]+=(xarr[ix+1]-xarr[ix])/2.0; ix++;
    for(; ix<np-1; ix++) Int_wi[ix]+=(xarr[ix+1]-xarr[ix-1])/2.0;
    Int_wi[ix]+=(xarr[ix]-xarr[ix-1])/2.0;

    return;
}

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
    double omin=xarr[0]*theta/300.0, omax=xarr.back()*theta*300.0;
    //double omin=xarr[0]*theta, omax=xarr.back()*theta;
    for(int i=0; i<npx; i++) KR[i].allocate_splines(nK);

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int i=0; i<npx; i++)
        KR[i].init(omin, xarr[i]*theta, omax, nK, theta, epsilon, epsilon, 2, stim);

    for(int i=0; i<npx; i++) //x
    {
        // reset matrix
        for(int j=0; j<npx; j++) Msc[i][j]=0.0;

        double omega_fac=Int_wi[i] * theta; // dnu' weight
        double om0=xarr[i]*theta;
        // diagonal element for reference
        Msc[i][i]= KR[i].Kernel(om0) * omega_fac;

        for(int j=i+1; j<npx; j++) //xp>x
        {
            double omega_fac=Int_wi[j] * theta; // dnu' weight
            double omp=xarr[j]*theta;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc[i][j]= KR[i].Kernel(omp) * omega_fac;

            if(abs(Msc[i][j]/Msc[i][i])<epsilon) break;
        }

        for(int j=i-1; j>=0; j--) //xp<x
        {
            double omega_fac=Int_wi[j] * theta; // dnu' weight
            double omp=xarr[j]*theta;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc[i][j]= KR[i].Kernel(omp) * omega_fac;

            if(abs(Msc[i][j]/Msc[i][i])<epsilon) break;
        }
    }

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

        double omega_fac=Int_wi[i] * theta; // dnu' weight
        double om0=xarr[i]*theta;
        // diagonal element for reference
        Msc[i][i]= kernel(om0, om0, theta) * omega_fac;

        for(int j=i+1; j<npx; j++) //xp>x
        {
            double omega_fac=Int_wi[j] * theta; // dnu' weight
            double omp=xarr[j]*theta;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc[i][j]= kernel(om0, omp, theta) * omega_fac;

            if(abs(Msc[i][j]/Msc[i][i])<epsilon) break;
        }

        for(int j=i-1; j>=0; j--) //xp<x
        {
            double omega_fac=Int_wi[j] * theta; // dnu' weight
            double omp=xarr[j]*theta;

            // P(nu-->nu') here for all pairs i == col and j == row
            Msc[i][j]= kernel(om0, omp, theta) * omega_fac;

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
                                       string type, int maxMom, bool stimMom)
{
    this->init(xearr, Int_wi, The, eps_thresh, type, maxMom, stimMom);
}

//==================================================================================================
// main setup routine directly producing sparse matrix. This version is useful for setting up
// multiple scattering matrices, where each core can handle one from the calling program
//==================================================================================================
void Msc_representation::init_serial(const vector<double> &xearr,
                                     const vector<double> &Int_wi,
                                     double The, double eps_thresh,
                                     string type, int maxMom, bool stimMom)
{
    string func=" Msc_representation::init_serial : ";

    if(verbosity_scat_matrix>0)
        cout << func + " setting up scattering matrix for The= " << The << endl;

    this->The=The;
    this->eps_thresh=eps_thresh;

    double (*kernel)(double omega0, double omega, double theta)=NULL;

    if(type=="exact") kernel=thermal_kernel_exact;
    else if(type=="SS_K") kernel=thermal_kernel_SS_K;
    else if(type=="SS_C") kernel=thermal_kernel_SS_C;
    else throw_error("Msc_representation::init_serial", "kernel type not available", 1);

    // create matrix
    int npx=xearr.size();
    vector<double> Msc_rows(npx); // rows of matrix;

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
                              string type, int maxMom, bool stimMom)
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
    CSpack_scattering_matrix::compute_scattering_matrix(xearr, The, Int_wi, Msc_full, type, eps_thresh);

    // saving step
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
    if(k<Sigmas.size()) return Sigmas[k][i];
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
                                               string type, int maxMom, bool stimMom)
{
    init(xearr, Int_wi,
         The_min, The_max, logdens_The,
         eps_thresh, eps_interpol,
         type, maxMom, stimMom);
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
                                   string type, int maxMom, bool stimMom)
{
    this->The_min=The_min;
    this->The_max=The_max;
    this->logdens_The=logdens_The;
    this->eps_interpol=eps_interpol;

    npThe=init_xarr(The_min, The_max, The_arr, logdens_The, 1);

    Msc_The.resize(npThe);
    Msc_full.resize(xearr.size(), vector<double>(xearr.size(), 0.0));
    Msc_sparse.clear();

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for(int iT=0; iT<npThe; iT++)
        Msc_The[iT].init_serial(xearr, Int_wi, The_arr[iT], eps_thresh, type, maxMom, stimMom);

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
        lgM  [ix]=log(Msc_The[iT+ix].Msc(i, j));

        if(lgM[ix]==0.0) return 0.0; // if any of the matrix elements = 0 don't interpolate
    }

    return exp(do_interpol(log(The), lgThe, lgM));
}

//==================================================================================================
double Msc_representation_Te :: Msc(int i, int j, double The)
{
    if(npThe==0) throw_error("Msc_representation_Te :: Msc", "data not set", 1);

    bool show=0;
    bool rangeok=check_xrange(The, The_arr, show, "The");
    if(!rangeok) throw_error("Msc_representation_Te :: Msc", "The outside of range", 2);

    // find index around The
    unsigned int iT=get_start_index_interpol(The, The_arr, 4);

    return Msc(i, j, iT, The);
}

//==================================================================================================
const ODE_solver_LA::ODE_solver_matrix& Msc_representation_Te :: Get_Msc(double The)
{
    if(fabs(The_curr/The-1.0)>eps_interpol)
    {
        The_curr=The;
        Msc_sparse.clear();

        if(npThe==0) throw_error("Msc_representation_Te :: Msc", "data not set", 1);

        bool show=0;
        bool rangeok=check_xrange(The, The_arr, show, "The");
        if(!rangeok) throw_error("Msc_representation_Te :: Msc", "The outside of range", 2);

        // find index around The
        unsigned int iT=get_start_index_interpol(The, The_arr, 4);

        int npx=Msc_full.size();

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
        // interpolation across Te
        for(int i=0; i<npx; i++)
            for(int j=0; j<npx; j++)
                Msc_full[i][j]=Msc(i, j, iT, The);
#ifdef OPENMP_ACTIVATED
#pragma omp barrier
#endif

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

    //CSpack_scattering_matrix::Integral_weights_trapz(vecx, vecw); // simplest rule
    CSpack_scattering_matrix::Integral_weights(vecx, vecw); // 5 points rule [more accurate]

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
