//==================================================================================================
//
// Integral weights using Lagrange interpolation polynomial
//
//==================================================================================================
// Integral int g(x) dx is written as sum_i wi gi. The function g(x) is expressed in terms of
// Lagrange polynomials. An extra weigt function f(x) can be introduced to write
// int f(x) g(x) dx = sum_i wfi gi.
//==================================================================================================
//
// Author: Jens Chluba
// first implementation: June 2021
// last modification   : June 2021
//
//==================================================================================================

#include "routines.h"
#include "Patterson.h"

#include "CSpack.h"

using namespace std;
using namespace CSpack_functions;
using namespace CSpack_kernels;

//==================================================================================================
// preamble
//==================================================================================================

//==================================================================================================
// main functions in namespace
//==================================================================================================
namespace CSpack_weights {

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
    Int_wi.clear();
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
    Int_wi.clear();
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
void Integral_weights_logx(vector<double> &xarr, vector<double> &Int_wi)
{
    vector<double> lgx=xarr;
    for(int i=0; i<(int)lgx.size(); i++) lgx[i]=log(lgx[i]);

    Integral_weights(lgx, Int_wi);

    for(int i=0; i<(int)lgx.size(); i++) Int_wi[i]*=xarr[i];

    return;
}

//==================================================================================================
double Lagrange_Polynomial(double x, int m, int order, const double *xak)
{
    // xak[0, 1, ... , order]
    // m is index of required term
    double lx=1.0;
    for(int k=0; k<=order; k++) lx*=(m!=k ? (x-xak[k])/(xak[m]-xak[k]) : 1.0);

    return lx;
}

//==================================================================================================
struct Lagrange_Polynomial_integral_data
{
    int m, order;
    const double *xak;

    // to allow including additional weight function
    void *pfw;
    double (*fw)(double x, void *pfw);

    // for checking minimal / maximal x within bin
    bool do_check_x;
    double checkmin, checkmax;

    Lagrange_Polynomial_integral_data()
    {
        m=order=0;
        xak=NULL;
        pfw=NULL;
        fw=NULL;
        do_check_x=0;
        checkmin=-1.0e+300; checkmax=1.0e+300;
    }
};

double dLagrange_Polynomial_weight(double x, void *p)
{
    Lagrange_Polynomial_integral_data &d=*(Lagrange_Polynomial_integral_data *)p;
    double fw=(d.fw==NULL ? 1.0 : d.fw(x, d.pfw));

    return fw*Lagrange_Polynomial(x, d.m, d.order, d.xak);
}

void Lagrange_Polynomial_weights(int i, int order, const vector<double> &xa,
                                 Lagrange_Polynomial_integral_data &LD,
                                 vector<double> &w,
                                 bool do_norm=1)
{
    // order == order of polynomial. This means order+1 points are used
    int nx=(int)xa.size();

    // define range round grid-point i
    double ishift=-order/2;
    while(i+ishift<0) ishift++;
    while(i+ishift+order>=nx) ishift--;

    // point to place in array that is required
    int istart=i+ishift;
    LD.xak=&xa[istart];

    // interval around bin is the same for all integrals
    double xl=(i==0    ? xa[0]     : (xa[i]+xa[i-1])/2.0);
    double xu=(i==nx-1 ? xa.back() : (xa[i]+xa[i+1])/2.0);
    double epsrel=1.0e-8, epsabs=1.0e-50;

    // check range
    double Kl=(LD.do_check_x ? LD.fw(LD.checkmin, LD.pfw) : xl);
    double Ku=(LD.do_check_x ? LD.fw(LD.checkmax, LD.pfw) : xu);
    double xll=max(xl, Kl), xul=min(xu, Ku);

    //if(xl!=xll) cout << xl << " " << xll << " L " << xul << endl;
    //if(xu!=xul) cout << xll << " " << xul << " U " << xu << endl;

    // to fix normalization [old / incorrect version]
    //double fwj=(do_norm==1 ? (LD.fw==NULL ? 1.0 : LD.fw(xa[i], LD.pfw)) : 1.0);

    for(LD.m=0; LD.m<=order; LD.m++)
    {
        // to fix normalization
        // [JC: this had a serios bug! Weight factor for each j has to be taken out...]
        double fwj=(do_norm==1 ? (LD.fw==NULL ? 1.0 : LD.fw(xa[istart+LD.m], LD.pfw)) : 1.0);

        w[istart+LD.m]+=Integrate_using_Patterson_adaptive(xll, xul, epsrel, epsabs,
                                                           dLagrange_Polynomial_weight, &LD)/fwj;
    }

    return;
}

//==================================================================================================
// compute one weight coefficient
//==================================================================================================
void Lagrange_Polynomial_weights(int i, int order,
                                 const vector<double> &xa, vector<double> &w,
                                 double (*fw)(double x, void *pfw), void *pfw,
                                 bool do_norm)
{
    Lagrange_Polynomial_integral_data LD;
    LD.order=order;
    LD.fw=fw;
    LD.pfw=pfw;
    LD.do_check_x=1;

    Lagrange_Polynomial_weights(i, order, xa, LD, w, do_norm);
}

//==================================================================================================
// explicit computation of weights
//==================================================================================================
void compute_all_Lagrange_Polynomial_weights(int order, const vector<double> &xa, vector<double> &w)
{
    int nx=(int)xa.size();
    w.clear();
    w.resize(nx, 0.0);

    Lagrange_Polynomial_integral_data LD;
    LD.order=order;

    for(int i=0; i<nx; i++) Lagrange_Polynomial_weights(i, order, xa, LD, w, 0);

    return;
}

//==================================================================================================
// compute weights using extra weight function
//==================================================================================================
void compute_all_Lagrange_Polynomial_weights(int order, const vector<double> &xa, vector<double> &w,
                                             double (*fw)(double x, void *pfw), void *pfw)
{
    int nx=(int)xa.size();
    w.clear();
    w.resize(nx, 0.0);

    Lagrange_Polynomial_integral_data LD;
    LD.order=order;
    LD.fw=fw;
    LD.pfw=pfw;

    for(int i=0; i<nx; i++) Lagrange_Polynomial_weights(i, order, xa, LD, w, 1);

    return;
}

}

//==================================================================================================
//==================================================================================================

