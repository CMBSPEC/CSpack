//==================================================================================================
// Author: Jens Chluba
// first implementation: Jan 2010
// last modification   : Aug 2021
//==================================================================================================
// Aug  5 th, 2021: made integrator backward compatible
// Mar  30th, 2020: generalized integrator to allow for parameters.
// Oct  20th, 2011: One can now give pmax (refinement level) to the integrator.
// June 14th, 2011: Nested table will only be updated if more abscissae and weights are needed.

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>

#include "Definitions.h"
#include "Chebyshev_Int.h"
#include "routines.h"

using namespace std;

//==================================================================================================
//
// parameters for integrator
//
//==================================================================================================
static int Chebyshev_pmax=16;
static bool Chebyshev_verbosity=1;

void no_Chebyshev_alert(){ Chebyshev_verbosity=0; }
void Chebyshev_alert(){ Chebyshev_verbosity=1; }

//==================================================================================================
//
// simple test functions
//
//==================================================================================================
double Chebyshev_test_fcn(double x, void *pars)
{
    return 1.0/(1.0+0.1*sin(10.0*x));
    //return 1.0/(1.0+x*x);
    //return 1.0/(exp(x)-1.0)/x;
    //return x*x*x*x/(exp(x)-1.0);
    //return x*x*x/(exp(x)-1.0);
    //return x*x/(exp(x)-1.0);
    //return x*x;
}

double Chebyshev_test_fcn_log(double lgx, void *pars)
{
    double x=exp(lgx);
    return x*Chebyshev_test_fcn(x, pars);
}

//==================================================================================================
//
// functions for mapping of f(x) to interval mu=(-1, 1)
//
//==================================================================================================
struct Chebyshev_vars
{
    double lambda, a, b, c;
    void *pars;

    double (*fptr)(double);
    double (*fptr_pars)(double, void *pars);

    double eval_func(double x)
    { return fptr_pars==NULL ? fptr(x) : fptr_pars(x, pars); }

    Chebyshev_vars()
    {
        pars=NULL;
        fptr=NULL;
        fptr_pars=NULL;
    }
};

//==================================================================================================
double Chebyshev_Int_f1(double mu, Chebyshev_vars& v)
{
    double dx_dmu=v.c/pow(v.a-mu, 2);
    double x=v.lambda*(v.b+mu)/(v.a-mu);
    return dx_dmu*v.eval_func(x);
}

//==================================================================================================
double Chebyshev_Int_f2(double mu, Chebyshev_vars& v)
{
    double dx_dmu=v.a;
    double x=v.a*mu+v.b;
    return dx_dmu*v.eval_func(x);
}

//==================================================================================================
//
// computations of weights & abscissae; little optimization was done here
//
//==================================================================================================
vector<vector<double> > Chebyshev_xi, Chebyshev_wi;

//==================================================================================================
double sum_Chebyshev_weights(int n, double theta_k)
{
    double r=0.0;
    for(int m=(n+1)/2; m>=1; m--) r+=sin((2*m-1)*theta_k)/(2*m-1);
    return r*4.0*sin(theta_k)/(n+1);
}

//==================================================================================================
void compute_nested_Chebyshev_table(int p, vector<vector<double> > &xi, vector<vector<double> > &wi,
                                    int mess=0)
{
    //==================================================================================
    // checks whether computation is necessary
    //==================================================================================
    if(p<=(int)wi.size()-1) return;

    //==================================================================================
    vector<double> dum;
    double theta_k, tt;

    if(mess>=1) cout << " Updating nested Table for Chebychev integration rules with p= " << p << endl;

    //==================================================================================
    // empty element for p=0
    //==================================================================================
    if(wi.size()==0)
    {
        wi.clear();
        xi.clear();
        wi.push_back(dum);
        xi.push_back(dum);
    }

    //==================================================================================
    for(int pl=(int)wi.size(); pl<=p; pl++)
    {
        int nmax=(int)pow(2.0, pl)-1;

        if(mess>=1)
        {
            cout.precision(16);
            cout << "\n nmax " << nmax << endl;
            cout << " wi " << endl;
        }

        //==============================================================================
        dum.clear();
        for(int k=1; k<nmax/2+2; k++)
        {
            theta_k=PI*k/(nmax+1);
            tt=sum_Chebyshev_weights(nmax, theta_k);
            dum.push_back(tt);
        }

        if(mess>=1) for(int i=0; i<(int)dum.size(); i++) cout << dum[i] << ", ";
        wi.push_back(dum);

        if(mess>=1) cout << "\n xi " << endl;

        //==============================================================================
        dum.clear();
        if(pl==1) dum.push_back(0.0);
        else
        {
            for(int k=1; k<nmax/2+1; k+=2)
            {
                theta_k=PI*k/(nmax+1);
                tt=cos(theta_k);
                dum.push_back(tt);
            }
        }

        if(mess>=1) for(int i=0; i<(int)dum.size(); i++) cout << dum[i] << ", ";
        xi.push_back(dum);

        if(mess>=1) wait_f_r("\n");
    }
    if(mess>=1) wait_f_r("\n");
    //==================================================================================

    return;
}

//==================================================================================================
void setup_nested_Chebyshev_table(int pmax, int mess)
{
    for(int p=0; p<=pmax; p++)
        compute_nested_Chebyshev_table(p, Chebyshev_xi, Chebyshev_wi, mess);
    
    return;
}

//==================================================================================================
//
// Functions to compute the integral
//
//==================================================================================================
void compute_integral_function_Chebyshev(int p, Chebyshev_vars& v,
                                         vector<double> &Chebyshev_fvals,
                                         double (*f)(double, Chebyshev_vars&),
                                         int *neval, double *r)
{
    int debug=0;

    //==============================================================================
    *neval=(int)pow(2.0, p)-1;
    *r=0.0;
    compute_nested_Chebyshev_table(p, Chebyshev_xi, Chebyshev_wi);

    //==============================================================================
    // only add new points
    //==============================================================================
    int index_f=0, nxi=Chebyshev_xi[p].size();
    vector<double> loc_fvals(2*nxi);

    //==============================================================================
    // copy first point
    //==============================================================================
    double mu=Chebyshev_xi[p][0];
    loc_fvals[0]=f(mu, v)+f(-mu, v);

    //==============================================================================
    // integral points
    //==============================================================================
    for(int k=1; k<nxi; k++)
    {
        //------------------------------------------
        // copy old point
        //------------------------------------------
        loc_fvals[2*k-1]=Chebyshev_fvals[index_f++];

        //------------------------------------------
        // calc new point
        //------------------------------------------
        mu=Chebyshev_xi[p][k];
        loc_fvals[2*k]=f(mu, v)+f(-mu, v);
    }

    //==============================================================================
    // finish last point (computed f(mu=0) outside)
    //==============================================================================
    loc_fvals[2*nxi-1]=Chebyshev_fvals[index_f];

    //==============================================================================
    // save things for next call
    //==============================================================================
    Chebyshev_fvals=loc_fvals;
    loc_fvals.clear();

    //==============================================================================
    // compute integral
    //==============================================================================
    *r=DC_sumprod(Chebyshev_fvals, Chebyshev_wi[p]);

    if(debug>=1) cout << "\n order " << p << " neval= " << *neval << " " << *r << endl;
    if(debug>=2) wait_f_r();

    return;
}

//==================================================================================================
// Driver that loops over the different Chebyshev rules
//
// a : lower bound of integral
// b : upper bound of integral (constraint: a<b)
//
// epsrel: requested relative accuracy
// epsabs: requested absolute accuracy
//
// *fptr : pointer to function f(x, pars) [pars can be set optionally]
// *neval: on exit it will contain the number of function evals
// *r: on exit it will be equal to the best estimate of the integral
//
// *pars: user parameters to be passed on in function call. Inside the function a recast is assumed
//  pmax: explicitly state maximal level of 2^pmax for function evaluations
//
// mess: optional message label for debugging
//==================================================================================================
int Integrate_using_Chebyshev(double a, double b, int pmax, double epsrel, double epsabs,
                              Chebyshev_vars &CV, int *neval, double *r, string mess)
{
    int debug=0;
    bool error_check=1;
    bool show_error_mess=Chebyshev_verbosity;
    
    if(pmax<0){ error_check=0; pmax=-pmax; }
    if(epsrel<0){ show_error_mess=0; epsrel=-epsrel; }
        
        //==============================================================================================
        *r=0.0;
    *neval=0;
    if(a>=b) return 2;
    
    if(debug>=2) cout << " Entering chebyshev integrator xl= " << a << " xu= " << b << endl;
    
    //==============================================================================================
    // parameters for transformation formula
    //==============================================================================================
    CV.a=(b-a)*0.5;
    CV.b=(b+a)*0.5;
    
    if(debug>=3) cout << " a = " << CV.a << " b= " << CV.b << endl;
    
    //==============================================================================================
    double fi=Chebyshev_Int_f2(0.0, CV);
    vector<double> fvals(1, fi);
    
    //==============================================================================
    // Integral for 1-point formula
    //==============================================================================
    *r=2.0*fi;
    *neval=1;
    if(debug>=2) cout << " f(mu=0.0)= " << fi << endl << endl;
    
    double r1, Dr;
    for(int it=2; it<=pmax; it++)
    {
        compute_integral_function_Chebyshev(it, CV, fvals, Chebyshev_Int_f2, neval, &r1);
        Dr= *r-r1;
        *r=r1;
        
        //==========================================================================
        // check error
        //==========================================================================
        if(fabs(Dr)<=max(epsabs, fabs( *r)*epsrel) && error_check)
        {
            if(debug>=1) cout << " result= " << *r << " neval= " << *neval << endl;
            
            return 0;
        }
    }
    
    if(!error_check) return 1;
    
    //==============================================================================================
    // reaching below here means the intgration failed
    //==============================================================================================
    if(show_error_mess)
        cout << " Integral did not converge using Chebyshev rules up to pmax= " << pmax << "."
             << "\n Other parameters/variables were " << a << " " << b
             << " " << *r << " " << *neval << " mess= " << mess << endl;

    return 1;
}

//==================================================================================================
int Integrate_using_Chebyshev(double a, double b, int pmax, double epsrel, double epsabs,
                              double (*fptr)(double, void *pars), void *pars,
                              int *neval, double *r, string mess)
{
    Chebyshev_vars CV;
    CV.fptr_pars=fptr;
    CV.pars=pars;

    return Integrate_using_Chebyshev(a, b, pmax, epsrel, epsabs, CV, neval, r, mess);
}

int Integrate_using_Chebyshev(double a, double b, double epsrel, double epsabs,
                              double (*fptr)(double, void *pars), void *pars,
                              int *neval, double *r, string mess)
{
    return Integrate_using_Chebyshev(a, b, Chebyshev_pmax, epsrel, epsabs, fptr, pars, neval, r, mess);
}

//==================================================================================================
int Integrate_using_Chebyshev(double a, double b, int pmax, double epsrel, double epsabs,
                              double (*fptr)(double, void *pars),
                              int *neval, double *r, string mess)
{
    Chebyshev_vars CV;
    CV.fptr_pars=fptr;

    return Integrate_using_Chebyshev(a, b, pmax, epsrel, epsabs, CV, neval, r, mess);
}

int Integrate_using_Chebyshev(double a, double b, double epsrel, double epsabs,
                              double (*fptr)(double, void *pars),
                              int *neval, double *r, string mess)
{
    return Integrate_using_Chebyshev(a, b, Chebyshev_pmax, epsrel, epsabs, fptr, neval, r, mess);
}

//==================================================================================================
int Integrate_using_Chebyshev(double a, double b, int pmax, double epsrel, double epsabs,
                              double (*fptr)(double), int *neval, double *r, string mess)
{
    Chebyshev_vars CV;
    CV.fptr=fptr;

    return Integrate_using_Chebyshev(a, b, pmax, epsrel, epsabs, CV, neval, r, mess);
}

int Integrate_using_Chebyshev(double a, double b, double epsrel, double epsabs,
                              double (*fptr)(double), int *neval, double *r, string mess)
{
    return Integrate_using_Chebyshev(a, b, Chebyshev_pmax, epsrel, epsabs, fptr, neval, r, mess);
}

//==================================================================================================
// Driver that loops over the different Chebyshev rules
//
// a : lower bound of integral
// x0: point inside the integration region (a< x0 < (a+b)/2); it allows
//     deforming the integration domain. If set outside this range x0=(3*a+b)/4.
// b : upper bound of integral (constraint: a<b)
//
// epsrel: requested relative accuracy
// epsabs: requested absolute accuracy
//
// *fptr : pointer to function f(x, pars) [pars can be set optionally]
// *neval: on exit it will contain the number of function evals
// *r: on exit it will be equal to the best estimate of the integral
//
// *pars: user parameters to be passed on in function call. Inside the function a recast is assumed
//
// mess: optional message label for debugging
//==================================================================================================
int Integrate_using_Chebyshev(double a, double x0, double b, double epsrel, double epsabs,
                              Chebyshev_vars &CV, int *neval, double *r, string mess)
{
    int debug=0;
    bool show_error_mess=Chebyshev_verbosity;
    
    if(epsrel<0){ show_error_mess=0; epsrel=-epsrel; }

    //==============================================================================================
    *r=0.0;
    *neval=0;
    if(a>=b) return 2;

    if(debug>=2) cout << " Entering chebyshev integrator xl= " << a << " xu= " << b
                      << " x0= " << x0 << endl;

    //==============================================================================================
    if(x0<=a || x0>=(a+b)/2.0*(1.0-1.0e-2) ) x0=(3.0*a+b)/4.0;            // in this case Lambda ~ 1

    //==============================================================================================
    // parameters for transformation formula
    //==============================================================================================
    CV.a=(b-a)/(a+b-2.0*x0);
    CV.b=x0*(b-a)/((a+b)*x0-2.0*a*b);
    CV.lambda=((a+b)*x0-2.0*a*b)/(a+b-2.0*x0);  // == a/b
    CV.c=CV.lambda*(CV.a+CV.b);

    if(debug>=3) cout << " a = " << CV.a << " b= " << CV.b
                      << " c= "  << CV.c << " lam= " << CV.lambda
                      << " x0= " << x0 << endl;

    if(CV.lambda==0.0) cout << " Integrate_using_Chebyshev:: check x0 " << x0
                           << " a= " << a << " b= " << b << endl;

    //==============================================================================================
    double fi=Chebyshev_Int_f1(0.0, CV);
    vector<double> fvals(1, fi);

    //==============================================================================
    // Integral for 1-point formula
    //==============================================================================
    *r=2.0*fi;
    *neval=1;
    if(debug>=2) cout << " f(mu=0.0)= " << fi << endl << endl;

    double r1, Dr;
    for(int it=2; it<=Chebyshev_pmax; it++)
    {
        compute_integral_function_Chebyshev(it, CV, fvals, Chebyshev_Int_f1, neval, &r1);
        Dr= *r-r1;
        *r=r1;

        //==========================================================================
        // check error
        //==========================================================================
        if(fabs(Dr)<=max(epsabs, fabs( *r)*epsrel))
        {
            if(debug>=1) cout << " result= " << *r << " neval= " << *neval << endl;

            return 0;
        }
    }

    //==============================================================================================
    // reaching below here means the intgration failed
    //==============================================================================================
    if(show_error_mess)
        cout << " Integral did not converge using Chebyshev rules up to pmax= " << Chebyshev_pmax << "."
             << "\n Other parameters/variables were " << a << " " << x0 << " " << b
             << " " << *r << " " << *neval << " mess= " << mess << endl;

    return 1;
}

//==================================================================================================
int Integrate_using_Chebyshev(double a, double x0, double b, double epsrel, double epsabs,
                              double (*fptr)(double, void *pars), void *pars,
                              int *neval, double *r, string mess)
{
    Chebyshev_vars CV;
    CV.fptr_pars=fptr;
    CV.pars=pars;

    return Integrate_using_Chebyshev(a, x0, b, epsrel, epsabs, CV, neval, r, mess);
}

//==================================================================================================
int Integrate_using_Chebyshev(double a, double x0, double b, double epsrel, double epsabs,
                              double (*fptr)(double, void *pars),
                              int *neval, double *r, string mess)
{
    Chebyshev_vars CV;
    CV.fptr_pars=fptr;

    return Integrate_using_Chebyshev(a, x0, b, epsrel, epsabs, CV, neval, r, mess);
}

//==================================================================================================
int Integrate_using_Chebyshev(double a, double x0, double b, double epsrel, double epsabs,
                              double (*fptr)(double),
                              int *neval, double *r, string mess)
{
    Chebyshev_vars CV;
    CV.fptr=fptr;

    return Integrate_using_Chebyshev(a, x0, b, epsrel, epsabs, CV, neval, r, mess);
}

//==================================================================================================
//==================================================================================================
