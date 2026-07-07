//======================================================================================
// Author: Jens Chluba
// first implementation: Jan 2002
// last modification: Feb 2018
//
// Purpose: collection of several simple routines
//======================================================================================
// Feb 2018: added Gaussian function (why not earlier??) and tanh-switch
// Oct 2013: added Ei(x) function
// Jul 2012: fixed a memory problem with gsl-1.15 compatibility of splines.
// Jun 2012: fixed a memory issue with spline setup routines.
// Jan 2012: added simple routines to load tables of data and create splines
// Dec 2011: added routines for Gamma and incomplete Gamma functions

#ifndef ROUTINES_H
#define ROUTINES_H

#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include <tuple>
#include <functional>
#include <algorithm>

#include <gsl/gsl_spline.h>

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
//#include <boost/math/interpolators/cubic_b_spline.hpp>

using namespace std;

//======================================================================================
// Special functions
//======================================================================================
double Gaussian(double x, double x0, double sigma);     // simple normalized Gaussian
double tanh_switch(double x, double x0, double sigma);  // simple tanh-switch

double Dawson_Int(const double &x); // Dawson integral; Based on Antia
double erf_JC(const double &x);     // Error function for a real argument
double erfc_JC(const double &x);    // complementary Error function for a real argument
double Gamma_JC(const double &x);                    // Gamma-function
// JC: fixed switching of variables [Feb 2023]
double Gamma_JC(const double &x, const double &a);   // incomplete Gamma-function
double Ei_JC(double x);             // Ei(x) = - E_1(x) = - int_{-x}^infty e^{-t}/t dt
double ex_Ei_JC(double x);          // exp(x) Ei(x)

double scaled_BesselK1(double x);   // modified Bessel function exp(x) K1(x)
double scaled_BesselK2(double x);   // modified Bessel function exp(x) K2(x)

double hyper_geom_2F1(double a, double b, double c, double x);

double max(double a, double b);
double min(double a, double b);

double factorial(int n);
double log10factorial(int n);
double log10factorial_full(int n);
double factorial_corrfac(int n);


double one_minus_exp_mx(double x);                   // == 1-exp(-x)
double nbb_func(double x);                           // == 1/(exp(x)-1)
double nbbp1_func(double x);                         // == 1+nbb(x)

//======================================================================================
// approximations for (1+del, alpha)-1 for alpha>0
//======================================================================================
double special_func_photon_integrals(double del, double alpha);

//======================================================================================
// checking for nan
//======================================================================================
bool isnan_JC(double a);

// simple error message
void throw_error(string funcname, string message, int k=1);
void throw_error(string funcname, string message, double val, int k=1);

//======================================================================================
// routines for interpolation; based in GSL
//======================================================================================

double Get_xmin_spline(int memindex, string mess);
double Get_xmax_spline(int memindex, string mess);
vector<double> Get_xmin_max_spline(int memindex, string mess);

int calc_spline_coeffies_JC(int nxi, const double *za, const double *ya,
                            string variable="");

int set_trivial_spline(string variable="");

void update_spline_coeffies_JC(int memindex, int nxi,
                               const double *za, const double *ya,
                               string variable="");

double Get_xmin_spline(int memindex, string mess="");
double Get_xmax_spline(int memindex, string mess="");
vector<double> Get_xmin_max_spline(int memindex, string mess="");

double calc_spline_JC(double x, int memindex, string mess="");

double calc_spline_deriv_JC(double x, int memindex, string mess="");

void free_spline_JC(int memindex, string mess="");

void free_all_splines_JC(int verbosity=0);

void show_spline_memory();


//======================================================================================
// load tables of data (added Jan 2012)
//======================================================================================
vector<vector<double> > load_data_from_file(string fname, int cols, bool showmess=0);

void load_data_from_file(string fname, int cols,
                         vector<int> &spline_mem_indices,
                         bool logx, bool logy, bool showmess=0);

void load_data_from_file_loglin(string fname, int cols, vector<int> &spline_mem_indices);
void load_data_from_file_loglog(string fname, int cols, vector<int> &spline_mem_indices);


//======================================================================================
// npol-1 is degree of the interpolating polynomial
//======================================================================================
void polint_JC(const double *xa, const double *ya, int na, const double x, int npol,
               double *y, double *dy);

void polint_JC(const double *xa, const double *ya, int na, const double x, int &istart,
               int npol, double *y, double *dy);

//======================================================================================
// simple grid setup
//======================================================================================
void init_xarr(double x0, double xm, double *xarr, int npts, int method_flag, int mess_flg);
void init_xarr(double x0, double xm, double *xarr, int npts, int method_flag);

int init_xarr(double x0, double xm, vector<double> &xarr,
              int npts, int method_flag, bool mess);

int init_xarr_dens(double x0, double xm, vector<double> &xarr, int logdens, bool mess=0);

//======================================================================================
void wait_f_r();
void wait_f_r(string mess);
void wait_f_r(int num);
void wait_f_r(double num);

//======================================================================================
void locate_JC(const double xx[], unsigned long n, double x, unsigned long *j);
void hunt(const double xx[], unsigned long n, double x, unsigned long *jlo);

//======================================================================================
// i/o modules
//======================================================================================
string int_to_string(int i);
string int_to_string(int i, int ni);

//======================================================================================
// root-finding methods
//======================================================================================
double find_root(double (* func)(double *), double x1, double x2, double xacc);
double find_root_brent(double (* func)(double *), double x1, double x2, double xacc);
double find_root_brent(double (* func)(double *, void *p), void *p,
                       double x1, double x2, double xacc, string call="");

//======================================================================================
// for xmgrace output
//======================================================================================
void plot_xy_function_linear(const vector<double> &xarr, const vector<double> &yarr);
void plot_xy_function_log   (const vector<double> &xarr, const vector<double> &yarr);
void plot_xy_function_loglog(const vector<double> &xarr, const vector<double> &yarr);

//======================================================================================
// Wigner 3J symbol (added 17.05.2011)
//======================================================================================
double Wigner_3J_Symbol(int j1, int j2, int j3, int m1, int m2, int m3);

//======================================================================================
// divide and conquer sum
//======================================================================================
double DC_sum(const double *yarr, int M);
double DC_sum(const vector<double> &yarr);

//======================================================================================
// divide and conquer sum of product
//======================================================================================
double DC_sumprod(const double *yarr, const double *zarr, int M);
double DC_sumprod(const vector<double> &yarr, const vector<double> &zarr);


//======================================================================================
// for checking, creating and deleting directories
//======================================================================================
bool check_if_dir_exists(string dir);
void create_directory_if_it_does_not_exist(string dir, bool message=1);
void change_to_executable(string fname);

//======================================================================================
complex<double> LnGamma_func(complex<double> z);  // complex number z

// Lambda == Gamma(a) Gamma(b) / Gamma(c) Gamma(d)
complex<double> LnLambda_func(complex<double> a, complex<double> b,
                              complex<double> c, complex<double> d);

double Harmonic_number_Re(double x);              // Re[H(i x)] is computed, x>0
double Harmonic_number_Im(double x);              // Im[H(i x)] is computed, x>0


//======================================================================================
// Check variable range
//======================================================================================
bool check_xrange  (double x  , const vector<double> &xa  , bool show=1, string mess="x");
bool check_lgxrange(double lgx, const vector<double> &lgxa, bool show=1, string mess="x");

//======================================================================================
// get index with x[i] < x assuming there is a constant grid
// grid is in ascending order
//======================================================================================
unsigned int get_start_index_interpol(double x, const vector<double> &xa, int npol);

// ========================================================================== //
// ============================== TK Functions ============================== //
// ========================================================================== //

double polyval(double x, vector<double> coeffs);

tuple<double, double> find_parabolic_peak(double x1, double y1, double x2, double y2, double x3, double y3);

tuple<double, double> find_parabolic_peak(tuple<double, double> p1, tuple<double, double> p2, tuple<double, double> p3);

vector<size_t> arg_local_extrema(vector<double> x_vals, vector<double> y_vals);

tuple<vector<double>, vector<double>> find_local_extrema(vector<double> x_vals, vector<double> y_vals);

vector<double> apply_func(const vector<double>&, function<double(double)>);

vector<double> log_arr(vector<double> my_arr);

vector<double> log_abs_arr(vector<double> my_arr);

vector<double> exp_arr(vector<double> my_arr);

vector<double> reverse(const vector<double>& my_arr);

bool does_file_exist(string fname);

void output_1d_data(vector<double> xvals, vector<double> yvals, string fname);

tuple<vector<double>, vector<double>> read_1d_data(string fname);

void output_grid(vector<vector<double>> data, vector<double> xvals, vector<double> yvals, string fname);

tuple<vector<double>, vector<double>, vector<vector<double>>> read_grid(string fname);

double simple_integrate(vector<double> zs, vector<double> ys, double zmin=-1.0, double zmax=-1.0);

vector<double> get_col(vector<vector<double>> matrix, int col);

vector<vector<double>> transpose_matrix(vector<vector<double>> mat);

vector<double> get_chebyshev_spacing(double, double, int);

double patterson_in_regs(double a, double b, double epsrel, double epsabs,
                         double (*fptr)(double, void *p), void *p, int nregs, int mode);

double patterson_in_regs_parallel(double a, double b, double epsrel, double epsabs,
                                  double (*fptr)(double, void *p), void *p, int nregs, int mode);

double chebyshev_in_regs(double a, double b, int power,
                       double (*fptr)(double, void *p), void *p, int nregs, int mode);

double chebyshev_in_regs_parallel(double a, double b, int power,
                                double (*fptr)(double, void *p), void *p, int nregs, int mode);

// ========================================================================== //
// ============================== Interp Proxy ============================== //
// ========================================================================== //

class interp_1d
{
public:
    virtual ~interp_1d() {};

    virtual double evaluate(const double&) = 0;
    virtual double evaluate_deriv(const double&) = 0;

    virtual void update(const vector<double>&, const vector<double>&) = 0;
};

class gsl_loglog_interp_1d: public interp_1d
{
private:
    int mem_index;

public:
    gsl_loglog_interp_1d(): mem_index{set_trivial_spline("gsl_loglog_interp_1d")} {};
    gsl_loglog_interp_1d(const vector<double>&, const vector<double>&);
    gsl_loglog_interp_1d(const vector<double>&, function<double(double)>);
    gsl_loglog_interp_1d(double, double, int, function<double(double)>);
    ~gsl_loglog_interp_1d() { free_spline_JC(mem_index); }

    double evaluate(const double&);
    double evaluate_deriv(const double&);

    void update(const vector<double>&, const vector<double>&);
    void update(const vector<double>&, function<double(double)>);
    void update(double, double, int, function<double(double)>);
};

class gsl_loglin_interp_1d: public interp_1d
{
private:
    int mem_index;

public:
    gsl_loglin_interp_1d(): mem_index{set_trivial_spline("gsl_loglin_interp_1d")} {};
    gsl_loglin_interp_1d(const vector<double>&, const vector<double>&);
    gsl_loglin_interp_1d(const vector<double>&, function<double(double)>);
    gsl_loglin_interp_1d(double, double, int, function<double(double)>);
    ~gsl_loglin_interp_1d() { free_spline_JC(mem_index); }

    double evaluate(const double&);
    double evaluate_deriv(const double&);

    void update(const vector<double>&, const vector<double>&);
    void update(const vector<double>&, function<double(double)>);
    void update(double, double, int, function<double(double)>);

};

// ========================================================================== //
// =============================== Dict output ============================== //
// ========================================================================== //

class dict_creator
{
private:
    string separator;
    stringstream ss;
public:
    dict_creator(): separator{" : "}, ss{} {}
    dict_creator(const string& sep): separator{sep}, ss{} {}
    ~dict_creator() { ss.clear(); };

    void add_entry(const string& name, const double& value);
    void add_entry(const string& name, const string& value);
    void add_entry(const string& name, const vector<double>& value);
    void add_entry(const string& name, const dict_creator& dict);

    void write_file(const string& fname);
};

// ========================================================================== //
// =========================== 2D bspline Routine =========================== //
// ========================================================================== //

class bspline_2D
{
private:
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> _splines;
    double _xmin, _xstep, _ymin, _ystep;
    int _nx, _ny;
public:
    bspline_2D(vector<vector<double>>&, double, double, double, double);
    ~bspline_2D(){ _splines.clear(); }

    double operator()(double, double);

    boost::math::interpolators::cardinal_cubic_b_spline<double> get_spline_at_x(double);
    boost::math::interpolators::cardinal_cubic_b_spline<double> get_spline_at_y(double);
};

//======================================================================================
// omp functionality
//======================================================================================
int get_thread_num();
int get_max_threads();
void set_num_threads(int numthreads);

//======================================================================================
// searches of vectors etc
//======================================================================================
vector<double> sort_vector_ascending(vector<double> vec);
int get_index_of_min(vector<double> &vec);
int get_index_of_max(vector<double> &vec);
double get_min_value(vector<double> &vec);
double get_max_value(vector<double> &vec);

#endif
//======================================================================================
//======================================================================================
