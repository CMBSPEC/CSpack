#ifndef SS_GSTAR_H
#define SS_GSTAR_H

#include <vector>
#include <string>

using namespace std;

class SS_gstar
{
public:
    double unit_conversion;

    vector<double> a_coeffs{1.0, 1.11724, 3.12672e-1, -4.68049E-02, -2.65004E-02, -1.19760E-03, 1.82812E-04, 1.36436E-04, 8.55051E-05, 1.22840E-05, 3.82259E-07, -6.87035E-09};
    vector<double> b_coeffs{1.43382E-02, 1.37559E-02, 2.92108E-03, -5.38533E-04, -1.62496E-04, -2.87906E-05, -3.84278E-06, 2.78776E-06, 7.40342E-07, 1.17210E-07, 3.72499E-09, -6.74107E-11};
    vector<double> c_coeffs{1, 6.07869E-01, -1.54485E-01, -2.24034E-01, -2.82147E-02, 2.90620E-02, 6.86778E-03, -1.00005E-03, -1.69104E-04, 1.06301E-05, 1.69528E-06, -9.33311E-08};
    vector<double> d_coeffs{7.07388E+01, 9.18011E+01, 3.31892E+01, -1.39779E+00, -1.52558E+00, -1.97857E-02, -1.60146E-01, 8.22615E-05, 2.02651E-02, -1.82134E-05, 7.83943E-05, 7.13518E-05};
    vector<double> frho_coeffs{1.0, 1.03757, 0.508630, 0.0893988};
    vector<double> brho_coeffs{1.0, 1.03317, 0.398264, 0.0648056};
    vector<double> fs_coeffs{1.0, 1.03400, 0.456426, 0.0595248};
    vector<double> bs_coeffs{1.0, 1.03397, 0.342548, 0.0506182};
    vector<double> sfit_coeffs{1.0, 1.034, 0.456426, 0.0595249};

    double mass_e = 511.0e-6*1.0e+9;
    double mass_mu = 0.1056*1.0e+9;
    double mass_pi0 = 0.135*1.0e+9;
    double mass_pipm = 0.140*1.0e+9;
    double mass_1 = 0.5*1.0e+9;
    double mass_2 = 0.77*1.0e+9;
    double mass_3 = 1.2*1.0e+9;
    double mass_4 = 2.0*1.0e+9;

    SS_gstar(string unit="eV");

    ~SS_gstar();

    double get_gstar_rho(double temp);

    double get_gstar_s(double temp);

    double f_rho(double x);

    double b_rho(double x);

    double f_s(double x);

    double b_s(double x);

    double sfit(double x);
};

#endif
