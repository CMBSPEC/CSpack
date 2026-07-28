
#include "saikawa_shirai_gstar.h"
#include "routines.h"


SS_gstar::SS_gstar(string unit)
{
    if (unit == "eV"){unit_conversion = 1.0; return;}
    if (unit == "MeV"){unit_conversion = 1.0e-6; return;}
    if (unit == "GeV"){unit_conversion = 1.0e-9; return;}
    throw_error("SS_gstar::SS_gstar()", "Did not recognise unit " + unit, 1);
}

SS_gstar::~SS_gstar()
{
    a_coeffs.clear();
    b_coeffs.clear();
    c_coeffs.clear();
    d_coeffs.clear();
    frho_coeffs.clear();
    brho_coeffs.clear();
    fs_coeffs.clear();
    bs_coeffs.clear();
    sfit_coeffs.clear();
}


double SS_gstar::get_gstar_rho(double temp)
{
    temp *= unit_conversion;
    if (temp >= 120.0e6)
        return polyval(log(temp*1.0e-9), a_coeffs)/polyval(log(temp*1.0e-9), b_coeffs);

    return 2.030 + \
           1.353*pow(sfit(mass_e/temp), 4.0/3.0) + \
           3.495*f_rho(mass_e/temp) + \
           3.446*f_rho(mass_mu/temp) + \
           1.05*b_rho(mass_pi0/temp) + \
           2.08*b_rho(mass_pipm/temp) + \
           4.165*b_rho(mass_1/temp) + \
           30.55*b_rho(mass_2/temp) + \
           89.4*b_rho(mass_3/temp) + \
           8209.0*b_rho(mass_4/temp);
}

double SS_gstar::get_gstar_s(double temp)
{
    temp *= unit_conversion;
    if (temp >= 120.0e6)
        return get_gstar_rho(temp)/(1.0 + polyval(log(temp*1.0e+9), c_coeffs)/polyval(log(temp*1.0e+9), d_coeffs));

    return 2.008 + \
           1.923*sfit(mass_e/temp) + \
           3.442*f_s(mass_e/temp) + \
           3.468*f_s(mass_mu/temp) + \
           1.034*b_s(mass_pi0/temp) + \
           2.068*b_s(mass_pipm/temp) + \
           4.16*b_s(mass_1/temp) + \
           30.55*b_s(mass_2/temp) + \
           90.0*b_s(mass_3/temp) + \
           6209.0*b_s(mass_4/temp);
}

double SS_gstar::f_rho(double x)
{
    return exp(-1.04855*x)*polyval(x, frho_coeffs);
}

double SS_gstar::b_rho(double x)
{
    return exp(-1.03149*x)*polyval(x, brho_coeffs);
}

double SS_gstar::f_s(double x)
{
    return exp(-1.04190*x)*polyval(x, fs_coeffs);
}

double SS_gstar::b_s(double x)
{
    return exp(-1.03365*x)*polyval(x, bs_coeffs);
}

double SS_gstar::sfit(double x)
{
    return 1.0 + 7.0/4.0*exp(-1.0419*x)*polyval(x, sfit_coeffs);
}
