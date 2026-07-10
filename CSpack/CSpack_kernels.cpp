//==================================================================================================
//  Created by Abir Sarkar on 15/11/2019 and modified by JC. These functions are based on
//  Sarkar, Chluba and Lee, MNRAS, 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#include "routines.h"
#include "Patterson.h"
#include "Compton_Kernel.h"

#include "CSpack.h"

using namespace std;
using namespace CSpack_functions;

//==================================================================================================
//
// COMPUTATION OF THE KERNEL
//
//==================================================================================================
double G_func(double omega, double omega0, double p0, double om0star, double omstar, double kappa)
{
    double common_fact = 3.0/(8.0*p0*gamma_f(p0)*pow(omega0,2));
    double lp=lambda_p(p0,omega0);
    double lm=lambda_m(p0,omega);

    double arg1 = pow(kappa/omstar ,2)*lp;
    double arg2 = pow(kappa/om0star,2)*lm;

    double fact1 = 2.0+pow((omstar-om0star)/(omega*omega0), 2)*(1.0+omega*omega0);
    double fact2 = sfunc(arg1)/ omstar     - sfunc(arg2)/ om0star;
    double fact3 = ffunc(arg1)/(omstar*lp) - ffunc(arg2)/(om0star*lm);

    return common_fact*kappa*(fact1+2.0*fact2+(1.0+omega*omega0)*fact3);
}

namespace CSpack_kernels {

//==================================================================================================
// EXACT KERNEL
//==================================================================================================
double kernel_exact(double omega0, double p0, double omega)
{
    double kerval=0.0;

    double omegab  = omegabar (omega0, p0, omega);
    double omegab0 = omega0bar(omega0, p0, omega);
    double p = pfunc_sc(omega0, p0, omega);
    double kappa1 = kappa(omega0, omega, p0, p);
    double kappa2 = kappa(omega, omega0, p, p0);

    double omega_min=omegamin(omega0,p0);
    double omega_cr=omegacrit(omega0,p0);
    double omega_max=omegamax(omega0,p0);

    if (omega0 > p0)
    {
        if(omega_min <= omega && omega < omega_cr)
            kerval = G_func(omega, omega0, p0, omegab0, omegab, kappa1);

        else if(omega_cr <= omega && omega <= omega0)
            kerval = G_func(omega, omega0, p0, omega, omega0, p0);

        else if(omega0 < omega && omega <= omega_max)
            kerval = G_func(omega, omega0, p0, omega0, omega, p);
    }
    else
    {
        if(omega_min <= omega && omega  < omega0)
            kerval = G_func(omega, omega0, p0, omegab0, omegab, kappa1);

        else if(omega0 <= omega && omega <=omega_cr)
            kerval = G_func(omega, omega0, p0, omegab, omegab0, kappa2);

        else if(omega_cr < omega && omega <= omega_max)
            kerval = G_func(omega, omega0, p0, omega0, omega, p);
    }

    return kerval;
}

//====================================================================================================================
// APPROXIMATE KERNELS
//====================================================================================================================
double kernel_recoil(double omega0, double p0, double omega)
{
    double delta = (omega0-omega)/omega;
    double kerval = 3.0/(8.0*pow(omega0,2))*(1.0 + pow(delta,2)/(1.0+delta) + pow((1.0-delta/omega0),2));
    return kerval;
}

double kernel_doppler(double omega0, double p0, double omega)
{
    double g0 = gamma_f(p0);
    double t = omega/omega0;
    double fact1 = (3.0+2.0*pow(p0, 2))/(2.0*p0)*(abs(log(t))-2.0*asinh(p0)) + (3.0+3.0*pow(p0, 2)+pow(p0, 4))/g0;
    double fact2 = 1.0+(10.0+8.0*pow(p0, 2)+4.0*pow(p0, 4))*t+pow(t,2);
    double kerval = 3.0/(8.0*omega0)*((1.0+t)/pow(p0, 5)*fact1 - abs(1.0-t)/(4.0*pow(p0, 6)*t)*fact2);
    return kerval;
}

double kernel_ur(double omega0, double p0, double omega)
{
    double g0 = gamma_f(p0);
    double G = 4.0*omega0*g0;
    double q = omega/(G*(g0-omega));
    double kerval = 3.0/(4.0*pow(g0,2)*omega0)*(2.0*q*log(q)+(1.0+2.0*q+pow(G*q, 2)/(2.0*(1.0 + G*q)))*(1.0-q));
    return kerval;
}

kernel_ptr Get_kernel_pointer(string type, string calling_func)
{
    if(type=="exact") return kernel_exact;
    else if(type=="recoil") return kernel_recoil;
    else if(type=="doppler") return kernel_doppler;
    else if(type=="ur") return kernel_ur;
    else throw_error(calling_func, "choose type 'exact', 'recoil', 'doppler', 'ur'", 1);

    return NULL;
}

double kernel_all(double omega0, double p0, double omega, string type)
{
    kernel_ptr K=Get_kernel_pointer(type, "kernel_all");
    return K(omega0, p0, omega);
}

//==================================================================================================
//
// COMPUTATION OF THE THERMALLY AVERAGED KERNEL
//
//==================================================================================================

//==================================================================================================
// LIMITS OF INTEGRATION
//==================================================================================================
double lower_limit(double omega0, double omega)
{
    double lim=0.0;
    if(omega <= omega0)
    {
        if(omega <= omega0/(1.0+2.0*omega0))
        {
            //lim = ( (omega0-omega)*sqrt((1.0+omega*omega0)/(omega*omega0)) - (omega0+omega) )/2.0;
            double N=2.0*( (omega0-omega)*sqrt((1.0+omega*omega0)/(omega*omega0))+(omega0+omega) );
            lim = (omega/omega0+omega0/omega-4.0*omega*omega0-2.0)/N;
        }

        else lim = 0.0;
    }

    else if(omega > omega0 && omega0 <= 0.5)
    {
        if(omega <= omega0/(1.0-2.0*omega0))
            lim = sqrt( (omega-omega0)*(omega-omega0+2.0) );

        else
        {
            //lim = ( (omega-omega0)*sqrt((1.0+omega*omega0)/(omega*omega0)) + (omega0+omega) )/2.0;
            double N=2.0*( (omega0-omega)*sqrt((1.0+omega*omega0)/(omega*omega0))+(omega0+omega) );
            lim = -(omega/omega0+omega0/omega-4.0*omega*omega0-2.0)/N;
        }
    }

    else if(omega > omega0 && omega0 > 0.5)
    {
        lim = sqrt( (omega-omega0)*(omega-omega0+2.0) );
    }

    return lim;
}

double lower_limit_dop(double omega0, double omega)
{
    double lim=0.0;
    // Modified limits for Doppler-dominated kernel integrals
    if(omega <= omega0) lim = (omega0 - omega)/(2.0*sqrt(omega*omega0));
    else lim = (omega - omega0)/(2.0*sqrt(omega*omega0));

    return lim;
}

//==================================================================================================
// FOR ALL OF THE ABOVE IN ONE MODULE
//==================================================================================================
struct Integration_data
{
    double omega0, omega;
    double theta;
    kernel_ptr K;
    
    Integration_data()
    {
        omega0=omega=theta=0.01;
        K=NULL;
    }
};

double integrand_p_all(double lp, void *q)
{
    Integration_data *d=(Integration_data *)q;
    double p=exp(lp);
    double fact= mb_dist_func(p, d->theta) * pow(p, 3) * d->K(d->omega0, p, d->omega);
    return fact;
}
 
double thermal_kernel_all(double omega0, double omega, double theta, string type)
{
    if(type=="SS_K") return PK_Kernel(omega0/const_h_mec2, omega/const_h_mec2, theta)/const_h_mec2;
    else if(type=="SS_C") return P_Compton(omega0/const_h_mec2, omega/const_h_mec2, theta)/const_h_mec2;

    double epsrel=1.0e-9, epsabs=1.0e-100;
    double pmax = sqrt( (theta*log(1.0e-30) - 2.0)*theta*log(1.0e-30) );
    double pb=pbar(theta);

    Integration_data d;
    d.omega=omega; d.omega0=omega0; d.theta=theta;
    d.K=Get_kernel_pointer(type, "thermal_kernel_all");

    double a, b;
    if (type=="doppler") a = lower_limit_dop(omega0, omega);
    else a = lower_limit(omega0, omega);

    a=max(pb*1.0e-12, a);
    b=max(pmax, pb*20.0);

    double r=Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs, integrand_p_all, &d);

//    double la=log(a), lb=min(log(a*4.0), log(b)), r=0.0;
//    r=Integrate_using_Patterson_adaptive(la, lb, epsrel, epsabs, integrand_p_all, &d);
//
//    la=lb; lb=min(log(a*8.0), log(b));
//    r+=Integrate_using_Patterson_adaptive(la, lb, epsrel, epsabs, integrand_p_all, &d);
//
//    la=lb; lb=min(log(a*16.0), log(b));
//    r+=Integrate_using_Patterson_adaptive(la, lb, epsrel, epsabs, integrand_p_all, &d);
//
//    la=lb; lb=log(b);
//    r+=Integrate_using_Patterson_adaptive(la, lb, epsrel, epsabs, integrand_p_all, &d);
//
    return r/mb_dist_norm(theta);
}
 
//==================================================================================================
double thermal_kernel_exact(double omega0, double omega, double theta)
{
    return thermal_kernel_all(omega0, omega, theta, "exact");
}
 
double thermal_kernel_recoil(double omega0, double omega, double theta)
{
    return thermal_kernel_all(omega0, omega, theta, "recoil");
}

double thermal_kernel_doppler(double omega0, double omega, double theta)
{
    return thermal_kernel_all(omega0, omega, theta, "doppler");
}

double thermal_kernel_ur(double omega0, double omega, double theta)
{
    return thermal_kernel_all(omega0, omega, theta, "ur");
}

//==================================================================================================
// kernels from Sazonov & Sunyaev 2000
//==================================================================================================
double thermal_kernel_SS_K(double omega0, double omega, double theta)
{
    return thermal_kernel_all(omega0, omega, theta, "SS_K");
}

double thermal_kernel_SS_C(double omega0, double omega, double theta)
{
    return thermal_kernel_all(omega0, omega, theta, "SS_C");
}

}

//==================================================================================================
// neutrino scattering KERNELS
//==================================================================================================
namespace CSpack_kernels_nu {

double gL =0.731;
double gR =gL-0.5;
double alpha_norm=gL*gL+gR*gR-gL*gR;
double alpha_LR=gL*gR/alpha_norm;
double beta_LR=gR*gR/alpha_norm;

//==================================================================================================
double I20_omega14(double pt, double lambda)
{
    return -pow(lambda, 5)/80.0+pow(lambda, 3)*pt*pt/24.0-lambda*pt*pt*pt*pt/16.0;
}

double I01_omega14(double pt, double lambda, double kappa, double Sigma)
{
    return pow(lambda, 3)/48.0-lambda*(pt*pt-kappa)/8.0-Sigma/16.0;
}

double I11_omega14(double pt, double lambda, double kappa, double Sigma)
{
    return -pow(lambda, 5)/160.0+pow(lambda, 3)*(3.0*pt*pt-2.0*kappa)/96.0
    -lambda*pt*pt*(pt*pt-kappa)/16.0-(pt*pt+lambda*lambda)*Sigma/32.0;
}

double I02_omega14(double pt, double lambda, double kappa, double Sigma, double omega1, double omega3)
{
    double Delta2=pow(omega1-omega3, 2), o1o3=omega1*omega3, pt2_kap=pt*pt-kappa;
    double special_term=(lambda==0.0 ? 0.0 : (Sigma*Sigma-4.0*pt*pt*pt*pt*Delta2)/lambda );
    
    return -3.0*pow(lambda, 5)/640.0+pow(lambda, 3)*(3.0*pt2_kap+Delta2)/96.0
    -lambda*(9.0*lambda*Sigma+24.0*(pt2_kap-o1o3)*o1o3 + 2.0*(5.0*pt*pt+3.0)*Delta2 )/64.0
    - 3.0*Sigma*pt2_kap/32.0 + special_term/128.0;
}

//==================================================================================================
double I_alpha(double pt, double lambda, double kappa, double Sigma)
{
    return I20_omega14(pt, lambda)-I01_omega14(pt, lambda, kappa, Sigma);
}

double I_beta(double pt, double lambda, double kappa, double Sigma, double omega1, double omega3)
{
    return I02_omega14(pt, lambda, kappa, Sigma, omega1, omega3)-2.0*I11_omega14(pt, lambda, kappa, Sigma);
}

//==================================================================================================
vector<double> Get_lambda12_lu(double omega1, double p2, double omega3, string zone)
{
    double rho=omega3/omega1;
    double p4=pfunc_sc(omega1, p2, omega3);
    
    if(zone=="I") return {1.0/(gamma_f(p2)+p2), rho*(gamma_f(p4)+p4)};
    
    else if(zone=="II")
    {
        double p2p=gamma_f(p2)+p2;
        if(p2<omega1) return {1.0/p2p, p2p};
        else if(p2>omega1) return {rho/(gamma_f(p4)+p4), p2p};
    }
    
    else if(zone=="III") return {rho/(gamma_f(p4)+p4), rho*(gamma_f(p4)+p4)};
    
    return {0, 0};
}

//==================================================================================================
vector<double> Get_lambda_lu(double omega1, double p2, double omega3, string zone)
{
    double p4=pfunc_sc(omega1, p2, omega3);
    
    if(zone=="I") return {p2+omega1, fabs(p4-omega3)};
    
    else if(zone=="II")
    {
        if(p2<omega1) return {p2+omega1, fabs(p2-omega1)};
        else if(p2>omega1) return {p4+omega3, fabs(p2-omega1)};
    }
    
    else if(zone=="III") return {p4+omega3, fabs(p4-omega3)};
    
    return {0, 0};
}

//==================================================================================================
vector<double> Get_Sigma_lu(double omega1, double p2, double omega3, string zone)
{
    // Sigma = (p2+omega1)*(p2-omega1)*(p4+omega3)*(p4-omega3)/lambda
    double p4=pfunc_sc(omega1, p2, omega3);
    double s=(p4>omega3 ? 1.0 : -1.0);
    
    // p4>omega3 ?
    if(zone=="I") return {(p2-omega1)*(p4+omega3)*(p4-omega3), (p2+omega1)*(p2-omega1)*(p4+omega3)*s};
    
    else if(zone=="II")
    {
        if(p2<omega1) return {(p2-omega1)*(p4+omega3)*(p4-omega3), -(p2+omega1)*(p4+omega3)*(p4-omega3)};
        else if(p2>omega1) return {(p2+omega1)*(p2-omega1)*(p4-omega3), (p2+omega1)*(p4+omega3)*(p4-omega3)};
    }
    
    else if(zone=="III") return {(p2+omega1)*(p2-omega1)*(p4-omega3), (p2+omega1)*(p2-omega1)*(p4+omega3)*s};
    
    return {0, 0};
}

//==================================================================================================
double kernel_exact(double omega1, double p2, double omega3)
{
    double g2=gamma_f(p2);
    double gt=g2+omega1, pt=pfunc(gt);
    double kappa=gt*(omega1+omega3)-2.0*omega1*omega3;
    
    string zone=Get_zone(omega1, p2, omega3);
    vector<double> lambda=Get_lambda_lu(omega1, p2, omega3, zone);
    vector<double> Sigma =Get_Sigma_lu (omega1, p2, omega3, zone);
    
    // comment: all functions scaled by omega1^4
    double I20 =I20_omega14(pt, lambda[1])
    -I20_omega14(pt, lambda[0]);
    
    double I_a =I_alpha(pt, lambda[1], kappa, Sigma[1])
    -I_alpha(pt, lambda[0], kappa, Sigma[0]);
    
    double I_b =I_beta (pt, lambda[1], kappa, Sigma[1], omega1, omega3)
    -I_beta (pt, lambda[0], kappa, Sigma[0], omega1, omega3);
    
    return (I20+alpha_LR*I_a+beta_LR*I_b)/pow(omega1, 4)/g2/p2;
}
}

//==================================================================================================
// thermally-averaged kernels over Fermi-Dirac distribution
//==================================================================================================
namespace CSpack_kernels_FD {

struct Integration_data
{
    double omega0, omega;
    double theta, mue, norm;
    kernel_ptr K;
    int add_FB{0};
};

//==================================================================================================
double f_FD(double p, double theta, double mue)
{
    double gamma=sqrt(1.0+p*p);
    double Dg=p*p/(1.0+gamma); // == gamma-1
    // [Comment: since we normalize to int p^2 f dp the possibly very
    //  small factor exp(-(1.0-mue)/theta) we analytically cancelled]
    return exp(-Dg/theta)/(1.0+exp( (mue-1.0-Dg)/theta ) );
}

double f_FD_blocking(double p, double theta, double mue)
{
    double gamma=sqrt(1.0+p*p);
    double em=exp(-(gamma-mue)/theta);
    return 1.0/(1.0+em);
}

double f_FD_blocking_rel(double x)
{
    double em=exp(-x);
    return 1.0/(1.0+em);
}

double f_Bose_stimul_rel(double x)
{
    double em=exp(-x);
    return 1.0/(1.0-em);
}

double integrand_norm_FD(double lp, void *q)
{
    Integration_data &d=*(Integration_data *)q;
    double p=exp(lp);
    return f_FD(p, d.theta, d.mue) * pow(p, 3);
}

double norm_FD(double theta, double mue)
{
    double epsrel=1.0e-9, epsabs=1.0e-50;
    double pb=pbar(theta);
    
    Integration_data d{0, 0, theta, mue, 1.0, NULL, 0};
    
    double a=max(pb*1.0e-12, 1.0e-16);
    double b=pb*1.0e+2;
    
    double r=Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs, integrand_norm_FD, &d);
    
    return r;
}

//==================================================================================================
double integrand_p_kernel(double lp, void *q)
{
    Integration_data &d=*(Integration_data *)q;
    double p=exp(lp);
    double fact= f_FD(p, d.theta, d.mue) * pow(p, 3) * d.K(d.omega0, p, d.omega);
    
    // blocking / stimulation factors for final state particles
    double fB_fac=(d.add_FB>0 ? f_FD_blocking(pfunc_sc(d.omega0, p, d.omega), d.theta, d.mue) : 1.0);
    fB_fac*=(d.add_FB==2 ? f_FD_blocking_rel(d.omega/d.theta) : 1.0); // neutrino final state blocking
    fB_fac*=(d.add_FB==3 ? f_Bose_stimul_rel(d.omega/d.theta) : 1.0); // photon final state stimulation
    
    return fact/d.norm * fB_fac;
}

double thermal_kernel_FD(double omega0, double omega, double theta, double mue,
                         kernel_ptr K, int add_FB)
{
    double epsrel=1.0e-9, epsabs=1.0e-50;
    double pmax = sqrt( (theta*log(1.0e-30) - 2.0)*theta*log(1.0e-30) );
    double pb=pbar(theta);
    double norm=norm_FD(theta, mue);
    
    Integration_data d{omega0, omega, theta, mue, norm, K, add_FB};
    
    double a=CSpack_kernels::lower_limit(omega0, omega), b;
    
    a=max(pb*1.0e-12, a);
    b=max(pmax, pb*20.0);
    
    double r=Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs, integrand_p_kernel, &d);
    
    return r;
}

//==================================================================================================
void output_thermal_kernel(string fname, int np,
                           vector<double> omega0,
                           double theta, double mue,
                           kernel_ptr K, int add_FB)
{
    double om_l=omega0[0]/theta*1.0e-3, om_u=omega0.back()/theta*1.0e+3;
    vector<double> oarr(np);
    init_xarr(om_l, om_u, &oarr[0], np, 1, 0);

    ofstream ofile;
    ofile.open(fname.c_str());
    ofile.precision(10);

    ofile << "# FD norm= " << norm_FD(theta, mue) << endl;
    ofile << "# omega1 = ";
    for(int io=0; io<(int)omega0.size(); io++) ofile << omega0[io] << " ";
    ofile << endl;
        
    for(int k=0; k<np; k++)
    {
        ofile << oarr[k]*theta << " ";
        for(int io=0; io<(int)omega0.size(); io++)
            ofile << thermal_kernel_FD(omega0[io], oarr[k]*theta, theta, mue, K, add_FB) << " ";
        
        ofile << endl;
    }

    ofile.close();
}

}

//==================================================================================================
void output_kernel(string fname, int np,
                   double omega0, vector<double> p0,
                   kernel_ptr K)
{
    double om_l=omegamin(omega0, p0.back()), om_u=omegamax(omega0, p0.back());
    vector<double> oarr(np);
    init_xarr(om_l, om_u, &oarr[0], np, 1, 0);

    ofstream ofile;
    ofile.open(fname.c_str());
    ofile.precision(10);

    ofile << "# p2 = ";
    for(int ip=0; ip<(int)p0.size(); ip++) ofile << p0[ip] << " ";
    ofile << endl;
        
    for(int k=0; k<np; k++)
    {
        ofile << oarr[k] << " ";
        for(int ip=0; ip<(int)p0.size(); ip++)
            ofile << K(omega0, p0[ip], oarr[k]) << " ";
        
        ofile << endl;
    }

    ofile.close();
}

//==================================================================================================
void output_kernel(string fname, int np,
                   double omega0, double p0,
                   kernel_ptr K)
{
    output_kernel(fname, np, omega0, vector<double>{p0}, K);
}

//==================================================================================================
//==================================================================================================
