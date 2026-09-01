//==================================================================================================
// Created by Abir Sarkar on 15/11/2019 and modified by JC. These functions are based on
// Sarkar, Chluba and Lee, MNRAS, 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3705S/abstract)
//==================================================================================================

#include <gsl/gsl_sf_bessel.h>

#include "routines.h"
#include "Patterson.h"

#include "CSpack.h"

using namespace std;
using namespace CSpack_functions;
using namespace CSpack_kernels;

namespace CSpack_kernel_moments {

//==================================================================================================
//
// ANALYTICAL EXPRESSIONS OF THE MOMENTS [using v2 of arXiv:1905.00868]
//
//==================================================================================================
double compute_exact_moments_analytical(double omega0, double p0, int l)
{
    if(p0<=5.0e-2 && omega0 <=5.0e-2) return compute_nr_moments_analytical(omega0, p0, l);

    double moment = 0.0;
    double g0 = gamma_f(p0), o=omega0, o2=omega0*omega0;
    double ap=alphap(omega0, p0), am=alpham(omega0, p0);
    double lam=1.0 + 4.0*(g0+omega0)*omega0;

    // always need zeroth moment
    double fact1 =  ( 4.0*g0+9.0*o+2.0*g0*o2 )/(4.0*p0*o2) * log(ap/am);
    double fact2 = -( 1.0+1.0/lam-(1.0-2.0/o2)*log(lam) )/2.0;
    double fact3 =  f_moment(omega0, p0)/(p0*o);
    double moment0=(fact1+fact2+fact3)*3.0/(8.0*g0*o);

    if(l==0) moment=moment0;
    else if(l==1)
    {
        moment=(2.0*g0-o)/(2.0*o) * moment0;

        fact1 =(4.0*g0+(8.0-o2)*o)*log(lam);
        fact1-=(1.0+4.0*p0*p0+(5.0*g0+35.0/6.0*o+g0*o2)*o)/p0 * log(ap/am);
        fact1*=3.0/(32.0*g0*pow(o2, 2));
        fact2=-( 63.0+1.0/lam/lam+8.0*lam-(70.0-(6.0-4.0/lam)/lam)*o2 )/(64.0*g0*o*o2);

        moment+=fact1+fact2;
    }
    else if(l==2)
    {
        moment=(2.0+g0*(2.0*g0-o))/(2.0*o2) * moment0;

        fact1 =(2.0-(7.0-o2)*o2)/o * log(lam);
        fact1-=(155.0-90.0*o2)/(24.0*p0) * log(ap/am);
        fact1*=3.0/(16.0*g0*pow(o2, 2));
        fact2 =(4.0*g0+3.0*(2.0-o2)*o) * log(lam);
        fact2-=( 2.0*(6.0*g0+o)*g0 + (25.0+9.0*g0*o-6.0*o2)*o2 )/(3.0*p0) * log(ap/am);
        fact2*=3.0/(32.0*pow(o, 5));
        fact3 =( 162.0-(2.0-(5.0-1.0/lam)/lam)/lam -(141.0+23.0*lam)*lam )/(512.0*g0*pow(o, 5));
        fact3+=( 243.0+(12.0+1.0/(lam*lam))/lam+54.0*lam )/(64.0*g0*o*o2);
        fact3-=( 114.0-(3.0-(4.0+1.0/lam)/lam)/lam )/(32.0*g0*o);

        moment+=fact1+fact2+fact3;
    }

    return moment;
}

//==================================================================================================
// expressions assume omega0 && p0 << 1
//==================================================================================================
double compute_nr_moments_analytical(double omega0, double p0, int l)
{
    double moment = 0.0;

    if(l==0)
    {
        double fact11 = 1.0-(2.0-(26.0/5.0-(133.0/10.0-(1144.0/35.0-544.0/7.0*omega0)*omega0)*omega0)*omega0)*omega0;
        double fact12 = (5.0/3.0-(52.0/5.0-931.0/20.0*omega0)*omega0)*omega0;
        double fact13 = 7.0/12.0*omega0;
        moment = fact11 - fact12*pow(p0, 2) + fact13*pow(p0, 4);
    }

    else if(l==1)
    {
        double fact21 = -omega0*(1.0-(21.0/5.0-(147.0/10.0-(1616.0/35.0-940.0/7.0*omega0)*omega0)*omega0)*omega0);
        double fact22 = 4.0/3.0-(47.0/6.0-(189.0/5.0-9551.0/60.0*omega0)*omega0)*omega0;
        double fact23 = 553.0/120.0*omega0;
        moment = fact21 + fact22*pow(p0, 2) - fact23*pow(p0, 4);
    }

    else if(l==2)
    {
        double fact31 = pow(omega0, 2)*(7.0/5.0-(44.0/5.0-(1364.0/35.0-1020.0/7.0*omega0)*omega0)*omega0);
        double fact32 = 2.0/3.0-(42.0/5.0-(161.0/3.0-1886.0/7.0*omega0)*omega0)*omega0;
        double fact33 = 14.0/5.0-763.0/25.0*omega0;
        moment = fact31 + fact32*pow(p0, 2) + fact33*pow(p0, 4);
    }

    return moment;
}

//==================================================================================================
double compute_ur_moments_analytical(double omega0, double p0, int l)
{
    double moment = 0.0;
    double g1 = 4.0*p0*omega0;

    if(l==0)
    {
        double fact11 =(6.0*log(g1)-3.0)/(4.0*g1) - (42.239-(27.0-6.0*log(g1))*log(g1))/(2.0*pow(g1, 2));
        double fact12 =(39.0+24.0*log(g1))/(2.0*pow(g1, 3)) + (9.0 - 6.0*log(g1))*pow(omega0, 2)/(pow(g1, 3));

        moment = fact11+fact12;
    }

    else if(l==1)
    {
        moment = 7.0/(4.0*pow(g1,3)*pow(omega0,2))
                 + (18.0 - 36.0*log(g1) - pow(omega0,2)*(17.0 - 6.0*log(g1)))/pow(g1,3)
                 + (11.0 - 6.0*log(g1))/(4.0*g1) - (11.0 - 6.0*log(g1))/(16.0*pow(omega0,2))
                 + (12.0 + 9.0*log(g1))/(2.0*pow(g1*omega0,2))
                 + (94.739 - (54.0 - 6.0*log(g1))*log(g1))/(4.0*pow(g1,2))
                 - (26.870 - (18.0 - 3.0*log(g1))*log(g1))/(4.0*g1*pow(omega0,2));
    }

    else if(l==2)
    {
        double fact31 = (29.0-12.0*log(g1))/(8.0*g1) + (45.0-33.0*log(g1))/(4.0*pow(g1, 2))
                        + (529.0-300.0*log(g1))/(8.0*pow(g1, 3))
                        +(64.989-(43.0-6.0*log(g1))*log(g1))/(32.0*pow(omega0, 4));

        double fact32 = (109.0+96.0*log(g1))/(64.0*pow(omega0, 4)*g1) + 45.0/(64.0*pow(omega0, 4)*pow(g1, 2))
                        - 35.0/(192.0*pow(omega0, 4)*pow(g1, 3)) - ((29.0-12.0*log(g1))*g1)/(128.0*pow(omega0, 4));

        double fact33 = (65.0-24.0*log(g1))/(32.0*pow(omega0, 2))
                       + (214.74-(119.0-6.0*log(g1))*log(g1))/(16.0*pow(omega0, 2)*g1)
                       - (192.69+(21.0+36.0*log(g1))*log(g1))/(8.0*pow(omega0, 2)*pow(g1, 2))
                       + 23.0/(2.0*pow(omega0, 2)*pow(g1, 3)) + ((41.0-12.0*log(g1))*pow(omega0, 2))/(2.0*pow(g1, 3));

        moment = fact32 + fact33 - fact31;
    }

    return moment;

}

//==================================================================================================
double compute_recoil_moments_analytical(double omega0, double p0, int l)
{
    double moment = 0.0;
    double x = 1.0+2.0*omega0;

    if(l==0)
    {
        moment = 3.0*(1.0-x+15.0*pow(x,2)+pow(x,3))/(8.0*pow(1.0-x,2)*pow(x,2))
                +3.0*(3.0+6.0*x-pow(x,2))*log(x)/(4.0*pow(1.0-x,3));
    }
    else if(l==1)
    {
        moment =(2.0-5.0*x-3.0*pow(x,2)-71.0*pow(x,3)+5.0*pow(x,4))/(8.0*pow(1.0-x,2)*pow(x,3))
                -3.0*(7.0+6.0*x-pow(x,2))*log(x)/(4.0*pow(1.0-x,3));
    }
    else if(l==2)
    {
        moment = (3.0-11.0*x+12.0*pow(x,2)+28.0*pow(x,3)+177.0*pow(x,4)-17.0*pow(x,5))/(16.0*pow(1.0-x,2)*pow(x,4))
                +3.0*(11.0+6.0*x-pow(x,2))*log(x)/(4.0*pow(1.0-x,3));
    }

    return moment;
}


//==================================================================================================
double compute_doppler_moments_analytical(double omega0, double p0, int l)
{
    double moment = 0.0;

    if(l==0) moment = 1.0;
    else if(l==1) moment = 4.0/3.0*pow(p0, 2);
    else if(l==2) moment = 2.0/3.0*pow(p0, 2) + 14.0/5.0*pow(p0, 4);

    return moment;
}

//==================================================================================================
//
// explicit INTEGRATION OF THE MOMENTS
//
//==================================================================================================
double integrand_moment(double lp0, void *p)
{
    double *d=(double *)p;
    double p0=exp(lp0);
    return pow(p0, 3)*mb_dist_func(p0, d[2])*compute_exact_moments_analytical(d[0], p0, int(d[1]));
}
 
//==================================================================================================
double moment_Int_therm(double omega0, int l, double theta)
{
    double epsrel=1.0e-8, epsabs=1.0e-100;
    
    double d[3]={omega0, double(l), theta};

    double a=sqrt(2.0*theta)*1.0e-8;
    double f=30.0, b=sqrt((2.0+f*theta)*f*theta);
    double r=Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs, integrand_moment, &d);
    return r/mb_dist_norm(theta);
}
 
//==================================================================================================
double moment_Int_app1(double omega0, int l, double theta)
{
    double moment = 0.0;

    if(l==0)
    {
        double fact11 = 1.0 - 2.0*omega0 + (26.0*pow(omega0, 2))/5.0 - (133.0*pow(omega0, 3))/10.0 + (1144.0*pow(omega0, 4))/35.0 - 544.0/7.0*pow(omega0, 5);
        double fact12 = (5.0/3.0 - 52.0/5.0*omega0 +931.0/20.0*pow(omega0, 2))*omega0*(3.0*theta + (15.0*pow(theta, 2))/2.0 + (45.0*pow(theta, 3))/8.0 - (45.0*pow(theta, 4))/8.0);
        double fact13 = 7.0/12.0*omega0*(15.0*pow(theta, 2) + 90.0*pow(theta, 3) + 225.0*pow(theta, 4));
        moment = fact11 + fact12 + fact13;
    }

    else if(l==1)
    {
        double fact21 = omega0*(-1.0 + 21.0/5.0*omega0 - 147.0/10.0*pow(omega0, 2) + 1616.0/35.0*pow(omega0, 3) - 940.0/7.0*pow(omega0, 4));

        double fact22 =  (4.0/3.0 - 47.0/6.0*omega0 +189.0/5.0*pow(omega0, 2) - 9551.0/60.0*pow(omega0, 3) )*(3.0*theta + (15.0*pow(theta, 2))/2.0 + (45.0*pow(theta, 3))/8.0 - (45.0*pow(theta, 4))/8.0);

        double fact23 = 553.0/120.0*omega0*(15.0*pow(theta, 2) + 90.0*pow(theta, 3) + 225.0*pow(theta, 4));

        moment = fact21 + fact22 + fact23;
    }

    else if(l==2)
    {
        double fact31 = 1.0/5.0*pow(omega0, 2)*(7.0 - 44.0*omega0 + 1364.0/7.0*pow(omega0, 2));

        double fact32 = (2.0/3.0 - 42.0/5.0*omega0 +161.0/3.0*pow(omega0, 2) )*(3.0*theta + (15.0*pow(theta, 2))/2.0 + (45.0*pow(theta, 3))/8.0 - (45.0*pow(theta, 4))/8.0);

        double fact33 = 7.0/5.0*(2.0 - 109.0/5.0*omega0)*(15.0*pow(theta, 2) + 90.0*pow(theta, 3) + 225.0*pow(theta, 4));

        moment = fact31 + fact32 + fact33;
    }

    return moment;
}

double moment_Int_app2(double omega0, int l, double theta)
{
    double moment = 0.0;

    if(l==0)
    {
        double fact11 = 1.0 - 2.0*omega0 + (26.0*pow(omega0, 2))/5.0 - (133.0*pow(omega0, 3))/10.0 + (1144.0*pow(omega0, 4))/35.0 - 544.0/7.0*pow(omega0, 5);
        double fact12 = (5.0/3.0 - 52.0/5.0*omega0 +931.0/20.0*pow(omega0, 2))*omega0*(3.0*theta*gsl_sf_bessel_Kn(3, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta));
        double fact13 = 7.0/12.0*omega0*(15.0*pow(theta, 2)*gsl_sf_bessel_Kn(4, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta));
        moment = fact11 + fact12 + fact13;

    }

    else if(l==1)
    {
        double fact21 = omega0*(-1.0 + 21.0/5.0*omega0 - 147.0/10.0*pow(omega0, 2) + 1616.0/35.0*pow(omega0, 3) - 940.0/7.0*pow(omega0, 4));

        double fact22 =  (4.0/3.0 - 47.0/6.0*omega0 +189.0/5.0*pow(omega0, 2) - 9551.0/60.0*pow(omega0, 3) )*(3.0*theta*gsl_sf_bessel_Kn(3, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta));

        double fact23 = 553.0/120.0*omega0*(15.0*pow(theta, 2)*gsl_sf_bessel_Kn(4, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta));

        moment = fact21 + fact22 + fact23;
    }

    else if(l==2)
    {
        double fact31 = 1.0/5.0*pow(omega0, 2)*(7.0 - 44.0*omega0 + 1364.0/7.0*pow(omega0, 2));

        double fact32 = (2.0/3.0 - 42.0/5.0*omega0 +161.0/3.0*pow(omega0, 2) )*(3.0*theta*gsl_sf_bessel_Kn(3, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta));

        double fact33 = 7.0/5.0*(2.0 - 109.0/5.0*omega0)*(15.0*pow(theta, 2)*gsl_sf_bessel_Kn(4, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta));

        moment = fact31 + fact32 + fact33;
    }

    return moment;
}

double moment_Int_app3(double omega0, int l, double theta)
{
    double moment = 0.0;

    if(l==0)
    {
        double fact11 = ((6.0*omega0*(2.0 + omega0*(1.0+ omega0)*(8.0+ omega0)) + 3.0*pow(1.0+2.0*omega0, 2)*(-2.0 + (-2.0 + omega0)*omega0)*log(1.0+2.0*omega0))/(8.0*pow(omega0, 3)*pow(1.0+2.0*omega0, 2)));

        double fact12 = 1.0/(16.0*pow(omega0, 3)*pow(1.0+2.0*omega0, 4))*(2.0*omega0*(-6.0 + omega0*(-50.0 + omega0*(-157.0 + omega0*(-231.0 + 4.0*(-31.0 + omega0)*omega0)))) - pow(1.0+2.0*omega0, 4)*(-6.0 + omega0*(-8.0 + 3.0*omega0))*log(1.0+2.0*omega0));

        double fact13 = (1.0/(320.0*pow(omega0, 3)*pow(1.0+2.0*omega0, 6)))*(2.0*omega0*(90.0 + omega0*(1118.0 + omega0*(5803.0 +omega0*(16173.0 + 4.0*omega0*(6505.0 + 2.0*omega0*(2749.0 + 18.0*(50.0 - 3.0*omega0)*omega0)))))) + pow(1.0+2.0*omega0, 6)*(-90.0 + omega0*(-128.0 + 45.0*omega0))*log(1.0+2.0*omega0)) ;

        double fact14 = (2*omega0*(-1050 + omega0*(-17286 + omega0*(-124715 + omega0*(-515569 + 4*omega0*(-335073 + omega0*(-564583 + 4*omega0*(-147701 + 3*omega0*(-28637 - 6516*omega0 + 500*pow(omega0,2))))))))) -
                         3*pow(1 + 2*omega0,8)*(-350 + omega0*(-512 + 175*omega0))*log(1 + 2*omega0))/(4480.*pow(omega0,3)*pow(1 + 2*omega0,8));

        double fact15 = (2*omega0*(66150 + omega0*(1355154 + omega0*(12506901 + omega0*(68506187 + 4*omega0*(61722883 + 4*omega0*(38314734 + omega0*(66494503 + 2*omega0*(39458041 + 6*omega0*(5051887 + 4*omega0*(554448 + (93556 - 8575*omega0)*omega0)))))))))) + 3*pow(1 + 2*omega0,10)*(-22050 + omega0*(-32768 + 11025*omega0))*log(1 + 2*omega0))/(322560.*pow(omega0,3)*pow(1 + 2*omega0,10));

        moment = fact11 + fact12*(3.0*theta*gsl_sf_bessel_Kn(3, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta)) + fact13*(15.0*pow(theta, 2)*gsl_sf_bessel_Kn(4, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta))
        + fact14*(105.0*pow(theta, 3)*gsl_sf_bessel_Kn(5, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta))
        + fact15*(945.0*pow(theta, 4)*gsl_sf_bessel_Kn(6, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta));

    }

    else if(l==1)
    {
        double fact21 = 1.0/(8.0*pow(omega0, 3)*pow((1.0+2.0*omega0), 3))*(2.0*omega0*(-9.0 + omega0*(-51.0 + omega0*(-93.0 + omega0*(-51.0 + 10.0*omega0)))) - 3.0*(-3.0 + omega0)*(1.0 + omega0)*pow((1.0+2.0*omega0), 3)*log(1.0+2.0*omega0));

        double fact22 =  1.0/(48.0*pow(omega0, 3)*pow((1.0+2.0*omega0), 5))*(2.0*omega0*(39.0 + omega0*(363.0 + omega0*(1353.0 + omega0*(2429.0 + 2.0*omega0*(1033.0 + 2.0*(119.0 - 54.0*omega0)*omega0))))) + 3.0*pow((1.0+2.0*omega0), 5)*(-13.0 + omega0*(-4.0 + 3.0*omega0))*log(1.0+2.0*omega0));

        double fact23 = (1.0/(960.0*pow(omega0, 3)*pow((1.0+2.0*omega0), 7)))*(2.0*omega0*(-621.0 + omega0*(-8265.0 +  omega0*(-46659.0 + omega0*(-146219.0 + 2.0*omega0*(-133811.0 + 4.0*omega0*(-35475.0 + 2.0*omega0*(-9328.0 + omega0*(-803.0 + 882.0*omega0)))))))) -  3.0*pow((1.0+2.0*omega0), 7)*(-207.0 + omega0*(-64.0 + 45.0*omega0))*log(1.0+2.0*omega0)) ;

        double fact24 = 1.0/(4480.0*pow(omega0, 3)*(pow((1.0+2.0*omega0), 9)))*(2.0*omega0*(2475.0 + omega0*(42843.0 + omega0*(327681.0 + omega0*(1453225.0 +  2.0*omega0*(2065957.0 + 2.0*omega0*(1923657.0 + 2.0*omega0*(1165025.0 +  2.0*omega0*(429379.0 - 2.0*omega0*(-77103.0 + 2.0*omega0*(641.0 + 3550.0*omega0)))))))))) +3.0*pow((1.0 + 2.0*omega0),9)*(-825.0 + omega0*(-256.0 + 175.0*omega0))*log(1.0+2.0*omega0));

        double fact25 = 1.0/(322560.0*pow(omega0, 3)*(pow((1.0+2.0*omega0), 11)))*(2.0*omega0*(-158025.0 + omega0*(-3367677.0 +  omega0*(-32498767.0 + omega0*(-187382039.0 + 2.0*omega0*(-358705111.0 + 16.0*omega0*(-59864829.0 + omega0*(-112814093.0 +                                                           omega0*(-149156447.0 + 2.0*omega0*(-66797207.0 + 4.0*omega0*(-9238139.0 + omega0*(-2325491.0 +                                                                                                     omega0*(278941.0 + 227850.0*omega0)))))))))))) - 3.0*pow((1.0+2.0*omega0), 11)*(-52675.0 + omega0*(-16384.0 +  11025.0*omega0))*log(1.0+2.0*omega0));

        double fact26 = 1.0/(7096320.0*pow(omega0, 3)*pow((1.0+2.0*omega0), 13))*(2.0*omega0*(3155355.0 + omega0*(79866915.0 +omega0*(930559785.0 + omega0*(6606246545.0 +  2.0*omega0*(15938769653.0 + 2.0*omega0*(27602941431.0 + 2.0*omega0*(35275247881.0 + 4.0*omega0*(16758038141.0 + 2.0*omega0*(11783481765.0 + 2.0*omega0*(6004356693.0 +                                             2.0*omega0*(2111809098.0 + omega0*(913154373.0 - 2.0*omega0*(-78295993.0 + 20537338.0*omega0 +                                                                                                                 9128700.0*pow(omega0, 2)))))))))))))) + 15.0*pow((1.0+2.0*omega0), 13)*(-210357.0 + omega0*(-65536.0 + 43659.0*omega0))*log(1.0+2.0*omega0));

        double fact27 = 1.0/(123002880.0*pow(omega0, 3)*pow((1.0+2.0*omega0), 15))*(-2.0*omega0*(50426145.0 + omega0*(1478086845.0 +   omega0*(20179156575.0 + omega0*(170170774055.0 +  2.0*omega0*(495586237887.0 +  4.0*omega0*(527923872009.0 + 2.0*omega0*(849867503616.0 + omega0*(2103279963465.0 + 2.0*omega0*(2010759838845.0 +                                    8.0*omega0*(369996928128.0 + omega0*(413617080339.0 + omega0*(341012380617.0 +                                                             10.0*omega0*(19540335589.0 - 12.0*omega0*(-561883523.0 +2.0*omega0*(-29042882.0 +                                                                                                                                               21.0*omega0*(758075.0 +232078.0*omega0)))))))))))))))) + 15.0*pow((1.0+2.0*omega0), 15)*(3361743.0 + (1048576.0 - 693693.0*omega0)*omega0)*log(1.0+2.0*omega0)) ;



        moment = fact21 + fact22*(3.0*theta*gsl_sf_bessel_Kn(3, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta)) + fact23*(15.0*pow(theta, 2)*gsl_sf_bessel_Kn(4, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta)) + fact24*(105.0*pow(theta, 3)*gsl_sf_bessel_Kn(5, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta))+ fact25*(945.0*pow(theta, 4)*gsl_sf_bessel_Kn(6, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta))+0.0*fact26*(10395.0*pow(theta, 5)*gsl_sf_bessel_Kn(7, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta))+0.0*fact27*(135135.0*pow(theta, 6)*gsl_sf_bessel_Kn(8, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta));
    }

    else if(l==2)
    {
        double fact31 = (2*omega0*(12 + omega0*(90 + omega0*(247 + omega0*(283 + 92*omega0 - 34*pow(omega0,2))))) + 3*pow(1 + 2*omega0,4)*(-4 + (-2 + omega0)*omega0)*log(1 + 2*omega0))/(8.*pow(omega0,3)*pow(1 + 2*omega0,4));

        double fact32 = (2*omega0*(-12 + omega0*(-192 + omega0*(-1252 + omega0*(-4303 + omega0*(-8421 + 2*omega0*(-4510 + omega0*(-2167 + 60*omega0*(1 + 5*omega0)))))))) - 3*pow(1 + 2*omega0,6)*(-4 - 20*omega0 + 3*pow(omega0,3))*log(1 + 2*omega0))/
        (48.*pow(omega0,4)*pow(1 + 2*omega0,6));

        double fact33 =(2*omega0*(96 + omega0*(1716 + omega0*(13484 + omega0*(61379 + omega0*(175765 + 2*omega0*(163730 + omega0*(193501 + 8*omega0*(16555 - 2*omega0*(-2260 + omega0*(598 + 399*omega0)))))))))) + 3*pow(1 + 2*omega0,8)*(-32 - 92*omega0 + 15*pow(omega0,3))*log(1 + 2*omega0))/(320.*pow(omega0,4)*pow(1 + 2*omega0,8));

        double fact34 = (-2*omega0*(3840 + omega0*(82500 + omega0*(800780 + omega0*(4627465 + omega0*(17736627 + 2*omega0*(23440006 + omega0*(43475789 + 4*omega0*(13949295 + omega0*(11804815 + 4*omega0*(1435472 - 3*omega0*(-65649 + 4*omega0*(12381 + 4775*omega0)))))))))))) + 15*pow(1 + 2*omega0,10)*(256 + 636*omega0 - 105*pow(omega0,3))*log(1 + 2*omega0))/(13440.*pow(omega0,4)*pow(1 + 2*omega0,10));

        double fact35 = (2*omega0*(86016 + omega0*(2177868 + omega0*(25347028 + omega0*(179373341 + omega0*(860553267 + 2*omega0*(1478906198 + omega0*(3723868903 + 8*omega0*(866765493 + 4*omega0*(295888811 + omega0*(287960691 + omega0*(186563217 + 4*omega0*(16551110 + omega0*(193709 - 6*omega0*(369956 + 101675*omega0)))))))))))))) +  21*pow(1 + 2*omega0,12)*(-4096 - 9500*omega0 + 1575*pow(omega0,3))*log(1 + 2*omega0))/(322560.*pow(omega0,4)*pow(1 + 2*omega0,12));

        double fact36 = (2*omega0*(-589824 + omega0*(-17239428 + omega0*(-234450156 + omega0*(-1966810425 + omega0*(-11375616843 + 2*omega0*(-24016263678 + omega0*(-76456610433 + 4*omega0*(-46541597019 +                                                                                            omega0*(-87141699323 + 8*omega0*(-15596232756 + omega0*(-16768479985 + 8*omega0*(-1626839227 + omega0*(-833964401 + omega0*(-210853762 + omega0*(30546243 + 4*omega0*(9474881 + 2030805*omega0)))))))))))))))) - 9*pow(1 + 2*omega0,14)*(-65536 + 245*omega0*(-596 + 99*pow(omega0,2)))*log(1 + 2*omega0))/(2.36544e6*pow(omega0,4)*pow(1 + 2*omega0,14));

        double fact37 = ((2*omega0*(86507520 + omega0*(2869466820 + omega0*(44690457020 + omega0*(433871845495 + omega0*(2939505455601 +2*omega0*(7371620191386 + omega0*(28333571906145 - 8*omega0*(-10648491003843 +2*omega0*(-12645802083838 + omega0*(-23804572737506 + omega0*(-35375606677467 +                                                                                                                                                                                       8*omega0*(-5124164752725 + omega0*(-4515663102261 + 2*omega0*(-1440052715886 + 5*omega0*(-118378636427 + 8*omega0*(-2398153085 + 6*omega0*(171746135 + omega0*(110693206 + 19476765*omega0)))))))                                                                                                                                                                   ))))))))))))/pow(1 + 2*omega0,16) + 165*(-524288 + 1323*omega0*(-860 + 143*pow(omega0,2)))*log(1 + 2*omega0))/(3.6900864e8*pow(omega0,4));

        moment = fact31 + fact32*(3.0*theta*gsl_sf_bessel_Kn(3, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta)) + fact33*(15.0*pow(theta, 2)*gsl_sf_bessel_Kn(4, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta)) + fact34*(105.0*pow(theta, 3)*gsl_sf_bessel_Kn(5, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta)) + fact35*(945.0*pow(theta, 4)*gsl_sf_bessel_Kn(6, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta))+0.0*fact36*(10395.0*pow(theta, 5)*gsl_sf_bessel_Kn(7, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta))+0.0*fact37*(135135.0*pow(theta, 6)*gsl_sf_bessel_Kn(8, 1.0/theta))/(gsl_sf_bessel_Kn(2, 1.0/theta));
    }

    return moment;
}

//==================================================================================================
//
// compute moments numerically
//
//==================================================================================================
struct Integration_data_moments
{
    double omega0;   // omega0
    double theta;    // temperature
    int k;           // order of moment
    double (*thermal_kernel_ptr)(double, double, double);
    
    Integration_data_moments()
    {
        thermal_kernel_ptr=NULL;
        k=0;
    }
};

double integrand_moments_all(double lgomega, void *q)
{
    Integration_data_moments *d=(Integration_data_moments *)q;
    
    double omega=exp(lgomega);
    double K=d->thermal_kernel_ptr(d->omega0, omega, d->theta);
    double Dnu_nuk=pow(omega/d->omega0-1.0, d->k);
    
    return omega*K*Dnu_nuk;
}
 
//==================================================================================================
double moment_2D_Int_therm_all(double omega0, int k, double theta, string type)
{
    Integration_data_moments d;
    d.omega0=omega0;
    d.theta=theta;
    d.k=k;
    d.thermal_kernel_ptr=Get_kernel_pointer(type, "moment_2D_Int_therm_all");

    // rough estimate for width of averaged kernel
    double Dnu_nu2=sqrt(2.0*theta+7.0/5.0*omega0*omega0);
    
    double a=max(omega0*(1.0-30.0*Dnu_nu2), omega0*1.0e-8);
    double b=min(omega0*(1.0+30.0*Dnu_nu2), omega0*1.0e+8);

    double epsrel=1.0e-8, epsabs=1.0e-50;
    return Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs,
                                              integrand_moments_all, &d);
}

//==================================================================================================
//
// compute moments numerically
//
//==================================================================================================
struct Integration_2Ddata_moments
{
    double omega0, p0;      // omega0 and p0
    double theta;           // temperature
    int k;                  // order of moment
    kernel_ptr K;
    bool use_stim;

    Integration_2Ddata_moments()
    {
        K=NULL;
        k=0;
        use_stim=0;
    }
};

//==================================================================================================
double integrand_moments_all_o(double lgomega, void *q)
{
    Integration_2Ddata_moments *d=(Integration_2Ddata_moments *)q;
    
    double omega=exp(lgomega);
    double K=d->K(d->omega0, d->p0, omega);
    double Dnu_nuk=pow(omega/d->omega0-1.0, d->k);
    double stim=(d->use_stim ? one_minus_exp_mx(d->omega0/d->theta)/one_minus_exp_mx(omega/d->theta) : 1.0);

    return omega*K*Dnu_nuk*stim;
}
 
double integrand_moments_all_p(double lgp0, void *q)
{
    Integration_2Ddata_moments *d=(Integration_2Ddata_moments *)q;
    d->p0=exp(lgp0);
    
    double a=max(1.0e-16, omegamin(d->omega0, d->p0));
    double b=omegamax(d->omega0, d->p0);
    double fac=mb_dist_func(d->p0, d->theta) * pow(d->p0, 3);

    double epsrel=1.0e-9, epsabs=1.0e-100;
    return fac*Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs,
                                                  integrand_moments_all_o, &(*d));
}
 
//==================================================================================================
double moment_2D_Int_therm_all_II(double omega0, int k, double theta, string type, bool stim)
{
    Integration_2Ddata_moments d;
    d.omega0=omega0;
    d.theta=theta;
    d.k=k;
    d.use_stim=stim;
    d.K=Get_kernel_pointer(type, "moment_2D_Int_therm_all_II");

    double a=sqrt(2.0*theta)*1.0e-8;
    double f=30.0, b=sqrt((2.0+f*theta)*f*theta);

    double epsrel=1.0e-8, epsabs=1.0e-50;
    double r=Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs,
                                                integrand_moments_all_p, &d);
    return r/mb_dist_norm(theta);
}

//==================================================================================================
//
// compute moment 2Sigma2-Sigma1 numerically
//
//==================================================================================================
double G_integrand_moments_all_o(double lgomega, void *q)
{
    Integration_2Ddata_moments *d=(Integration_2Ddata_moments *)q;
    
    double omega=exp(lgomega);
    double K=d->K(d->omega0, d->p0, omega);
    double Dnu_nuk=(2.0*omega/d->omega0-3.0)*(omega/d->omega0-1.0);
    
    return omega*K*Dnu_nuk;
}
 
double G_integrand_moments_all_p(double lgp0, void *q)
{
    Integration_2Ddata_moments *d=(Integration_2Ddata_moments *)q;
    d->p0=exp(lgp0);
    
    double a=max(1.0e-16, omegamin(d->omega0, d->p0));
    double b=omegamax(d->omega0, d->p0);
    double fac=mb_dist_func(d->p0, d->theta) * pow(d->p0, 3);

    double epsrel=1.0e-8, epsabs=1.0e-50;
    return fac*Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs,
                                                  G_integrand_moments_all_o, &(*d));
}
 
//==================================================================================================
double G_moment_2D_Int_therm_all_II(double omega0, double theta, string type)
{
    Integration_2Ddata_moments d;
    d.omega0=omega0;
    d.theta=theta;
    d.K=Get_kernel_pointer(type, "G_moment_2D_Int_therm_all_II");

    double a=sqrt(2.0*theta)*1.0e-8;
    double f=30.0, b=sqrt((2.0+f*theta)*f*theta);

    double epsrel=1.0e-8, epsabs=1.0e-50;
    double r=Integrate_using_Patterson_adaptive(log(a), log(b), epsrel, epsabs,
                                                G_integrand_moments_all_p, &d);
    return r/mb_dist_norm(theta);
}

//==================================================================================================
//
// compute diffusion coefficients
//
//==================================================================================================
void compute_all_FP_coefficients(double omega0, double theta, string type,
                                 double &A, double &B, double &C)
{
    double domega0=omega0*1.0e-3;
    double x=omega0/theta;

    double Sigma1p1=moment_2D_Int_therm_all_II(omega0+domega0, 1, theta, type);
    double Sigma1=moment_2D_Int_therm_all_II(omega0, 1, theta, type);
    double Sigma1m1=moment_2D_Int_therm_all_II(omega0-domega0, 1, theta, type);

    double Sigma2p1=moment_2D_Int_therm_all_II(omega0+domega0, 2, theta, type);
    double Sigma2=moment_2D_Int_therm_all_II(omega0, 2, theta, type);
    double Sigma2m1=moment_2D_Int_therm_all_II(omega0-domega0, 2, theta, type);

    double dSigma1_domega0=(Sigma1p1-Sigma1m1)/(2.0*domega0);
    double dSigma2_domega0=(Sigma2p1-Sigma2m1)/(2.0*domega0);
    double d2Sigma2_domega02=(Sigma2p1-2.0*Sigma2+Sigma2m1)/(domega0*domega0);
    
    A=x*x*Sigma2/2.0/theta;
    B=x*(4.0*Sigma2-Sigma1+omega0*dSigma2_domega0)/theta;
    C=(3.0*(2.0*Sigma2-Sigma1)
          +omega0*(4.0*dSigma2_domega0-dSigma1_domega0)
          +omega0*omega0*d2Sigma2_domega02/2.0)/theta;

//    C=(3.0*(2.0*Sigma2-Sigma1))/theta;
//    C=(omega0*omega0*d2Sigma2_domega02)/theta;
//    C=omega0*(4.0*dSigma2_domega0-dSigma1_domega0)/theta;

    return;
}

}

namespace CSpack_kernel_moments_numerical {

struct Integration_1Ddata_moments
{
    double omega1{1}, p2{1};   // omega1 and p2
    int m{0};                  // order of moment
    double omin{1};
    kernel_ptr K{NULL};
    //bool use_stim{0};
};

double Integrand_moments_omega3(double lomega3, void *q)
{
    Integration_1Ddata_moments& d=*(Integration_1Ddata_moments *)q;

    double omega3=exp(lomega3);
    double K=d.K(d.omega1, d.p2, omega3);
    double Dnu_nuk=pow(omega3/d.omega1-1.0, d.m);
    
    return omega3*K*Dnu_nuk;
}

double Sigma_func(double omega1, double p2, int m, kernel_ptr K)
{
    Integration_1Ddata_moments D={omega1, p2, m, 0, K};
        
    double a=max(1.0e-16, omegamin(D.omega1, D.p2));
    double b=max(1.0e-16, omegamax(D.omega1, D.p2));
    double c=max(1.0e-16, omegacrit(D.omega1, D.p2));

    double epsrel=1.0e-8, epsabs=1.0e-50;
    double r=0.0;

    double la=log(a), lb=log(min(c, omega1));
    r+=Integrate_using_Patterson_adaptive(la, lb, epsrel, epsabs, Integrand_moments_omega3, &D);

    la=lb; lb=log(omega1);
    r+=Integrate_using_Patterson_adaptive(la, lb, epsrel, epsabs, Integrand_moments_omega3, &D);

    la=lb; lb=log(max(omega1, min(c, b)));
    r+=Integrate_using_Patterson_adaptive(la, lb, epsrel, epsabs, Integrand_moments_omega3, &D);

    la=lb; lb=log(b);
    r+=Integrate_using_Patterson_adaptive(la, lb, epsrel, epsabs, Integrand_moments_omega3, &D);

    return r;
}

}

//==================================================================================================
void output_kernel_moments(string fname, int np,
                           vector<double> omega0,
                           int m, kernel_ptr K)
{
    double p_l=1.0e-2, p_u=1.0e+3;
    vector<double> parr(np);
    init_xarr(p_l, p_u, &parr[0], np, 1, 0);

    ofstream ofile;
    ofile.open(fname.c_str());
    ofile.precision(10);

    ofile << "# omega1 = ";
    for(int io=0; io<(int)omega0.size(); io++) ofile << omega0[io] << " ";
    ofile << endl;
        
    for(int k=0; k<np; k++)
    {
        ofile << parr[k] << " ";
        
        for(int io=0; io<(int)omega0.size(); io++)
        {
            double norm=1.0;
            if(k>0) norm=CSpack_kernel_moments_numerical::Sigma_func(omega0[io], parr[k]*omega0[io], 0, K);
            
            ofile << CSpack_kernel_moments_numerical::Sigma_func(omega0[io], parr[k]*omega0[io], m, K)/norm << " ";
        }
        
        ofile << endl;
    }

    ofile.close();
}

//==================================================================================================
void output_kernel_moments(string fname, int np,
                           double omega0,
                           int m, kernel_ptr K)
{
    output_kernel_moments(fname, np, vector<double>{omega0}, m, K);
}

//==================================================================================================
//==================================================================================================
