//==================================================================================================
// ODE_solver routines
//
// Purpose: solve a system of stiff ODEs with moderate dimension (n<~100).
//
// Basic Aspects: The solver is based on a 6th order Gears method (implicit && stiffly-stable) with
// variable time-step. The linear algebra system is solved iteratively using a biconjugate gradient
// method. A guess for the next solution is extrapolated from the previous solutions;
//
// Author: Jens Chluba with contributions from Geoff Vasil at CITA
//
// First implementation: 12/01/2010
// Last modification   : 22/03/2021
//==================================================================================================
// 22/03/2021: moved LA-routines and ODE_solver_matrix to a separate file
// 01/03/2020: rework of setup parts and parameter communication [JC]
// 06/01/2020: significant improvement of Matrix routines and more general solver [JC]
// 31/03/2019: added Jacobian matrix setup check, which is useful for debugging [JC]
// 07/02/2019: tidied up solver structures using pointers [JC]
// 26/01/2019: added some step-size change settings + cleaned up the old routine [JC]
// 16/04/2018: added simple interpolation routine using the stored solutions
// 15/10/2013: added adaptive stepsize control to PDE+ODE+ADE routine
//
// 03/01/2012: modified PDE + ODE solver to allow additional OAEs.
//
// 28/06/2011: added option that allows to solve PDE+ODE systems, where the boundary conditions
//             are communicated as algebraic equations;
//
// 15/06/2011: set max order to 5; This was making the code a bit faster;
//             tidied the code a bit;
//--------------------------------------------------------------------------------------------------
// Solver details: https://ui.adsabs.harvard.edu/abs/2010MNRAS.407..599C/abstract
// Applications  : https://ui.adsabs.harvard.edu/abs/2011MNRAS.412..748C/abstract
//                 https://ui.adsabs.harvard.edu/abs/2012MNRAS.419.1294C/abstract
//==================================================================================================

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "ODE_solver.h"
#include "ODE_solver_LA.h"
#include "routines.h"

using namespace std;
using namespace ODE_solver_LA;

#ifdef OPENMP_ACTIVATED
//#define OPENMP_ACTIVATED_ODE
#endif

//=======================================================================================
namespace ODE_solver
{
    //===================================================================================
    //
    // Jacobian and ODE function part
    //
    //===================================================================================

    //===================================================================================
    // setup for indices of equations
    //===================================================================================
    void set_equation_indices(int nPDE, int nODE, int nOAE,
                              ODE_Solver_data &SD)
    {
        SD.nPDE=nPDE;
        SD.nODE=nODE;
        SD.nOAE=nOAE;
        SD.neq=nPDE+nODE+nOAE;

        SD.index_ODE=nPDE;
        SD.index_OAE=SD.index_ODE+nODE;

        return;
    }

    //===================================================================================
    // right hand side of the differential equation system dX/dz = g(z, X)
    //===================================================================================
    void ODE_solver_f(int neq, double z, double *y, double *f,
                      int col, const ODE_Solver_data &SD)
    {
        //==================================================================
        // evaluate the rhs of system  dy/dz == g(t,y)
        //==================================================================
        if(SD.p==NULL)
        {
           if(SD.fcn_col_ptr!=NULL) SD.fcn_col_ptr(&neq, &z, &y[0], &f[0], col);
           else throw_error("ODE_solver_f", "ODE-solver fcn not set", 1);
        }
        else
        {
            if(SD.fcn_col_para_ptr!=NULL) SD.fcn_col_para_ptr(&neq, &z, &y[0], &f[0], col, SD.p);
            else throw_error("ODE_solver_f", "ODE-solver fcn/parameters not set", 1);
        }

        return;
    }

    void ODE_solver_f(int neq, double z, double *y, double *f,
                      const ODE_Solver_data &SD)
    {
        ODE_solver_f(neq, z, y, f, -1, SD);
        return;
    }

    //===================================================================================
    //
    // convert to Jacobian
    //
    //===================================================================================
    void ODE_solver_convert_jac(int col, double Dz, double *r, ODE_Solver_data &SD)
    {
        // PDE + ODE part; for OAE part and boundary conditions nothing needs to be done
        for(int k=1; k<SD.index_ODE-1; k++) r[k]=-Dz*r[k];
        for(int k=SD.index_ODE; k<SD.index_OAE; k++) r[k]=-Dz*r[k];

        //===========================================================================
        // add unity to the i==j element for inner points of PDE and all ODEs
        //===========================================================================
        if( (col>0 && col<SD.index_ODE-1) ||
            (col>=SD.index_ODE && col<SD.index_OAE) ) r[col]+=1.0;

        return;
    }

    //===================================================================================
    //
    // to get the Jacobian from user
    //
    //===================================================================================
    void ODE_solver_jac_user(int neq, int col, double z, double Dz,
                             double *y, double *r, ODE_Solver_data &SD)
    {
        ODE_solver_f(neq, z, &y[0], &r[0], col, SD);

        //===============================================================================
        // define jacobian
        //===============================================================================
        ODE_solver_convert_jac(col, Dz, r, SD);

        return;
    }

    //===================================================================================
    //
    // to get numerical Jacobian of system
    //
    //===================================================================================
    inline double ODE_solver_df_dx_2point(const double &fp1, const double &fm1,
                                          const double &h, double eps)
    {
        if(fm1==0.0) return (fp1-fm1)/(2.0*h);
        double dum=fp1/fm1-1.0;
        return fabs(dum)<=eps ? 0.0 : dum*fm1/(2.0*h);
    }

    inline double ODE_solver_df_dx_4point(const double &fp2, const double &fp1,
                                          const double &fm1, const double &fm2,
                                          const double &h)
    { return (8.0*(fp1-fm1)+fm2-fp2)/(12.0*h); }

    //===================================================================================
    void ODE_solver_jac_2p(int neq, int col, double z, double Dz,
                           double *y, double *r, ODE_Solver_data &SD)
    {
        //===============================================================================
        // r[i] is reset here
        // y[i] should contain the solution at z
        //===============================================================================
        double y0=y[col], Dyj=y[col]*1.0e-5, eps=1.0e-14;
        if(y0==0.0) Dyj=1.0e-14;

        vector<double> fp1(neq), f0(neq), fm1(neq);

        //===============================================================================
        // derivatives with respect to Xj.
        // A two-point formula is used. It has accuracy O(h^2)
        //===============================================================================
        // get f(yj+Dyj)
        y[col]=y0+Dyj;
        ODE_solver_f(neq, z, &y[0], &fp1[0], SD);
        // get f(yj-Dyj)
        y[col]=y0-Dyj;
        ODE_solver_f(neq, z, &y[0], &fm1[0], SD);
        // restore y again
        y[col]=y0;
        ODE_solver_f(neq, z, &y[0], &f0[0], SD); // required for consistency

        //===============================================================================
        // define numerical derivative
        //===============================================================================
#ifdef OPENMP_ACTIVATED_ODE
#pragma omp parallel for default(shared)
#endif
        for(int k=0; k<neq; k++) r[k]=ODE_solver_df_dx_2point(fp1[k], fm1[k], Dyj, eps);

        //===============================================================================
        // define jacobian
        //===============================================================================
        ODE_solver_convert_jac(col, Dz, r, SD);

        return;
    }

    //===================================================================================
    void ODE_solver_jac_4p(int neq, int col, double z, double Dz,
                           double *y, double *r, ODE_Solver_data &SD)
    {
        //===============================================================================
        // r[i] is reset here
        // y[i] should contain the solution at z
        //===============================================================================
        double y0=y[col], Dyj=y[col]*1.0e-5;
        if(y0==0.0) Dyj=1.0e-14;

        vector<double> fp1(neq), fp2(neq), f0(neq), fm1(neq), fm2(neq);

        //===============================================================================
        // derivatives with respect to Xj.
        // A four-point formula is used. It has accuracy O(h^4)
        //===============================================================================
        // get f(yj+2Dyj)
        y[col]=y0+2*Dyj;
        ODE_solver_f(neq, z, &y[0], &fp2[0], SD);
        // get f(yj+Dyj)
        y[col]=y0+Dyj;
        ODE_solver_f(neq, z, &y[0], &fp1[0], SD);
        // get f(yj-Dyj)
        y[col]=y0-Dyj;
        ODE_solver_f(neq, z, &y[0], &fm1[0], SD);
        // get f(yj-2Dyj)
        y[col]=y0-2*Dyj;
        ODE_solver_f(neq, z, &y[0], &fm2[0], SD);
        // restore y again
        y[col]=y0;
        ODE_solver_f(neq, z, &y[0], &f0[0], SD); // required for consistency

        //===============================================================================
        // define numerical derivative
        //===============================================================================
#ifdef OPENMP_ACTIVATED_ODE
#pragma omp parallel for default(shared)
#endif
        for(int k=0; k<neq; k++) r[k]=ODE_solver_df_dx_4point(fp2[k], fp1[k],
                                                              fm1[k], fm2[k], Dyj);

        //===============================================================================
        // define jacobian
        //===============================================================================
        ODE_solver_convert_jac(col, Dz, r, SD);

        return;
    }

    //===================================================================================
    void ODE_solver_compute_Jacobian_Matrix(double z, double Dz, int neq, double *y,
                                            ODE_solver_matrix &Jac,
                                            ODE_Solver_data &SD, int mess=0)
    {
        //===============================================================================
        // z is the current redshift
        // y[] should contain the solution for which the jacobian is needed
        // Jac[] should have dimension neq*neq
        //===============================================================================
        if(mess==1) cout << " entering full Jacobinan update. " << endl;

        Jac.clear();

        //===============================================================================
        // this vector will be used to store the jacobian after calls of f_i(X)
        //===============================================================================
        vector<double> Jij(neq);

        //===============================================================================
        // fill Matrix with values
        //===============================================================================
        void (*Jac_func)(int, int, double, double, double *, double *, ODE_Solver_data &);

        if(SD.num_jac) Jac_func=ODE_solver_jac_2p;
        else Jac_func=ODE_solver_jac_user;

        for(int J=0; J<neq; J++)
        {
            Jac_func(neq, J, z, Dz, y, &Jij[0], SD);

            for(int row=0; row<neq; row++)
                if(Jij[row]!=0.0) Jac.save_info(neq, J, row, Jij[row]);
        }

        if(mess==1) cout << " Number of non-zero elements = " << Jac.A.size()
                         << ". Full matrix has " << neq*neq << " elements " << endl;

        if((int)Jac.diags.size()!=neq)
            cout << " Hmmm. That should not happen... " << endl;

        return;
    }

    //===================================================================================
    //
    // this function has to be called to set/reset errors
    //
    //===================================================================================
    void ODE_Solver_set_errors(const ODE_solver_accuracies &tolsv,
                               ODE_Solver_data &ODE_Solver_info)
    {
        //===============================================================================
        // Copy absolute & relative accuracies
        //===============================================================================
        ODE_Solver_info.tols=tolsv;
        if(ODE_Solver_info.tols.neq==0)
            throw_error("ODE_Solver_set_errors", "check tolerance array dimensions", 1);

        // find minimal relative tolerance
        double minrtol=1.0;
        for(int k=0; k<ODE_Solver_info.tols.neq; k++)
            if(minrtol>ODE_Solver_info.tols.rel[k] && ODE_Solver_info.tols.rel[k]>0.0)
                minrtol=ODE_Solver_info.tols.rel[k];

        if(ODE_Solver_info.verbosity>=1)
            cout << " setting tolerance for Jacobian iteration to " << minrtol << endl;
        ODE_Solver_info.tolSol=minrtol;
    }

    //===================================================================================
    //
    // memory setup for ODE_solver_Solution and ODE_solver_accuracies
    //
    //===================================================================================
    void ODE_Solver_data::clear()
    {
        Stemp.clear();
        Snewptr->clear();
        Sguessptr->clear();
        for(int k=0; k<6; k++) Snptr[k]->clear();

        F.clear();
        dY.clear();
        tols.clear();
        Jacobian_ptr->clear();

        delete Snewptr;
        delete Sguessptr;
        for(int k=0; k<6; k++) delete Snptr[k];
        delete Jacobian_ptr;
    }

    //===================================================================================
    //
    // this function has to be called to set up the memory and the initial solution
    //
    //===================================================================================
    void ODE_Solver_setup(ODE_solver_Solution &Sz,
                          ODE_solver_accuracies &tols,
                          ODE_Solver_data &ODE_Solver_info)
    {
        int neq=Sz.y.size();

        if(ODE_Solver_info.verbosity>=1)
        {
            if(ODE_Solver_info.Snewptr==NULL) cout << "\n Setting up memory " << endl;
            else  cout << "\n Resetting memory " << endl;
            cout << " Total number of equations = " << neq << endl;
        }

        //-------------------------------------------------------------
        // create vectors & allocate memory
        //-------------------------------------------------------------
        ODE_Solver_info.F.resize(neq);
        ODE_Solver_info.dY.resize(neq);
        ODE_Solver_info.Fcorr1.resize(neq);
        ODE_Solver_info.Fcorr2.resize(neq);
        ODE_Solver_info.Sptr=NULL;

        if(ODE_Solver_info.Snewptr==NULL)
        {
            //wait_f_r("Memory");
            ODE_Solver_info.Snewptr=new ODE_solver_Solution;
            ODE_Solver_info.Sguessptr=new ODE_solver_Solution;
            for(int k=0; k<6; k++) ODE_Solver_info.Snptr[k]=new ODE_solver_Solution;
            ODE_Solver_info.Jacobian_ptr= new ODE_solver_matrix;
        }

        ODE_Solver_info.Stemp.y.resize(neq); ODE_Solver_info.Stemp.dy.resize(neq);
        //
        ODE_Solver_info.Snewptr->y.resize(neq); ODE_Solver_info.Snewptr->dy.resize(neq);
        ODE_Solver_info.Sguessptr->y.resize(neq); ODE_Solver_info.Sguessptr->dy.resize(neq);
        //
        for(int k=0; k<6; k++)
        {
            ODE_Solver_info.Snptr[k]->y.resize(neq);
            ODE_Solver_info.Snptr[k]->dy.resize(neq);
        }

        //===============================================================================
        // copy solution
        //===============================================================================
        ODE_Solver_info.Snptr[0]->z=Sz.z;
        ODE_Solver_info.Snptr[0]->y=Sz.y;
        ODE_Solver_info.Snptr[0]->dy=Sz.dy;

        ODE_solver_f(neq, ODE_Solver_info.Snptr[0]->z,
                     &ODE_Solver_info.Snptr[0]->y[0],
                     &ODE_Solver_info.Snptr[0]->dy[0],
                     ODE_Solver_info);

        //===============================================================================
        // other data
        //===============================================================================
        ODE_Solver_info.Dz_z_last=0;
        ODE_Solver_info.zstart_solver=ODE_Solver_info.Snptr[0]->z;
        ODE_Solver_info.order=1;
        ODE_Solver_info.count=1;
        ODE_Solver_info.Jac_is_set=0;
        //
        ODE_Solver_info.n_up=0;
        ODE_Solver_info.n_down=0;

        //===============================================================================
        // Copy absolute & relative accuracies
        //===============================================================================
        ODE_Solver_set_errors(tols, ODE_Solver_info);

        return;
    }

    //===================================================================================
    void ODE_Solver_set_up_solution_and_memory(ODE_solver_Solution &Sz,
                                               ODE_solver_accuracies &tols,
                                               ODE_Solver_data &ODE_Solver_info,
                                               void (*fcn_ptr)(int *neq, double *z,
                                                               double *y, double *f,
                                                               int col),
                                               bool num_jac)
    {
        ODE_Solver_info.fcn_col_ptr=fcn_ptr;
        ODE_Solver_info.num_jac=num_jac;
        ODE_Solver_setup(Sz, tols, ODE_Solver_info);

        return;
    }

    void ODE_Solver_set_up_solution_and_memory(ODE_solver_Solution &Sz,
                                               ODE_solver_accuracies &tols,
                                               ODE_Solver_data &ODE_Solver_info,
                                               void (*fcn_ptr)(int *neq, double *z,
                                                               double *y, double *f,
                                                               int col, void *p),
                                               bool num_jac)
    {
        if(ODE_Solver_info.p==NULL) throw_error("ODE_Solver_set_up_solution_and_memory",
                                                "set initial parameter pointer *p", 1);

        ODE_Solver_info.fcn_col_para_ptr=fcn_ptr;
        ODE_Solver_info.num_jac=num_jac;
        ODE_Solver_setup(Sz, tols, ODE_Solver_info);

        return;
    }

    //===================================================================================
    void ODE_Solver_set_up_solution_and_memory_Anisotropies(ODE_solver_Solution &Sz,
                                                            ODE_solver_accuracies &tols,
                                                            ODE_Solver_data &ODE_Solver_info,
                                                            void (*fcn_ptr)(int *neq, double *z,
                                                                            double *y, double *f,
                                                                            int col, void *p),
                                                            bool num_jac)
    {
        ODE_Solver_set_up_solution_and_memory(Sz, tols, ODE_Solver_info, fcn_ptr, num_jac);
        return;
    }

    //===================================================================================
    //
    // time-step using Gears-method
    // order can be ==1; 2; 3; 4; 5; 6;
    //
    //===================================================================================

    //===================================================================================
    // aux-functions alpha_i coefficients
    //===================================================================================
    double ODE_Solver_Gears_fk1(double rk)
    { return 2.0+rk; }

    double ODE_Solver_Gears_fk2(double r1, double rk)
    {
        double r=ODE_Solver_Gears_fk1(r1)*ODE_Solver_Gears_fk1(rk);
        return r-1.0;
    }

    double ODE_Solver_Gears_fk3(double r1, double r2, double rk)
    {
        double r=ODE_Solver_Gears_fk2(r1, r2)*ODE_Solver_Gears_fk1(rk);
        return r-ODE_Solver_Gears_fk1(r1+r2);
    }

    double ODE_Solver_Gears_fk4(double r1, double r2, double r3, double rk)
    {
        double r=ODE_Solver_Gears_fk3(r1, r2, r3)*ODE_Solver_Gears_fk1(rk);
        r-=ODE_Solver_Gears_fk1(r1+r2)*r3;
        return r-ODE_Solver_Gears_fk2(r1, r2);
    }

    inline double ODE_Solver_Gears_fk5(double r1, double r2, double r3,
                                       double r4, double rk)
    {
        double r=ODE_Solver_Gears_fk4(r1, r2, r3, r4)*ODE_Solver_Gears_fk1(rk);
        r-=(ODE_Solver_Gears_fk2(r1, r2+r3)+r2*r3)*r4;
        return r-ODE_Solver_Gears_fk3(r1, r2, r3);
    }

    //===================================================================================
    // alpha_i coefficients
    //===================================================================================
    double ODE_Solver_delta0(double r1, double r2, double r3, double r4, double r5,
                             double a1, double a2, double a3, double a4, double a5)
    { return 1.0 + a1*r1 + a2*r2 + a3*r3 + a4*r4 + a5*r5; }

    double ODE_Solver_alp0(double a1, double a2, double a3, double a4, double a5)
    { return 1.0 - a1 - a2 - a3 - a4 - a5; }

    double ODE_Solver_alp1(double r1, double r2, double r3, double r4, double r5,
                           double a2, double a3, double a4, double a5)
    { return -((1.0 + a2*r2*(2.0+r2) + a3*r3*(2.0+r3) + a4*r4*(2.0+r4)
                    + a5*r5*(2.0+r5))/(r1*(2.0+r1))); }

    double ODE_Solver_alp2(double r1, double r2, double r3, double r4, double r5,
                           double a3, double a4, double a5)
    {
        double t=(1.0+r1);
        return -( t*t
                 +a3*r3*(r1-r3)*ODE_Solver_Gears_fk2(r1, r3)
                 +a4*r4*(r1-r4)*ODE_Solver_Gears_fk2(r1, r4)
                 +a5*r5*(r1-r5)*ODE_Solver_Gears_fk2(r1, r5)
                )
                /( r2*(r1-r2)*ODE_Solver_Gears_fk2(r1, r2) );
    }

    double ODE_Solver_alp3(double r1, double r2, double r3, double r4, double r5,
                           double a4, double a5)
    {
        double t=(1.0+r1)*(1.0+r2);
        return -( t*t
                 +a4*r4*(r1-r4)*(r2-r4)*ODE_Solver_Gears_fk3(r1, r2, r4)
                 +a5*r5*(r1-r5)*(r2-r5)*ODE_Solver_Gears_fk3(r1, r2, r5)
                )
                /( r3*(r1-r3)*(r2-r3)*ODE_Solver_Gears_fk3(r1, r2, r3) ) ;
    }

    double ODE_Solver_alp4(double r1, double r2, double r3, double r4, double r5,
                           double a5)
    {
        double t=(1.0+r1)*(1.0+r2)*(1.0+r3);
        return -( t*t
                 +a5*r5*(r1-r5)*(r2-r5)*(r3-r5)*ODE_Solver_Gears_fk4(r1, r2, r3, r5)
                )
                /( r4*(r1-r4)*(r2-r4)*(r3-r4)*ODE_Solver_Gears_fk4(r1, r2, r3, r4) );
    }

    double ODE_Solver_alp5(double r1, double r2, double r3, double r4, double r5)
    {
        double t=(1.0+r1)*(1.0+r2)*(1.0+r3)*(1.0+r4);
        return -t*t/( r5*(r1-r5)*(r2-r5)*(r3-r5)*(r4-r5)
                      *ODE_Solver_Gears_fk5(r1, r2, r3, r4, r5) );
    }

    double ODE_Solver_alpi(int i, const double *r, const double *a)
    {
        if(i==2) return ODE_Solver_alp1(r[0], r[1], r[2], r[3], r[4], a[2], a[3], a[4], a[5]);
        if(i==3) return ODE_Solver_alp2(r[0], r[1], r[2], r[3], r[4], a[3], a[4], a[5]);
        if(i==4) return ODE_Solver_alp3(r[0], r[1], r[2], r[3], r[4], a[4], a[5]);
        if(i==5) return ODE_Solver_alp4(r[0], r[1], r[2], r[3], r[4], a[5]);
        if(i==6) return ODE_Solver_alp5(r[0], r[1], r[2], r[3], r[4]);

        return 0;
    }

    //===================================================================================
    // Gear's corrector (implicit)
    //===================================================================================
    double ODE_Solver_compute_ynp1(int order, ODE_solver_Solution &Snp1_new,
                                   const ODE_solver_Solution &Snp1_guess,
                                   ODE_solver_Solution *Sn[6])
    {
        if(order<1 || order>6)
        {
            cerr << " check order for ODE_Solver_compute_ynp1 " << endl;
            exit(0);
        }

        //===============================================================================
        // the structures yn, ynm1, ynm2 contain the information from the previous time
        // steps the structure ynp1 contains the current version of ynp1 and dynm1
        //===============================================================================
        double delta=1.0;
        double Dznp1=Snp1_guess.z-Sn[0]->z;
        //
        double a[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double Dzn[5]={0.0, 0.0, 0.0, 0.0, 0.0};
        double rho[5]={0.0, 0.0, 0.0, 0.0, 0.0};

        for(int k=2; k<=order; k++){ Dzn[k-2]=Sn[0]->z-Sn[k-1]->z; rho[k-2]=Dzn[k-2]/Dznp1; }

        for(int k=order; k>=2; k--) a[k-1]=ODE_Solver_alpi(k, rho, a);

        a[0]=ODE_Solver_alp0(a[1], a[2], a[3], a[4], a[5]);

        delta=ODE_Solver_delta0(rho[0], rho[1], rho[2], rho[3], rho[4],
                                a[1], a[2], a[3], a[4], a[5]);

        // add up different orders
        for(int i=0; i<(int)Snp1_new.y.size(); i++)
        {
            Snp1_new.y[i]=Dznp1*delta*Snp1_guess.dy[i];

            double Sa=0.0;
            for(int mo=order-1; mo>=0; mo--) Sa+=a[mo]*Sn[mo]->y[i];

            Snp1_new.y[i]+=Sa;
        }

        return delta;
    }

    //===================================================================================
    // beta_i coefficients for extrapolation with avariable step-size
    // These expression where derived by Geoff Vasil (CITA)
    //===================================================================================
    double ODE_Solver_beta0(double b1, double b2, double b3, double b4, double b5)
    { return 1.0 - b1 - b2 - b3 - b4 - b5; }

    double ODE_Solver_beta1(double r1, double r2, double r3, double r4, double r5,
                            double b2, double b3, double b4, double b5)
    { return -(1.0 + b2*r2 + b3*r3 + b4*r4 + b5*r5)/r1; }

    double ODE_Solver_beta2(double r1, double r2, double r3, double r4, double r5,
                            double b3, double b4, double b5)
    { return -( (1.0+r1) + b3*r3*(r1-r3) + b4*r4*(r1-r4) + b5*r5*(r1-r5) )/((r1-r2)*r2); }

    double ODE_Solver_beta3(double r1, double r2, double r3, double r4, double r5,
                            double b4, double b5)
    { return ( (1.0+r1)*(1.0+r2) + b4*r4*(r1-r4)*(r2-r4)
                                 + b5*r5*(r1-r5)*(r2-r5) )/((r1-r3)*r3*(r3-r2)); }

    double ODE_Solver_beta4(double r1, double r2, double r3, double r4, double r5,
                            double b5)
    { return -( (1.0+r1)*(1.0+r2)*(1.0+r3) + b5*r5*(r1-r5)*(r2-r5)*(r3-r5) )
              /((r3-r4)*r4*(r4-r1)*(r4-r2)); }

    double ODE_Solver_beta5(double r1, double r2, double r3, double r4, double r5)
    { return (1.0+r1)*(1.0+r2)*(1.0+r3)*(1.0+r4)/((r2-r5)*r5*(r5-r1)*(r5-r3)*(r5-r4)); }

    double ODE_Solver_betai(int i, const double *r, const double *b)
    {
        if(i==2) return ODE_Solver_beta1(r[0], r[1], r[2], r[3], r[4], b[2], b[3], b[4], b[5]);
        if(i==3) return ODE_Solver_beta2(r[0], r[1], r[2], r[3], r[4], b[3], b[4], b[5]);
        if(i==4) return ODE_Solver_beta3(r[0], r[1], r[2], r[3], r[4], b[4], b[5]);
        if(i==5) return ODE_Solver_beta4(r[0], r[1], r[2], r[3], r[4], b[5]);
        if(i==6) return ODE_Solver_beta5(r[0], r[1], r[2], r[3], r[4]);

        return 0;
    }

    //===================================================================================
    //
    // extrapolation using old function values
    //
    //===================================================================================
    void ODE_Solver_extrapolate_ynp1(int order, ODE_solver_Solution &Snp1,
                                     ODE_solver_Solution *Sn[6])
    {
        if(order<1 || order>6)
        {
            cerr << " check order for ODE_Solver_extrapolate_ynp1 " << endl;
            exit(0);
        }

        //================================================================================
        // the structures yn, ynm1,... contain the information from the previous time
        // steps the structure ynp1 contains the current version of ynp1 and dynm1
        //================================================================================
        double Dznp1=Snp1.z-Sn[0]->z;
        //
        double b[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double Dzn[5]={0.0, 0.0, 0.0, 0.0, 0.0};
        double rho[5]={0.0, 0.0, 0.0, 0.0, 0.0};

        for(int k=2; k<=order; k++){ Dzn[k-2]=Sn[0]->z-Sn[k-1]->z; rho[k-2]=Dzn[k-2]/Dznp1; }

        for(int k=order; k>=2; k--) b[k-1]=ODE_Solver_betai(k, rho, b);

        b[0]=ODE_Solver_beta0(b[1], b[2], b[3], b[4], b[5]);

        // add up different orders
        for(int i=0; i<(int)Sn[0]->y.size(); i++)
        {
            double Sa=0.0;
            for(int mo=order-1; mo>=0; mo--) Sa+=b[mo]*Sn[mo]->y[i];

            Snp1.y[i]=Sa;
        }

        return;
    }

    //===================================================================================
    //
    // polynomial interpolation using old function values
    //
    //===================================================================================
    void ODE_Solver_interpolate(int order,
                                ODE_solver_Solution *Sn[6],
                                ODE_solver_Solution &Sn_interpol)
    {
        if(order<1 || order>6)
            throw_error("ODE_Solver_interpolate", "check order", 1);

        double Dzn[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double rho[6]={1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        // interpolation coefficients
        for(int k=0; k<order; k++) Dzn[k]=Sn_interpol.z-Sn[k]->z;

        for(int k=0; k<order; k++)
            for(int j=0; j<order; j++)
                if(k!=j) rho[k]*=Dzn[j]/(Sn[k]->z-Sn[j]->z);

        // add up different orders of the polynomial
        for(int i=0; i<(int)Sn[0]->y.size(); i++)
        {
            double y=0.0, dy=0.0;
            for(int mo=order-1; mo>=0; mo--)
            {
                y +=rho[mo]*Sn[mo]->y [i];
                dy+=rho[mo]*Sn[mo]->dy[i];
            }

            Sn_interpol.y [i]=y;
            Sn_interpol.dy[i]=dy;
        }

        return;
    }

    //===================================================================================
    //
    // Functions for the ODE solver (error checking)
    //
    //===================================================================================
    inline double ODE_Solver_error_check(double z, double rtol, double atol, double y,
                                         double Dy, int mflag, int component,
                                         int &error_type, string mess="")
    {
        Dy=fabs(Dy);
        error_type=0;

        double err_aim=rtol*fabs(y), fac=1.0e+10;
        if(err_aim<=atol)
        {
            if(mflag>=3) cout << " Time-stepper: absolute error control for ( " << mess
                              << " / component " << component << " ) at z= "
                              << z << " and y= " << y << " |Dy|= " << Dy
                              << " Dy_rel= " << err_aim << " Dy_abs= " << atol << endl;
            err_aim=atol;
            error_type=1;
        }

        if(Dy!=0.0) fac=err_aim/Dy;

        return fac;
    }

    //===================================================================================
    //
    // error and step-size estimate
    //
    //===================================================================================
    int ODE_Solver_estimate_error_and_next_possible_step_size(int order, const double zin,
                                                              double z, double zend,
                                                              double &Dz_z, const double Dz_z_max,
                                                              ODE_Solver_data &ODE_Solver_info,
                                                              const ODE_solver_Solution &S,
                                                              ODE_solver_Solution *Sn[6],
                                                              int messflag=1)
    {
        //--------------------------------------------------------
        // even if higher order is used just go with fourth order
        // for error check (added June 15)
        //--------------------------------------------------------
        int loc_order=(int) min(max(order, 1), 6);

        //--------------------------------------------------------
        // S contains best approximation to order 'order'
        // *Stemp will contain best approximation to order 'order-1'
        //--------------------------------------------------------
        int neq=S.y.size();
        ODE_Solver_info.Stemp.z=z;
        ODE_Solver_extrapolate_ynp1(order, ODE_Solver_info.Stemp, Sn);

        //--------------------------------------------------------
        // compute difference
        //--------------------------------------------------------
#ifdef OPENMP_ACTIVATED_ODE
#pragma omp parallel for default(shared)
#endif
        for(int i=0; i<neq; i++) ODE_Solver_info.Stemp.y[i]-=S.y[i];

        //--------------------------------------------------------
        // estimate errors
        //--------------------------------------------------------
        double power, fac=1.0e+10, dum;
        int max_err_time=-1, err_type=0;
        //--------------------------------------------------------
        if(loc_order<=1) power=1.0;
        else power=1.0/loc_order;

        for(int l=0; l<neq; l++)
        {
            if(ODE_Solver_info.tols.rel[l]!=0.0)
            {
                dum=ODE_Solver_error_check(z,
                                           ODE_Solver_info.tols.rel[l],
                                           ODE_Solver_info.tols.abs[l], S.y[l],
                                           ODE_Solver_info.Stemp.y[l], messflag,
                                           l, err_type, "?");

                if(fac>dum){ fac=dum; max_err_time=l; }
            }
        }

        if(messflag>=2 && max_err_time!=-1)
        {
            cout << " Max error in component " << max_err_time
                 << " (of " << neq << ") hitting "
                 << (err_type==0 ? "relative error " : "absolute error ")
                 << (err_type==0 ?
                     ODE_Solver_info.tols.rel[max_err_time] :
                     ODE_Solver_info.tols.abs[max_err_time])
                 << endl;
        }

        //--------------------------------------------------------
        // change step-size;
        // limit change by different values, depending on success
        //--------------------------------------------------------
        // control of step-size changes. Old quoted setting were
        // numerically unstable. Most crucial was change of
        // 'up_thresh', which controls the minimally required
        // change upwards [added 26.01.19, JC]
        //--------------------------------------------------------
        double max_fac    =5.0; // was =100.0;
        double down_thresh=0.7; // was =0.9;
        double up_thresh  =1.3; // was =1.1;
        //--------------------------------------------------------
        double Dz_z_min=1.0e-10;
        double Dz_z_new=Dz_z*min(max_fac, pow(fac, power));

        Dz_z_new=min(max(Dz_z_new, Dz_z_min), Dz_z_max);

        if(messflag==1)
        {
            cout << " Estimated change: ";
            cout << " order= " << order << " estimate order= " << loc_order
                 << " current step-size: " << Dz_z
                 << " -- new suggested step-size: " << Dz_z_new << " -- Dz= " << zin*Dz_z_new
                 << " z= " << z
                 << " zend= " << zend << " Dz= " << z-zend << " zin= " << zin << " Dz= " << zin-z;

            if(Dz_z_new==Dz_z) cout << " -- no change " << endl;
            else if(Dz_z_new==Dz_z_max) cout << " -- maximal step-size reached " << endl;
            else if(Dz_z_new>Dz_z) cout << " -- increase by: " << Dz_z_new/Dz_z << endl;
            else cout << " -- decrease by: " << Dz_z_new/Dz_z << endl;
        }

        //--------------------------------------------------------
        // to put limits on the changes in Dz (not too small)
        //--------------------------------------------------------
        if(Dz_z*down_thresh>Dz_z_new)
        {
            Dz_z=Dz_z_new;

            ODE_Solver_info.n_down=(int) min(ODE_Solver_info.n_down+1, 5);

            //--------------------------------------------------------
            // if there were several subsequent decreases of the
            // step-size take some more drastic action
            //--------------------------------------------------------
            if(ODE_Solver_info.n_down>=3)
            {
                if(ODE_Solver_info.n_down>=5)
                {
                    Dz_z=Dz_z_min;
                    cout << "\n Reset to minimal step-size " << endl;
                }
                else
                {
                    Dz_z/=pow(2.0, ODE_Solver_info.n_down-1);
                    cout << "\n Checking to decreasing the step-size by additional factor of "
                         << pow(2.0, ODE_Solver_info.n_down-1) << endl;
                }
            }

            if(Dz_z*0.5>Dz_z_new) return 1;
        }
        else if(Dz_z*up_thresh<Dz_z_new)
        {
            Dz_z=min(Dz_z_new, Dz_z_max);

            ODE_Solver_info.n_down=0;

            return 0;
        }

        ODE_Solver_info.n_down=0;
        return 0;
    }

    //===================================================================================
    void ODE_solver_compute_Jacobian_Matrix_Anisotropies(double z, double Dz,
                                                         int neq, double *y,
                                                         ODE_solver_matrix &Jac,
                                                         ODE_Solver_data &SD,
                                                         int mess=0)
    {
        //===============================================================================
        // z is the current redshift
        // y[] should contain the solution for which the jacobian is needed
        // Jac[] should have dimension neq*neq
        //===============================================================================
        if(mess==1) cout << " entering full Jacobinan update. " << endl;

        if(Jac.A.size()==0)
        {
            Jac.A.resize(neq*neq);
            Jac.row.resize(neq*neq);
            Jac.col.resize(neq*neq);
            Jac.diags.resize(neq);
        }

        //===============================================================================
        // this vector will be used to store the jacobian after calls of f_i(X)
        //===============================================================================
        vector<double> Jij(neq);

        //===============================================================================
        // fill Matrix with values
        //===============================================================================
        void (*Jac_func)(int, int, double, double, double *, double *, ODE_Solver_data &);

        if(SD.num_jac) Jac_func=ODE_solver_jac_2p;
        else Jac_func=ODE_solver_jac_user;

        for(int J=0; J<neq; J++)
        {
            Jac_func(neq, J, z, Dz, y, &Jij[0], SD);

            for(int row=0; row<neq; row++)
            {
                Jac.A[J+row*neq]=Jij[row];
                Jac.col[J+row*neq]=J;
                Jac.row[J+row*neq]=row;

                if(J==row) Jac.diags[J]=(int)(Jac.A.size()-1);
            }
        }

        if(mess==1) cout << " Number of non-zero elements = " << Jac.A.size()
            << ". Full matrix has " << neq*neq << " elements " << endl;

        if((int)Jac.diags.size()!=neq)
            cout << " Hmmm. That should not happen... " << endl;

        return;
    }

    //===================================================================================
    //
    // setting up the equation system for matrix solver
    //
    //===================================================================================
    int ODE_Solver_compute_new_Solution(int order, const double zout, const double reltol,
                                        ODE_Solver_data &ODE_Solver_info,
                                        ODE_solver_Solution * &Snew,
                                        ODE_solver_Solution *Sn[6])
    {
        //==============================================================
        double zin=Sn[0]->z, Dz=zout-zin, delta=0;
        int neq=Sn[0]->y.size(), converged=0;

        //==============================================================
        // tolerances for Jacobian convergence
        //==============================================================
        double tolJac = reltol/2.0;

        //==============================================================
        int Jac_loops=0;
        int LA_error=0;
        int mess=0;

        //==============================================================
        // set initial values for iteration: extrapolate from z--> zout
        //==============================================================
        ODE_Solver_info.Sguessptr->z=Snew->z=zout;
        ODE_Solver_extrapolate_ynp1(order, *ODE_Solver_info.Sguessptr, Sn);

        do{
            Jac_loops++;

            if(LA_error!=0)
            {
                ODE_Solver_extrapolate_ynp1(order, *ODE_Solver_info.Sguessptr, Sn);
                LA_error=0;
                if(mess>=2) cout << " resetting Solver " << endl;
            }

            //--------------------------------------------------------
            // redshift z==zn --> compute f_i(yguess, z)
            //--------------------------------------------------------
            ODE_solver_f(neq, zout,
                         &ODE_Solver_info.Sguessptr->y[0],
                         &ODE_Solver_info.Sguessptr->dy[0],
                         ODE_Solver_info);

            //--------------------------------------------------------
            // after calling this function, ynew will contain y(yguess)
            // delta is the coefficient in front of h f(y); it depends
            // on the order that is used to compute things
            //--------------------------------------------------------
            delta=ODE_Solver_compute_ynp1(order, *Snew, *ODE_Solver_info.Sguessptr, Sn);

            //--------------------------------------------------------
            // this is rhs of J_F Dx = -F
            //--------------------------------------------------------
#ifdef OPENMP_ACTIVATED_ODE
#pragma omp parallel for default(shared)
#endif
            for(int i=0; i<neq; i++)
                ODE_Solver_info.F[i]=Snew->y[i]-ODE_Solver_info.Sguessptr->y[i];

            //--------------------------------------------------------
            // compute Jacobian
            //--------------------------------------------------------
            ODE_solver_compute_Jacobian_Matrix(zout, Dz*delta, neq,
                                               &ODE_Solver_info.Sguessptr->y[0],
                                               *ODE_Solver_info.Jacobian_ptr,
                                               ODE_Solver_info);

            //--------------------------------------------------------
            // solve equation for dY
            //--------------------------------------------------------
            if(ODE_Solver_info.LA_solver=="default" ||
               ODE_Solver_info.LA_solver=="Gaussian")
                LA_error=ODE_Solver_Solve_LA_system_Gaussian(*ODE_Solver_info.Jacobian_ptr,
                                                             ODE_Solver_info.F,
                                                             ODE_Solver_info.dY, 1);

            else if(ODE_Solver_info.LA_solver=="Gaussian-GSL")
                LA_error=ODE_Solver_Solve_LA_system_GSL(*ODE_Solver_info.Jacobian_ptr,
                                                        ODE_Solver_info.F,
                                                        ODE_Solver_info.dY, 1);

            else throw_error("ODE_Solver_compute_new_Solution",
                             "LA-solver not available", 1);

            //--------------------------------------------------------
            // if nan was produced return main
            //--------------------------------------------------------
            if(LA_error==10) return 10;

            //--------------------------------------------------------
            // update current solution
            //--------------------------------------------------------
            if(LA_error==0)
#ifdef OPENMP_ACTIVATED_ODE
#pragma omp parallel for default(shared)
#endif
                for(int i=0; i<neq; i++)
                    ODE_Solver_info.Sguessptr->y[i]+=ODE_Solver_info.dY[i];

            //--------------------------------------------------------
            // check convergence of iteration
            //--------------------------------------------------------
            converged=1;

            if(LA_error==0)
            {
                for(int k=0; k<ODE_Solver_info.tols.neq; k++)
                {
                    double rel=ODE_Solver_info.tols.rel[k];
                    double abs=ODE_Solver_info.tols.abs[k];
                    double dY=fabs(ODE_Solver_info.dY[k]);
                    double dy_rel=rel*fabs(ODE_Solver_info.Sguessptr->y[k]);

                    if(rel!=0.0 && dY>=max(dy_rel, abs))
                    {
                        converged=0;
                        if(mess>=2)
                            cout << " " << k << " " << dY << " " << dy_rel << " " << abs << endl;
                        break;
                    }
                }
            }

            if(converged==0)
            {
                if(!(Jac_loops%5))
                {
                    tolJac/=2.0;
                    cout << " Tightening Jac error setting at it-# " << Jac_loops << endl;
                }

                if(Jac_loops>25)
                {
                    converged=1;
                    if(mess>=2){ cout << " Not very happy :S " << endl; wait_f_r(); }
                    return 1;
                }
            }
        }
        while(converged==0);

        if(mess>=1) cout << "\n Number of Jacobian iterations " << Jac_loops
                         << " needed for propagation from z= " << zin
                         << " by Dz= " << zin-zout << endl;

        //--------------------------------------------------------
        // swap variables, so that Snew contains new solution
        //--------------------------------------------------------
        ODE_Solver_info.Sptr=Snew;
        Snew=ODE_Solver_info.Sguessptr;
        ODE_Solver_info.Sguessptr=ODE_Solver_info.Sptr;
        ODE_Solver_info.Sptr=NULL;

        return 0;
    }

    //===================================================================================
    //
    // do time-step
    //
    //===================================================================================
    int ODE_Solver_Solve_history(double zs, double zend,
                                 double Dz_in, double Dz_max,
                                 ODE_solver_Solution &Sz,
                                 ODE_Solver_data &ODE_Solver_info)
    {
        if(ODE_Solver_info.Snewptr==NULL)
        {
            cout << " ODE_Solver_Solve_history:: please set up the memory first! Exiting"
                 << endl;
            exit(0);
        }

        //===========================================================================
        // no PDEs, no OAEs
        //===========================================================================
        set_equation_indices(0, Sz.y.size(), 0, ODE_Solver_info);

        //===========================================================================
        // define integration direction
        //===========================================================================
        if(zs>zend) ODE_Solver_info.direction=1;
        else ODE_Solver_info.direction=-1;

        //===========================================================================
        // setup
        //===========================================================================
        double zout=zs, zin=zs;
        double Dz_z=max(Dz_in/zs, ODE_Solver_info.Dz_z_last), Dz_z_max;
        int redo_run=0, LAerror=0;
        int maxorder=5;
        int messflag=ODE_Solver_info.verbosity;

        if(zs!=0.0) Dz_z_max=fabs((zend-zs)/zs);
        else Dz_z_max=fabs(zend-zs);

        if(Dz_max!=0.0 && zs!=0.0) Dz_z_max=Dz_max/zs;

        Dz_z=min(Dz_z, Dz_z_max);

        //===========================================================================
        // for initial step reduce the accuracy
        //===========================================================================
        // make sure that max order is not set larger
        ODE_Solver_info.order=(int)min(ODE_Solver_info.order, maxorder);

        //===========================================================================
        // try some small step for very first integration
        //===========================================================================
        if(ODE_Solver_info.Dz_z_last==0.0) Dz_z=ODE_Solver_info.tolSol/2.0;

        do{
            if(ODE_Solver_info.direction==1)
            {
                if(zin<0.0) zout=max(zin*(1.0+Dz_z), zend);
                else if(zin>ODE_Solver_info.tolSol) zout=max(zin*(1.0-Dz_z), zend);
                else zout=max(-ODE_Solver_info.tolSol, zend);
            }
            else
            {
                if(zin<-ODE_Solver_info.tolSol) zout=min(zin*(1.0-Dz_z), zend);
                else if(zin>0.0) zout=min(zin*(1.0+Dz_z), zend);
                else zout=min(ODE_Solver_info.tolSol, zend);
            }

            LAerror=ODE_Solver_compute_new_Solution(ODE_Solver_info.order, zout,
                                                    ODE_Solver_info.tolSol,
                                                    ODE_Solver_info,
                                                    ODE_Solver_info.Snewptr,
                                                    ODE_Solver_info.Snptr);

            if(LAerror==1)
            {
                LAerror=0;
                Dz_z/=2.0;
                Dz_z=max(Dz_z, 1.0e-10);
                zout=zin;
                redo_run=1;
            }
            else if(LAerror==10){ return 10; }
            else if(LAerror==0)
            {
                //--------------------------------------------------------
                // estimate error and next possible stepsize
                //--------------------------------------------------------
                redo_run=0;
                redo_run=ODE_Solver_estimate_error_and_next_possible_step_size(ODE_Solver_info.order,
                                                                               zin, zout, zend, Dz_z,
                                                                               Dz_z_max, ODE_Solver_info,
                                                                               *ODE_Solver_info.Snewptr,
                                                                               ODE_Solver_info.Snptr,
                                                                               messflag);

                //--------------------------------------------------------
                // accepting the current step if 'redo_run==0'
                //--------------------------------------------------------
                if(redo_run==0)
                {
                    ODE_Solver_info.Sptr=ODE_Solver_info.Snptr[5];
                    for(int k=5; k>0; k--) ODE_Solver_info.Snptr[k]=ODE_Solver_info.Snptr[k-1];
                    ODE_Solver_info.Snptr[0]=ODE_Solver_info.Snewptr;
                    ODE_Solver_info.Snewptr=ODE_Solver_info.Sptr;
                    ODE_Solver_info.Sptr=NULL;
                    ODE_Solver_info.count++;

                    if(ODE_Solver_info.count>=ODE_Solver_info.order+1)
                    {
                        ODE_Solver_info.order=(int)min(ODE_Solver_info.order+1, maxorder);
                        ODE_Solver_info.count=1;
                    }
                    //
                    zin=zout;
                }
                else{ redo_run=1; ODE_Solver_info.Jac_is_set=0; }
            }
            else{ cout << " Unknown error in LA-solver accurred: " << LAerror << endl; exit(1); }
        }
        while(ODE_Solver_info.direction*zout>ODE_Solver_info.direction*zend);

        //===========================================================================
        // if run was accepted then Snptr contains new solution!!!
        //===========================================================================
        ODE_Solver_info.Dz_z_last=Dz_z;
        //
        Sz.z =ODE_Solver_info.Snptr[0]->z;
        Sz.y =ODE_Solver_info.Snptr[0]->y;
        Sz.dy=ODE_Solver_info.Snptr[0]->dy;

        return 0;
    }

    //===================================================================================
    //
    // setting up the equation system for matrix solver
    //
    //===================================================================================
    int ODE_Solver_compute_new_Solution_Anisotropies(int order, const double zout,
                                                     const double reltol,
                                                     ODE_Solver_data &ODE_Solver_info,
                                                     ODE_solver_Solution * &Snew,
                                                     ODE_solver_Solution *Sn[6], int n_Nu)
    {
        //==============================================================
        double zin=Sn[0]->z, Dz=zout-zin, delta=0;
        int neq=Sn[0]->y.size(), converged=0;

        //==============================================================
        // tolerances for Jacobian convergence
        //==============================================================
        double tolJac = reltol/3.0;

        //==============================================================
        int Jac_loops=0;
        int LA_error=0;
        int mess=0;

        //==============================================================
        // set initial values for iteration: extrapolate from z--> zout
        //==============================================================
        ODE_Solver_info.Sguessptr->z=Snew->z=zout;
        ODE_Solver_extrapolate_ynp1(order, *ODE_Solver_info.Sguessptr, Sn);

        do{
            Jac_loops++;

            if(LA_error!=0)
            {
                ODE_Solver_extrapolate_ynp1(order, *ODE_Solver_info.Sguessptr, Sn);
                LA_error=0;
                if(mess>=2) cout << " resetting Solver " << endl;
            }

            //--------------------------------------------------------
            // redshift z==zn --> compute f_i(yguess, z)
            //--------------------------------------------------------
            ODE_solver_f(neq, zout,
                         &ODE_Solver_info.Sguessptr->y[0],
                         &ODE_Solver_info.Sguessptr->dy[0],
                         ODE_Solver_info);

            //--------------------------------------------------------
            // after calling this function, ynew will contain y(yguess)
            // delta is the coefficient in front of h f(y); it depends
            // on the order that is used to compute things
            //--------------------------------------------------------
            delta=ODE_Solver_compute_ynp1(order, *Snew, *ODE_Solver_info.Sguessptr, Sn);

            //--------------------------------------------------------
            // this is rhs of J_F Dx = -F
            //--------------------------------------------------------
#ifdef OPENMP_ACTIVATED_ODE
#pragma omp parallel for default(shared)
#endif
            for(int i=0; i<neq; i++)
                ODE_Solver_info.F[i]=Snew->y[i]-ODE_Solver_info.Sguessptr->y[i];

            //--------------------------------------------------------
            // compute Jacobian
            //--------------------------------------------------------
            ODE_solver_compute_Jacobian_Matrix_Anisotropies(zout, Dz*delta, neq,
                                                            &ODE_Solver_info.Sguessptr->y[0],
                                                            *ODE_Solver_info.Jacobian_ptr,
                                                            ODE_Solver_info);

            //--------------------------------------------------------
            // solve equation for dY
            //--------------------------------------------------------
            if(ODE_Solver_info.LA_solver=="default")
                LA_error=ODE_Solver_Solve_LA_system_Anisotropies(*ODE_Solver_info.Jacobian_ptr,
                                                                 ODE_Solver_info.F,
                                                                 ODE_Solver_info.dY, n_Nu, 1);

            else if(ODE_Solver_info.LA_solver=="Gaussian")
                LA_error=ODE_Solver_Solve_LA_system_Gaussian(*ODE_Solver_info.Jacobian_ptr,
                                                             ODE_Solver_info.F,
                                                             ODE_Solver_info.dY, 1);

            else if(ODE_Solver_info.LA_solver=="Gaussian-GSL")
                LA_error=ODE_Solver_Solve_LA_system_GSL(*ODE_Solver_info.Jacobian_ptr,
                                                        ODE_Solver_info.F,
                                                        ODE_Solver_info.dY, 1);

            else throw_error("ODE_Solver_compute_new_Solution_Anisotropies",
                             "LA-solver not available", 1);

            //--------------------------------------------------------
            // if nan was produced return main
            //--------------------------------------------------------
            if(LA_error==10) return 10;

            //--------------------------------------------------------
            // update current solution
            //--------------------------------------------------------
            if(LA_error==0)
#ifdef OPENMP_ACTIVATED_ODE
#pragma omp parallel for default(shared)
#endif
                for(int i=0; i<neq; i++)
                    ODE_Solver_info.Sguessptr->y[i]+=ODE_Solver_info.dY[i];

            //--------------------------------------------------------
            // check convergence of iteration
            //--------------------------------------------------------
            converged=1;

            if(LA_error==0)
            {
                for(int k=0; k<ODE_Solver_info.tols.neq; k++)
                {
                    double rel=ODE_Solver_info.tols.rel[k];
                    double abs=ODE_Solver_info.tols.abs[k];
                    double dY=fabs(ODE_Solver_info.dY[k]);
                    double dy_rel=rel*fabs(ODE_Solver_info.Sguessptr->y[k]);

                    if(rel!=0.0 && dY>=max(dy_rel, abs))
                    {
                        converged=0;
                        if(mess>=2)
                            cout << " " << k << " " << dY << " " << dy_rel << " " << abs << endl;
                        break;
                    }
                }
            }

            if(converged==0)
            {
                if(!(Jac_loops%5))
                {
                    tolJac/=2.0;
                    cout << " Tightening Jac error setting at it-# " << Jac_loops << endl;
                }

                if(Jac_loops>25)
                {
                    converged=1;
                    if(mess>=2){ cout << " Not very happy :S " << endl; wait_f_r(); }
                    return 1;
                }
            }
        }
        while(converged==0);

        if(mess>=1) cout << "\n Number of Jacobian iterations " << Jac_loops
            << " needed for propagation from z= " << zin
            << " by Dz= " << zin-zout << endl;

        //--------------------------------------------------------
        // swap variables, so that Snew contains new solution
        //--------------------------------------------------------
        ODE_Solver_info.Sptr=Snew;
        Snew=ODE_Solver_info.Sguessptr;
        ODE_Solver_info.Sguessptr=ODE_Solver_info.Sptr;
        ODE_Solver_info.Sptr=NULL;

        return 0;
    }

    //===================================================================================
    //
    // do time-step
    //
    //===================================================================================
    int ODE_Solver_Solve_history_Anisotropies(double zs, double zend,
                                              double Dz_in, double Dz_max,
                                              ODE_solver_Solution &Sz,
                                              ODE_Solver_data &ODE_Solver_info, int n_Nu)
    {
        if(ODE_Solver_info.Snewptr==NULL)
        {
            cout << " ODE_Solver_Solve_history:: please set up the memory first! Exiting"
                 << endl;
            exit(0);
        }

        //===========================================================================
        // no PDEs, no OAEs
        //===========================================================================
        set_equation_indices(0, Sz.y.size(), 0, ODE_Solver_info);

        //===========================================================================
        // define integration direction
        //===========================================================================
        if(zs>zend) ODE_Solver_info.direction=1;
        else ODE_Solver_info.direction=-1;

        //===========================================================================
        // setup
        //===========================================================================
        double zout=zs, zin=zs;
        double Dz_z=max(Dz_in/zs, ODE_Solver_info.Dz_z_last), Dz_z_max;
        int redo_run=0, LAerror=0;
        int maxorder=5;
        int messflag=ODE_Solver_info.verbosity;

        if(zs!=0.0) Dz_z_max=fabs((zend-zs)/zs);
        else Dz_z_max=fabs(zend-zs);

        if(Dz_max!=0.0 && zs!=0.0) Dz_z_max=Dz_max/zs;

        Dz_z=min(Dz_z, Dz_z_max);

        //===========================================================================
        // for initial step reduce the accuracy
        //===========================================================================
        // make sure that max order is not set larger
        ODE_Solver_info.order=(int)min(ODE_Solver_info.order, maxorder);

        //===========================================================================
        // try some small step for very first integration
        //===========================================================================
        if(ODE_Solver_info.Dz_z_last==0.0) Dz_z=ODE_Solver_info.tolSol/2.0;

        do{
            if(ODE_Solver_info.direction==1)
            {
                if(zin<0.0) zout=max(zin*(1.0+Dz_z), zend);
                else if(zin>ODE_Solver_info.tolSol) zout=max(zin*(1.0-Dz_z), zend);
                else zout=max(-ODE_Solver_info.tolSol, zend);
            }
            else
            {
                if(zin<-ODE_Solver_info.tolSol) zout=min(zin*(1.0-Dz_z), zend);
                else if(zin>0.0) zout=min(zin*(1.0+Dz_z), zend);
                else zout=min(ODE_Solver_info.tolSol, zend);
            }

            LAerror=ODE_Solver_compute_new_Solution_Anisotropies(ODE_Solver_info.order, zout,
                                                                 ODE_Solver_info.tolSol,
                                                                 ODE_Solver_info,
                                                                 ODE_Solver_info.Snewptr,
                                                                 ODE_Solver_info.Snptr, n_Nu);

            if(LAerror==1)
            {
                LAerror=0;
                Dz_z/=2.0;
                Dz_z=max(Dz_z, 1.0e-10);
                zout=zin;
                redo_run=1;
            }
            else if(LAerror==10){ return 10; }
            else if(LAerror==0)
            {
                //--------------------------------------------------------
                // estimate error and next possible stepsize
                //--------------------------------------------------------
                redo_run=0;
                redo_run=ODE_Solver_estimate_error_and_next_possible_step_size(ODE_Solver_info.order,
                                                                               zin, zout, zend, Dz_z,
                                                                               Dz_z_max, ODE_Solver_info,
                                                                               *ODE_Solver_info.Snewptr,
                                                                               ODE_Solver_info.Snptr,
                                                                               messflag);

                redo_run=0;

                //--------------------------------------------------------
                // accepting the current step if 'redo_run==0'
                //--------------------------------------------------------
                if(redo_run==0)
                {
                    ODE_Solver_info.Sptr=ODE_Solver_info.Snptr[5];
                    for(int k=5; k>0; k--) ODE_Solver_info.Snptr[k]=ODE_Solver_info.Snptr[k-1];
                    ODE_Solver_info.Snptr[0]=ODE_Solver_info.Snewptr;
                    ODE_Solver_info.Snewptr=ODE_Solver_info.Sptr;
                    ODE_Solver_info.Sptr=NULL;
                    ODE_Solver_info.count++;

                    if(ODE_Solver_info.count>=ODE_Solver_info.order+1)
                    {
                        ODE_Solver_info.order=(int)min(ODE_Solver_info.order+1, maxorder);
                        ODE_Solver_info.count=1;
                    }
                    //
                    zin=zout;
                }
                else{ redo_run=1; ODE_Solver_info.Jac_is_set=0; }
            }
            else{ cout << " Unknown error in LA-solver accurred: " << LAerror << endl; exit(1); }
        }
        while(ODE_Solver_info.direction*zout>ODE_Solver_info.direction*zend);

        //===========================================================================
        // if run was accepted then Snptr contains new solution!!!
        //===========================================================================
        ODE_Solver_info.Dz_z_last=Dz_z;
        //
        Sz.z =ODE_Solver_info.Snptr[0]->z;
        Sz.y =ODE_Solver_info.Snptr[0]->y;
        Sz.dy=ODE_Solver_info.Snptr[0]->dy;

        return 0;
    }

    //===================================================================================
    //
    // Check users Jacobian using numerical derivatives
    //
    //===================================================================================
    void check_Jacobian_setup(int neq, int J, double z, double Dz, double *y,
                              const vector<double> &Jij, ODE_Solver_data &SD)
    {
        vector<double> JijN(neq);
        double eps=SD.check_Jacobian_eps;

        ODE_solver_jac_2p(neq, J, z, Dz, y, &JijN[0], SD);
        //ODE_solver_jac_4p(neq, J, z, Dz, y, &JijN[0], SD);

        for(int row=0; row<neq; row++)
//        if(fabs(JijN[row]-Jij[row])>fabs(eps*Jij[row]) && fabs(Jij[row])>1.0e-16)
        if(fabs(JijN[row]-Jij[row])>fabs(eps*Jij[row]) && fabs(JijN[row])>1.0e-16)
        {
            cout << " Element (" << row << ", " << J
                 << ") inconsistent at relative level >~" << eps << endl;
            cout << " User == " << Jij[row] << " || Numerical == " << JijN[row] << endl;
            cout << " neq == " << neq << " nODE == " << SD.nODE << endl;
            cout << " k_ODE == " << SD.index_ODE << " k_OAE == " << SD.index_OAE << endl;

            //wait_f_r();
        }

        return;
    }

    //===================================================================================
    //
    // Jacobian setup routine for PDE + ODE solver; it is assumed that jacobian setup is
    // performed by user, but here the structure is saved.
    //
    //===================================================================================
    void ODE_solver_compute_Jacobian_Matrix_PDE_ODE_OAE_full(double z, double Dz,
                                                             int neq, double *y,
                                                             ODE_solver_matrix &Jac,
                                                             ODE_Solver_data &SD, int mess=0)
    {
        //===============================================================================
        // z is the current redshift
        // y[] should contain the solution for which the jacobian is needed
        // Jac[] should have dimension neq*neq
        //===============================================================================
        if(mess==1) cout << " entering full Jacobinan update. " << endl;

        Jac.clear();

        //===============================================================================
        // this vector will be used to store the jacobian after calls of f_i(X)
        //===============================================================================
        vector<double> Jij(neq);

        //===============================================================================
        // fill Matrix with values (PDE part)
        //===============================================================================
        void (*Jac_func)(int, int, double, double, double *, double *, ODE_Solver_data &);

        if(SD.num_jac) Jac_func=ODE_solver_jac_2p;
        else Jac_func=ODE_solver_jac_user;

        for(int J=0; J<SD.index_ODE; J++)
        {
            Jac_func(neq, J, z, Dz, y, &Jij[0], SD);

            if(SD.check_jacobian && !SD.num_jac && z<=SD.check_z)
                check_Jacobian_setup(neq, J, z, Dz, y, Jij, SD);

            for(int row=0; row<SD.index_ODE; row++)
                if(Jij[row]!=0.0
                   // always save lower boundary condition
                   || (row<=2 && J<=4)
                   // always save upper boundary condition
                   || (row>=SD.index_ODE-3 && row<=SD.index_ODE-1 && SD.index_ODE-5<=J))

                    Jac.save_info(neq, J, row, Jij[row]);

            //===========================================================================
            // for ode part ALL elements are saved
            //===========================================================================
            if(SD.do_check_nonzero_ODE)
            {
                for(int row=SD.index_ODE; row<neq; row++)
                    if(Jij[row]!=0.0) Jac.save_info(neq, J, row, Jij[row]);
            }
            else for(int row=SD.index_ODE; row<neq; row++) Jac.save_info(neq, J, row, Jij[row]);
        }

        //===============================================================================
        // fill Matrix with values (ODE + OAE parts)
        //===============================================================================
        for(int J=SD.index_ODE; J<neq; J++)
        {
            Jac_func(neq, J, z, Dz, y, &Jij[0], SD);

            if(SD.check_jacobian && !SD.num_jac && z<=SD.check_z)
                check_Jacobian_setup(neq, J, z, Dz, y, Jij, SD);

            //===========================================================================
            // for ode part ALL elements are saved
            //===========================================================================
            if(SD.do_check_nonzero_ODE)
            {
                for(int row=0; row<neq; row++)
                    if(Jij[row]!=0.0) Jac.save_info(neq, J, row, Jij[row]);
            }
            else for(int row=0; row<neq; row++) Jac.save_info(neq, J, row, Jij[row]);
        }

        if(mess>=1) cout << " Number of non-zero elements = " << Jac.A.size()
                         << ". Full matrix has " << neq*neq << " elements."
                         << " The sparseness is " << 1.0*Jac.A.size()/(neq*neq) << endl;

        if((int)Jac.diags.size()!=neq) cerr << " Hmmm. That should not happen... " << endl;

        return;
    }

    //===================================================================================
    void ODE_solver_compute_Jacobian_Matrix_PDE_ODE_OAE(double z, double Dz,
                                                        int neq, double *y,
                                                        ODE_solver_matrix &Jac,
                                                        ODE_Solver_data &SD, int mess=0)
    {
        if(Jac.A.size()==0)
        {
            ODE_solver_compute_Jacobian_Matrix_PDE_ODE_OAE_full(z, Dz, neq, y, Jac, SD, mess);
            return;
        }

        //===============================================================================
        // z is the current redshift
        // y[] should contain the solution for which the jacobian is needed
        // Jac[] should have dimension neq*neq
        //===============================================================================
        if(mess==1) cout << " entering Jacobinan update using previous structure. " << endl;

        //===============================================================================
        // this vector will be used to store the jacobian after calls of f_i(X)
        //===============================================================================
        unsigned int i=0;
        vector<double> Jij(neq);

        //===============================================================================
        // fill Matrix with values
        //===============================================================================
        void (*Jac_func)(int, int, double, double, double *, double *, ODE_Solver_data &);

        if(SD.num_jac) Jac_func=ODE_solver_jac_2p;
        else Jac_func=ODE_solver_jac_user;

        for(int J=0; J<neq; J++)
        {
            Jac_func(neq, J, z, Dz, y, &Jij[0], SD);

            if(SD.check_jacobian && !SD.num_jac && z<=SD.check_z)
                check_Jacobian_setup(neq, J, z, Dz, y, Jij, SD);

            for(; i<Jac.A.size() && J==Jac.col[i]; i++) Jac.A[i]=Jij[Jac.row[i]];
        }

        if(mess==1) cout << " Number of non-zero elements = " << Jac.A.size()
                         << ". Full matrix has " << neq*neq << " elements " << endl;

        if((int)Jac.diags.size()!=neq)
            cout << " Hmmm. That should not happen... " << endl;

        //Jac.show_whole_matrix();
        //wait_f_r();

        return;
    }

    //===================================================================================
    //
    // setting up the equation system for matrix solver
    //
    //===================================================================================
    int ODE_Solver_compute_new_Solution_PDE_ODE_OAE(int order,
                                                    const double zout,
                                                    const double reltol,
                                                    ODE_Solver_data &ODE_Solver_info,
                                                    ODE_solver_Solution * &Snew,
                                                    ODE_solver_Solution *Sn[6])
    {
        //===============================================================================
        double zin=Sn[0]->z, Dz=zout-zin, delta=0;
        int neq=Sn[0]->y.size(), converged=0;

        //===============================================================================
        int Jac_loops=0, Jac_loops_min=1, Jac_loops_max=10;
        int mess=ODE_Solver_info.verbosity-1;

        //===============================================================================
        // set initial values for iteration: extrapolate from z --> zout to get yguess.
        // This is only based on previous solutions of the problem & does not depend on
        // whether one is dealing with ODE, PDE or OAE.
        //===============================================================================
        ODE_Solver_info.Sguessptr->z=Snew->z=zout;
        ODE_Solver_extrapolate_ynp1(order, *ODE_Solver_info.Sguessptr, Sn);

        do{
            Jac_loops++;

            //---------------------------------------------------------------------------
            // redshift z==zn --> compute f_i(yguess, z)
            //---------------------------------------------------------------------------
            ODE_solver_f(neq, zout,
                         &ODE_Solver_info.Sguessptr->y[0],
                         &ODE_Solver_info.Sguessptr->dy[0],
                         ODE_Solver_info);

            //---------------------------------------------------------------------------
            // after calling this function, ynew will contain y(yguess)
            // delta is the coefficient in front of h f(y); it depends
            // on the order that is used to compute things.
            //---------------------------------------------------------------------------
            // This is like evaluating Eq. (2) of CVD 2010 (1003.4928)
            // using the initial guess. Generally Snew.y != yguess, so
            // system has to be solved next.
            //---------------------------------------------------------------------------
            delta=ODE_Solver_compute_ynp1(order, *Snew, *ODE_Solver_info.Sguessptr, Sn);

            //---------------------------------------------------------------------------
            // this is rhs of J_F Dx = -F for PDE + ODE parts
            //---------------------------------------------------------------------------
            for(int i=0; i<ODE_Solver_info.index_OAE; i++)
                ODE_Solver_info.F[i]=Snew->y[i]-ODE_Solver_info.Sguessptr->y[i];

            //---------------------------------------------------------------------------
            // OAE part
            //---------------------------------------------------------------------------
            for(int i=ODE_Solver_info.index_OAE; i<neq; i++)
                ODE_Solver_info.F[i]=-ODE_Solver_info.Sguessptr->dy[i];

            //---------------------------------------------------------------------------
            // convert the boundary conditions of the PDE
            //---------------------------------------------------------------------------
            if(ODE_Solver_info.nPDE>2)
            {
                ODE_Solver_info.F[0]=-ODE_Solver_info.Sguessptr->dy[0];
                int ui=ODE_Solver_info.index_ODE-1;
                ODE_Solver_info.F[ui]=-ODE_Solver_info.Sguessptr->dy[ui];
            }

            //---------------------------------------------------------------------------
            // compute Jacobian
            //---------------------------------------------------------------------------
            ODE_solver_compute_Jacobian_Matrix_PDE_ODE_OAE(zout, Dz*delta, neq,
                                                           &ODE_Solver_info.Sguessptr->y[0],
                                                           *ODE_Solver_info.Jacobian_ptr,
                                                           ODE_Solver_info);

            for(int it=0; it<=ODE_Solver_info.max_it; it++)
            {
                //-----------------------------------------------------------------------
                // solve equation for dY
                //-----------------------------------------------------------------------
                if(ODE_Solver_info.LA_solver=="default" ||
                   ODE_Solver_info.LA_solver=="Gaussian")
                    ODE_Solver_Solve_LA_system_Gaussian(*ODE_Solver_info.Jacobian_ptr,
                                                        ODE_Solver_info.F,
                                                        ODE_Solver_info.dY, mess);

                else if(ODE_Solver_info.LA_solver=="BICG")
                    ODE_Solver_Solve_LA_system(*ODE_Solver_info.Jacobian_ptr,
                                               ODE_Solver_info.F,
                                               //ODE_Solver_info.dY, reltol/5.0, mess);
                                               ODE_Solver_info.dY, 1.0e-10, mess);

                else if(ODE_Solver_info.LA_solver=="Gaussian-GSL")
                    ODE_Solver_Solve_LA_system_GSL(*ODE_Solver_info.Jacobian_ptr,
                                                   ODE_Solver_info.F,
                                                   ODE_Solver_info.dY, mess);

                //-----------------------------------------------------------------------
                // get f(y)-f0(y)
                //-----------------------------------------------------------------------
                if(it<ODE_Solver_info.max_it)
                {
                    ODE_solver_f(neq, zout,
                                 &ODE_Solver_info.Sguessptr->y[0],
                                 &ODE_Solver_info.Fcorr1[0], -2,
                                 ODE_Solver_info);
                }

                //-----------------------------------------------------------------------
                // update current solution
                //-----------------------------------------------------------------------
#ifdef OPENMP_ACTIVATED_ODE
#pragma omp parallel for default(shared)
#endif
                for(int i=0; i<neq; i++)
                    ODE_Solver_info.Sguessptr->y[i]+=ODE_Solver_info.dY[i];

                //-----------------------------------------------------------------------
                // add correction from full matrix and solve again
                //-----------------------------------------------------------------------
                if(it<ODE_Solver_info.max_it)
                {
                    // compute delta fb(Y^(0))
                    ODE_solver_f(neq, zout,
                                 &ODE_Solver_info.Sguessptr->y[0],
                                 &ODE_Solver_info.Fcorr2[0], -2,
                                 ODE_Solver_info);

                    for(int i=0; i<neq; i++)
                    {
                        ODE_Solver_info.F[i] = ODE_Solver_info.Fcorr2[i]-ODE_Solver_info.Fcorr1[i];
                        ODE_Solver_info.F[i]*=-Dz*delta;

                        //cout << ODE_Solver_info.F[i] << " ";
                        //cout << ODE_Solver_info.Fcorr1[i] << " " << ODE_Solver_info.Fcorr2[i] << endl;
                        //ODE_Solver_info.F[i]=0.0;
                    }

                    //cout << endl;
                    //wait_f_r();
                }
            }

            //---------------------------------------------------------------------------
            // check convergence of iteration
            //---------------------------------------------------------------------------
            converged=1;

            // if(LA_error==0) JC: TODO: check if this may be needed and better...
            for(int k=0; k<ODE_Solver_info.tols.neq; k++)
            {
                double rel=ODE_Solver_info.tols.rel[k];
                double abs=ODE_Solver_info.tols.abs[k];
                double dY=fabs(ODE_Solver_info.dY[k]);
                double dy_rel=rel*fabs(ODE_Solver_info.Sguessptr->y[k]);

                if(rel!=0.0 && dY>=max(dy_rel, abs))
                {
                    converged=0;
                    if(mess>=2) // && zout < 1500)
                    {
                        cout << neq << " " << ODE_Solver_info.nODE << " ";

                        if(k<ODE_Solver_info.index_ODE) cout << " PDE: " << k << endl;
                        else cout << " ODE/OAE: " << k-ODE_Solver_info.index_ODE << endl;

                        cout << " errors= " << dY << " " << dy_rel << " " << abs << " ";
                        cout << ODE_Solver_info.order << endl;
                    }

                    break;
                }
            }

            if(Jac_loops<Jac_loops_min) converged=0;

            // JC: TODO: there currently is no tightening of the tol here. Is this optimal?
            if(Jac_loops>=Jac_loops_max)
            {
                converged=1;
                if(mess>=1) cout << " Reach maximal number of iterations ( "
                                 << Jac_loops_max << " )" << endl;
            }

            if(mess>=1) cout << "\033[2K" << " zout= " << zout
                             << " it= " << Jac_loops << endl << "\033[1A";
        }
        while(converged==0);

        if(mess>=1) cout << " one step done " << zin-zout << " " << zin << " " << zout << endl;

        if(mess>=1) cout << "\n Number of Jacobian iterations " << Jac_loops
                         << " needed for propagation from z= " << zin
                         << " by Dz= " << zin-zout
                         << " or Dz/z= " << zin/zout-1.0 << endl;

        //===============================================================================
        // swap variables, so that Snew contains new solution
        //===============================================================================
        ODE_Solver_info.Sptr=Snew;
        Snew=ODE_Solver_info.Sguessptr;
        ODE_Solver_info.Sguessptr=ODE_Solver_info.Sptr;
        ODE_Solver_info.Sptr=NULL;

        return 0;
    }

    //==============================================================================================
    // This routine has can only be called after the setup was done
    // (ODE_Solver_set_up_solution_and_memory)
    //
    // After the cal ODE_solver_Solution Sz will contain the solution at time zend; initially Sz
    // should contain the solution at zs. Previous solutions will also be stored by the routine.
    //
    // Furthermore, this routine allows to treat PDE + ODE + OAE systems. The PDE is assumed to come
    // first while the next nODE elements are the ODEs followed by OAE ordinary equations. The
    // boundary conditions are assumed to be given as algebraic equations in the setup routines.
    // This is the assumes structure: bc, 1, 2, ..., nPDE-3, nPDE-2, bc, nODE, nOAE.
    //
    //==============================================================================================
    // Apr 12th, 2019: overshooting allowed [BB & JC]
    // Oct 15th, 2013: added adaptive step-size control; Dz_in < Dz < Dz_max will be ensured;
    // to use this option make sure that NOT all relative errors are == 0
    // if this option is not required just set Dz_in==Dz_max (error control will have no effect)
    //==============================================================================================
    int ODE_Solver_Solve_history_PDE_ODE_OAE_run(double zs, double zend, double zcrit,
                                                 double Dz_in, double Dz_max,
                                                 ODE_Solver_data &ODE_Solver_info,
                                                 int nODE, int nOAE, int setmaxorder)
    {
        //===========================================================================
        // setup
        //===========================================================================
        double zout=zs, zin=zs;
        double Dz_z=max(Dz_in/zs, ODE_Solver_info.Dz_z_last), Dz_z_max;
        int LAerror=0;
        int maxorder=setmaxorder;
        int messflag=ODE_Solver_info.verbosity, nsteps=0;

        //===========================================================================
        // set initial step-size
        //===========================================================================
        if(zs!=0.0) Dz_z_max=fabs((zcrit-zs)/zs);
        else Dz_z_max=fabs(zcrit-zs);

        if(Dz_max!=0.0 && zs!=0.0) Dz_z_max=Dz_max/zs;

        if(Dz_z==0.0) cerr << " please give non-zero step-size" << endl;

        Dz_z=min(Dz_z, Dz_z_max);

        //===========================================================================
        // make sure that max order is not set larger
        //===========================================================================
        ODE_Solver_info.order=(int)min(ODE_Solver_info.order, maxorder);

        do{
            if(ODE_Solver_info.direction==1)
            {
                if(zin<0.0) zout=max(zin*(1.0+Dz_z), zcrit);
                else if(zin>ODE_Solver_info.tolSol) zout=max(zin*(1.0-Dz_z), zcrit);
                else zout=max(-ODE_Solver_info.tolSol, zcrit);
            }
            else
            {
                if(zin<-ODE_Solver_info.tolSol) zout=min(zin*(1.0-Dz_z), zcrit);
                else if(zin>0.0) zout=min(zin*(1.0+Dz_z), zcrit);
                else zout=min(ODE_Solver_info.tolSol, zcrit);
            }
            // avoid taking tiny last step [JC Aug 2018]
            if(abs(zcrit-zout)<=Dz_in*1.0e-2) zout=zcrit;

            LAerror=ODE_Solver_compute_new_Solution_PDE_ODE_OAE(ODE_Solver_info.order, zout,
                                                                ODE_Solver_info.tolSol,
                                                                ODE_Solver_info,
                                                                ODE_Solver_info.Snewptr,
                                                                ODE_Solver_info.Snptr);

            nsteps++;                  // count steps/attempts
            Dz_z=fabs(zout/zin-1.0);   // save actual Dz_z

            if(LAerror==1)
            {
                LAerror=0;
                Dz_z/=2.0;
                Dz_z=max(Dz_z, 1.0e-10);
                zout=zin;
            }
            else if(LAerror==0)
            {
                //--------------------------------------------------------
                // estimate error and next possible stepsize
                //--------------------------------------------------------
                ODE_Solver_estimate_error_and_next_possible_step_size(ODE_Solver_info.order,
                                                                      zin, zout, zcrit, Dz_z,
                                                                      Dz_z_max, ODE_Solver_info,
                                                                      *ODE_Solver_info.Snewptr,
                                                                      ODE_Solver_info.Snptr,
                                                                      messflag);

                Dz_z=max(Dz_z, Dz_in/zout);

                //--------------------------------------------------------
                // provide some feedback
                //--------------------------------------------------------
                if(messflag>1)
                {
                    cout << " Last attempted Dz_z = " << Dz_z << endl;
                    cout << " Maximal Dz_z = " << fabs(Dz_max/zin) << endl;
                    cout << " Minimal Dz_z = " << fabs(Dz_in/zin) << endl;
                }

                //--------------------------------------------------------
                // accepting the current step
                //--------------------------------------------------------
                ODE_Solver_info.Sptr=ODE_Solver_info.Snptr[5];
                for(int k=5; k>0; k--) ODE_Solver_info.Snptr[k]=ODE_Solver_info.Snptr[k-1];
                ODE_Solver_info.Snptr[0]=ODE_Solver_info.Snewptr;
                ODE_Solver_info.Snewptr=ODE_Solver_info.Sptr;
                ODE_Solver_info.Sptr=NULL;
                ODE_Solver_info.count++;

                if(ODE_Solver_info.count>=ODE_Solver_info.order+1)
                {
                    ODE_Solver_info.order=(int)min(ODE_Solver_info.order+1, maxorder);
                    ODE_Solver_info.count=1;
                }
                //
                zin=zout;
            }
            else{ cout << " Unknown error in LA-solver accurred: " << LAerror << endl; exit(1); }
        }
        while(ODE_Solver_info.direction*zout>ODE_Solver_info.direction*zend);

        //===========================================================================
        // save last step-size
        //===========================================================================
        ODE_Solver_info.Dz_z_last=Dz_z;

        if(messflag>1)
        {
            cout << " Number of steps/attempts = " << nsteps << endl;
            cout << " Next starting Dz_z = " << Dz_z << endl;
            cout << " maximal Dz_z = " << Dz_max << endl;
        }

        return 0;
    }


    //==============================================================================================
    // driver routine with overshooting up to z==zcrit
    //==============================================================================================
    int ODE_Solver_Solve_history_PDE_ODE_OAE(double zs, double zend, double zcrit,
                                             double Dz_in, double Dz_max,
                                             ODE_solver_Solution &Sz,
                                             ODE_Solver_data &ODE_Solver_info,
                                             int nODE, int nOAE, int setmaxorder)
    {
        if(ODE_Solver_info.Snewptr==NULL)
            throw_error("ODE_Solver_Solve_history_PDE_ODE_OAE",
                        "please set up the memory first! Exiting.", 1);

        if(ODE_Solver_info.verbosity>0) cout << " computing solution " << endl;

        //===========================================================================
        // storing number of PDEs
        //===========================================================================
        int neq=Sz.y.size();
        set_equation_indices(neq-nODE-nOAE, nODE, nOAE, ODE_Solver_info);

        //===========================================================================
        // define integration direction
        //===========================================================================
        double zlast=ODE_Solver_info.Snptr[0]->z;
        if(zs>zend) ODE_Solver_info.direction=1;
        else ODE_Solver_info.direction=-1;

        if(ODE_Solver_info.direction*zs<ODE_Solver_info.direction*zlast)
            throw_error("ODE_Solver_Solve_history_PDE_ODE_OAE",
                        "this should never happen", 2);

        //===========================================================================
        // check if run is needed [if order == 1 --> only one solution known]
        //===========================================================================
        if(ODE_Solver_info.direction*zend<ODE_Solver_info.direction*zlast ||
           ODE_Solver_info.order==1)
        {
            if(ODE_Solver_info.verbosity>0) cout << " calling solver " << endl;

            ODE_Solver_Solve_history_PDE_ODE_OAE_run(zlast, zend, zcrit,
                                                     Dz_in, Dz_max,
                                                     ODE_Solver_info,
                                                     nODE, nOAE, setmaxorder);
        }

        //===========================================================================
        // interpolate or just copy solution
        //===========================================================================
        if(zend==ODE_Solver_info.Snptr[0]->z)
        {
            if(ODE_Solver_info.verbosity>0) cout << " Copying solution " << endl;

            Sz.z=zend;
            Sz.y =ODE_Solver_info.Snptr[0]->y;
            Sz.dy=ODE_Solver_info.Snptr[0]->dy;
        }
        else
        {
            if(ODE_Solver_info.verbosity>0) cout << " Interpolating solution " << endl;

            ODE_Solver_Interpolate_Solution(zend, ODE_Solver_info, Sz);
        }

        return 0;
    }

    //==============================================================================================
    // driver routine without overshooting
    //==============================================================================================
    int ODE_Solver_Solve_history_PDE_ODE_OAE(double zs, double zend,
                                             double Dz_in, double Dz_max,
                                             ODE_solver_Solution &Sz,
                                             ODE_Solver_data &ODE_Solver_info,
                                             int nODE, int nOAE, int setmaxorder)
    {
        // call overshooting routine with zcrit == zend
        return ODE_Solver_Solve_history_PDE_ODE_OAE(zs, zend, zend, Dz_in, Dz_max,
                                                    Sz, ODE_Solver_info,
                                                    nODE, nOAE, setmaxorder);
    }

    //==============================================================================================
    // to obtain interpolated solution between zs and ze. The solver should have been run before
    // and the solution at ze should be available already.
    //==============================================================================================
    int ODE_Solver_Interpolate_Solution(double z,
                                        ODE_Solver_data &ODE_Solver_info,
                                        ODE_solver_Solution &Sn)
    {
        double zend=ODE_Solver_info.Snptr[0]->z;
        double zs  =ODE_Solver_info.Snptr[ODE_Solver_info.order-1]->z;

        if(!(ODE_Solver_info.direction*zs>=ODE_Solver_info.direction*z &&
             ODE_Solver_info.direction*z>=ODE_Solver_info.direction*zend))

            throw_error("ODE_Solver_Interpolate_Solution", "not for extrapolation!", 1);

        Sn.z=z;
        ODE_Solver_interpolate(ODE_Solver_info.order, ODE_Solver_info.Snptr, Sn);

        //-----------------------------------------------------------------------
        // call ODE setup on interpolated solution to get dy instead of interpol
        // [JC: this also ensures that the ODE settings of user are updated...]
        //-----------------------------------------------------------------------
        ODE_solver_f(Sn.y.size(), z, &Sn.y[0], &Sn.dy[0], ODE_Solver_info);

        if(ODE_Solver_info.verbosity>1)
        {
            cout << " zs= " << zs << " z= " << z << " zend= " << zend << endl;
            if(z<1.0e+4) wait_f_r();
        }

        return 0;
    }
}

//==================================================================================================
//==================================================================================================
