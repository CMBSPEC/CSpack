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
// Last modification   : 01/03/2020
//==================================================================================================
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

#ifndef _ODE_SOLVER_REC_H_
#define _ODE_SOLVER_REC_H_

#include <vector>
#include <string>
#include "ODE_solver_LA.h"

using namespace std;
using namespace ODE_solver_LA;

namespace ODE_solver
{
    //==============================================================================================
    // Structure that contains the solution y and its derivative dy at redshift z
    //==============================================================================================
    struct ODE_solver_Solution
    {
        int neq;
        double z;
        vector<double> y, dy;

        void allocate_memory(int neqv) { neq=neqv; y.resize(neq); dy.resize(neq); }
        void clear() { neq=0; z=0.0; y.clear(); dy.clear(); }
        ODE_solver_Solution(){ neq=0; z=0.0; }
        ODE_solver_Solution(int neqv){ neq=neqv; allocate_memory(neqv); }
    };

    //==============================================================================================
    struct ODE_solver_accuracies
    {
        int neq;
        vector<double> rel, abs;

        void allocate_memory(int neqv) { neq=neqv; rel.resize(neq); abs.resize(neq); }
        void clear() { neq=0; rel.clear(); abs.clear(); }
        ODE_solver_accuracies(){ neq=0; }
        ODE_solver_accuracies(int neqv){ neq=neqv; allocate_memory(neqv); }
    };

    //==============================================================================================
    struct ODE_Solver_data
    {
        int verbosity;

        void (*fcn_col_ptr)(int *neq, double *z, double *y, double *f, int col);
        void (*fcn_col_para_ptr)(int *neq, double *z, double *y, double *f, int col, void *p);

        bool num_jac;
        bool check_jacobian; // to check consistency of Jacobian
        double check_Jacobian_eps, check_z;
        bool do_check_nonzero_ODE; // check if ODE part entries are !=0.0
        string LA_solver;

        ODE_solver_Solution Stemp;
        ODE_solver_Solution *Sptr;
        ODE_solver_Solution *Snewptr, *Sguessptr;
        ODE_solver_Solution *Snptr[6];

        ODE_solver_accuracies tols;

        ODE_solver_matrix *Jacobian_ptr;

        vector<double> F, dY, Fcorr1, Fcorr2;

        double Dz_z_last;
        double zstart_solver;
        double tolSol;
        int order;
        int count;
        int Jac_is_set;
        // counting number of significant up and down-steps
        int n_up, n_down;
        double direction;
        //
        int nPDE; // number of pde points (including boundary condition)
        int nODE; // number of odes
        int nOAE; // number of ades
        int neq;  // total number of equations
        int index_ODE, index_OAE;
        // to pass on parameters for the setup
        void *p;
        int max_it; // for iterative solution of full matrix problem

        void reset_Solver_Info()
        {
            Snewptr=NULL; // this indicates a fresh start is happening
            return;
        }

        // default values
        ODE_Solver_data()
        {
            verbosity=0;
            reset_Solver_Info();
            p=NULL;
            nPDE=nODE=nOAE=neq=index_ODE=index_OAE=0;
            fcn_col_ptr=NULL;
            fcn_col_para_ptr=NULL;
            check_jacobian=0;    // to check consistency of Jacobian
            check_z=0.0;         // check performed at z<=check_z
            check_Jacobian_eps=0.0;
            do_check_nonzero_ODE=0;
            max_it=0;

            LA_solver="default"; // uses pre-selected default
                                 // other options (availability depends on selected ODE solver)
                                 // 'Gaussian': Gaussian elimination learning structure
                                 // 'Gaussian-GSL': Gaussian elimiation GSL (very slow)
                                 // 'BICG' : Bi-conjugate gradient method with sparse routine
        }

        void check_Jacobian(double z, double eps)
        { check_jacobian=1; check_z=z; check_Jacobian_eps=eps; }

        void uncheck_Jacobian(){ check_jacobian=0;}

        void clear_Jacobian(){ Jacobian_ptr->clear(); }

        void clear();
        void set_verbosity(int verb) { verbosity=verb; }
    };

    //==============================================================================================
    // ODE_Solver_set_up_solution_and_memory::
    //
    // This routine has to be called before the solver is called. In particular the memory is set
    // and the function pointer for the ODE evaluation is given.
    //
    // ODE_solver_Solution Sz.z should be set to the initial redshift;
    // ODE_solver_Solution rtol[] & atol[] should contain the relative error request;
    // ODE_solver_Solution Sz.y[] should contain the initial solution;
    // ODE_solver_Solution Sz.dy[] will be computed using fn_ptr;
    //
    // ODE_solver_accuracies tols.rel[] should contain the relative error request;
    // ODE_solver_accuracies tols.abs[] should contain the relative error request;
    //
    // fn_ptr declares the ODE system according to y'= f(z, y);
    //
    //==============================================================================================
    void ODE_Solver_set_up_solution_and_memory(ODE_solver_Solution &Sz,
                                               ODE_solver_accuracies &tols,
                                               ODE_Solver_data &ODE_Solver_info,
                                               void (*fcn_ptr)(int *neq, double *z,
                                                               double *y, double *f,
                                                               int col),
                                               bool num_jac=1);

    void ODE_Solver_set_up_solution_and_memory(ODE_solver_Solution &Sz,
                                               ODE_solver_accuracies &tols,
                                               ODE_Solver_data &ODE_Solver_info,
                                               void (*fcn_ptr)(int *neq, double *z,
                                                               double *y, double *f,
                                                               int col, void *p),
                                               bool num_jac=1);

    //==============================================================================================
    // to set/reset error requirements
    //==============================================================================================
    void ODE_Solver_set_errors(const ODE_solver_accuracies &tols, ODE_Solver_data &ODE_Solver_info);

    //==============================================================================================
    // This routine has can only be called after the setup was done
    // (ODE_Solver_set_up_solution_and_memory)
    //
    // After the cal ODE_solver_Solution Sz will contain the solution at time zend; initially Sz
    // should contain the solution at zs. Previous solutions will also be stored by the routine.
    //
    //==============================================================================================
    int ODE_Solver_Solve_history(double zs, double zend,
                                 double Dz_in, double Dz_max,
                                 ODE_solver_Solution &Sz,
                                 ODE_Solver_data &ODE_Solver_info);

    //==============================================================================================
    // This routine has can only be called after the setup was done
    // (ODE_Solver_set_up_solution_and_memory)
    //
    // After the cal ODE_solver_Solution Sz will contain the solution at time zend; initially Sz
    // should contain the solution at zs. Previous solutions will also be stored by the routine.
    //
    // This routine is tailored for the problem of CMB anisotropies (Dodelson formulation)!!!
    //
    //==============================================================================================
    void ODE_Solver_set_up_solution_and_memory_Anisotropies(ODE_solver_Solution &Sz,
                                                            ODE_solver_accuracies &tols,
                                                            ODE_Solver_data &ODE_Solver_info,
                                                            void (*fcn_ptr)(int *neq, double *z,
                                                                            double *y, double *f,
                                                                            int col, void *p),
                                                            bool num_jac=1);

    int ODE_Solver_Solve_history_Anisotropies(double zs, double zend,
                                              double Dz_in, double Dz_max,
                                              ODE_solver_Solution &Sz,
                                              ODE_Solver_data &ODE_Solver_info, int n_Nu);

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
    // Oct 15th, 2013: added adaptive step-size control; Dz_in < Dz < Dz_max will be ensured;
    // to use this option make sure that NOT all relative errors are == 0
    // if this option is not required just set Dz_in==Dz_max (error control will have no effect)
    //==============================================================================================
    int ODE_Solver_Solve_history_PDE_ODE_OAE(double zs, double zend,
                                             double Dz_in, double Dz_max,
                                             ODE_solver_Solution &Sz,
                                             ODE_Solver_data &ODE_Solver_info,
                                             int nODE, int nOAE, int setmaxorder=5);

    //----------------------------------------------------------------------------------------------
    // version with overshooting allow until z == zcrit (if zcrit==zend --> no overshooting)
    //----------------------------------------------------------------------------------------------
    int ODE_Solver_Solve_history_PDE_ODE_OAE(double zs, double zend, double zcrit,
                                             double Dz_in, double Dz_max,
                                             ODE_solver_Solution &Sz,
                                             ODE_Solver_data &ODE_Solver_info,
                                             int nODE, int nOAE, int setmaxorder=5);

    //==============================================================================================
    // to obtain interpolated solution between zs and ze. The solver should have been run before
    // and the solution at ze should be available already.
    //==============================================================================================
    int ODE_Solver_Interpolate_Solution(double z,
                                        ODE_Solver_data &ODE_Solver_info,
                                        ODE_solver_Solution &Sz);
}

#endif
//==================================================================================================
//==================================================================================================
