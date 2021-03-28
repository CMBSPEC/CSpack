//==================================================================================================
//
// main (sparse) matrix Linear Algebra routines for ODE-solver
//
// Author: Jens Chluba
//
// First implementation: 12/01/2010
// Last modification   : 22/03/2021
//==================================================================================================

#ifndef _ODE_SOLVER_LA_H_
#define _ODE_SOLVER_LA_H_

#include <vector>
#include <string>

using namespace std;

namespace ODE_solver_LA
{
    //==============================================================================================
    struct ODE_solver_el_info_ind_col
    {
        unsigned long index;
        int col;
        
        void show_entry()
        {
            cout << " ODE_solver_el_info_ind_col::show_entry :"
                 << " index = " << index << " col = " << col << " ";
        }
    };

    struct ODE_solver_el_info_ind_row
    {
        unsigned long index;
        int row;
        
        void show_entry()
        {
            cout << " ODE_solver_el_info_ind_row::show_entry :"
                 << " index = " << index << " row = " << row << " ";
        }
    };

    struct ODE_solver_matrix
    {
        int dim;
        vector<double> A;
        vector<int> row;    // index --> row
        vector<int> col;    // index --> col
        vector<int> diags;  // diag  --> index
        
        // additional ways to slice matrix
        vector<int> col_to_ind_row_diag;
        vector<int> row_to_ind_col_diag;
        vector<vector<ODE_solver_el_info_ind_row> > col_to_ind_row;  // col --> <index, row>
        vector<vector<ODE_solver_el_info_ind_col> > row_to_ind_col;  // row --> <index, col>
        
        ODE_solver_matrix(){ dim=0; }
        
        void clear();
        void save_info(int dimv, int c, int r, double Jij);
        
        unsigned long Get_element_index(int c, int r);
        double Get_element(int c, int r);

        void show_entry(int c, int r);
        void show_entries_col(int c);
        void show_entries_row(int c);
        void show_whole_matrix();
        void show_whole_matrix(vector<double > &b);
    };

    //==============================================================================================
    // some useful matrix routines for solving A x = b and multiplying A x to obtain c
    //==============================================================================================
    void A_times_b_is_c(const ODE_solver_matrix &Jac,
                        const vector<double> &b,
                        vector<double> &c);

    //===================================================================================
    // a dot b --> c
    //===================================================================================
    double dot(const vector<double> &a, const vector<double> &b);

    //===================================================================================
    // Bi-conjugate gradient method for iterative solve
    //===================================================================================
    int ODE_Solver_Solve_LA_system(ODE_solver_matrix &Jac,
                                   vector<double > &b,
                                   vector<double > &x,
                                   double tol,
                                   int verbose=0);

    //===================================================================================
    // sparse matrix Gaussian elimination
    //===================================================================================
    int ODE_Solver_Solve_LA_system_Gaussian(ODE_solver_matrix &M,
                                            vector<double > &b,
                                            vector<double > &x,
                                            int verbose=0);

    //===================================================================================
    // slowest solver using GSL Gaussian elimination
    //===================================================================================
    int ODE_Solver_Solve_LA_system_GSL(ODE_solver_matrix &M,
                                       vector<double > &bi,
                                       vector<double > &x,
                                       int verbose=0);

    //===================================================================================
    // Solving the Matrix Equation A*x == b for x for the cosmological perturbation
    // equations. The sparesness and band structure are used for the high l-modes.
    // COMMENT: call of 'ODE_Solver_Solve_LA_system_Gaussian' is probably equivalent
    //===================================================================================
    int ODE_Solver_Solve_LA_system_Anisotropies(ODE_solver_matrix &M,
                                                vector<double > &bi,
                                                vector<double > &x,
                                                int n_Nu,
                                                int verbose=0);
}

#endif
//==================================================================================================
//==================================================================================================
