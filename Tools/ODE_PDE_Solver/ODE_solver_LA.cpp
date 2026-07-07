//==================================================================================================
//
// main (sparse) matrix Linear Algebra routines for ODE-solver
//
// Author: Jens Chluba
//
// First implementation: 12/01/2010
// Last modification   : 22/03/2021
//==================================================================================================

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "ODE_solver_LA.h"
#include <gsl/gsl_linalg.h>
#include "routines.h"

using namespace std;

//=======================================================================================
namespace ODE_solver_LA
{
    //===================================================================================
    //
    // Linear-Algebra part
    //
    //===================================================================================

    //===================================================================================
    // All these routines are for simple linear algebra operations assuming that the
    // dimension of the problem is small (i.e. comparable to a few like in RECFAST)
    //===================================================================================

    //===================================================================================
    // scalar product c=a*b
    //===================================================================================
    double dot(const vector<double> &a, const vector<double> &b)
    {
        double scalar=0.0;
        for(unsigned int i=0; i<a.size(); i++) scalar+=a[i]*b[i];
        return scalar;
    }

    //===================================================================================
    // norm of vector a
    //===================================================================================
    double norm(const vector<double> &a){ return sqrt(dot(a, a)); }

    //===================================================================================
    // compute c=A*b (i.e. matrix times vector)
    //===================================================================================
    void A_times_b_is_c(const ODE_solver_matrix &Jac,
                        const vector<double> &b,
                        vector<double> &c)
    {
        //------------------------------------------------------
        // here it is assumed that A is symmetric with dimension
        // equal to b & c; c is overwritten
        //------------------------------------------------------
        for(unsigned int j=0; j<c.size(); j++) c[j]=0.0;

//        for(unsigned int i=0; i<Jac.A.size(); i++)
//            c[Jac.row[i]]+= Jac.A[i] * b[Jac.col[i]];

#ifdef OPENMP_ACTIVATED
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
        for(int r=0; r<Jac.dim; r++)
        {
            int n_entries=Jac.row_to_ind_col[r].size();
            for(int i=0; i<n_entries; i++)
                c[r]+=Jac.A[Jac.row_to_ind_col[r][i].index] * b[Jac.row_to_ind_col[r][i].col];
        }

        return;
    }

    //===================================================================================
    // get inverse diagonal elements of matrix A
    //===================================================================================
    void Get_inverse_diags(const ODE_solver_matrix &Jac, vector<double> &d)
    {
        //-----------------------------------------------------------
        // here it is assumed that A is symmetric with dimension dim
        //-----------------------------------------------------------
        for(int i=0; i<Jac.dim; i++)
        {
            d[i]=Jac.A[Jac.diags[i]];

            if(d[i]==0)
            {
                cerr << " error in preconditioner for element: " << i << endl;
                exit(0);
            }
            else d[i]=1.0/d[i];
        }

        return;
    }

    //===================================================================================
    //
    // ODE_solver_matrix routines
    //
    //===================================================================================
    void ODE_solver_matrix::clear()
    {
        dim=0;
        A.clear(); row.clear(); col.clear(); diags.clear();
        col_to_ind_row.clear(); row_to_ind_col.clear();
        col_to_ind_row_diag.clear(); row_to_ind_col_diag.clear();
    }

    void ODE_solver_matrix::save_info(int dimv, int c, int r, double Jij)
    {
        dim=dimv;
        A.push_back(Jij);
        col.push_back(c);
        row.push_back(r);

        if(c==r) diags.push_back((int)(A.size()-1));

        // save column and row info
        ODE_solver_el_info_ind_row cdum;
        ODE_solver_el_info_ind_col rdum;

        cdum.index=rdum.index=A.size()-1;
        cdum.row=r; rdum.col=c;

        if(col_to_ind_row.size()==0) col_to_ind_row.resize(dim);
        col_to_ind_row[c].push_back(cdum);

        if(row_to_ind_col.size()==0) row_to_ind_col.resize(dim);
        row_to_ind_col[r].push_back(rdum);

        // save diagional info
        if(c==r)
        {
            col_to_ind_row_diag.push_back(col_to_ind_row[c].size()-1);
            row_to_ind_col_diag.push_back(row_to_ind_col[r].size()-1);
        }
    }

    unsigned long ODE_solver_matrix::Get_element_index(int c, int r)
    {
        int n_entries=col_to_ind_row[c].size();

        for(int i=0; i<n_entries; i++)
            if(col_to_ind_row[c][i].row==r)
                return col_to_ind_row[c][i].index;

        return dim*dim+10;  // dim^2 - 1 is the maximal index that should appear! (Thanks Abir Sarkar!)
    }

    double ODE_solver_matrix::Get_element(int c, int r)
    {
        unsigned long j=Get_element_index(c, r);
        return (j==(unsigned long)dim*dim+10 ? 0.0 : A[j]);
    }

    void ODE_solver_matrix::show_entry(int c, int r)
    {
        cout << " ODE_solver_matrix::show_entry : col = " << c << " row= " << r << endl;

        unsigned long j=Get_element_index(c, r);

        if(j==(unsigned long)dim*dim+10) cout << " no entry (==0.0) " << endl;
        else cout << scientific << " Jij = " << A[j] << endl;
    }

    void ODE_solver_matrix::show_entries_col(int c)
    {
        int n_entries=col_to_ind_row[c].size();

        cout << " ODE_solver_matrix::show_entries_col : col = " << c << endl;

        for(int i=0; i<n_entries; i++)
        {
            col_to_ind_row[c][i].show_entry();
            cout << scientific << " Jij = " << A[col_to_ind_row[c][i].index] << endl;
        }

        cout << endl;
    }

    void ODE_solver_matrix::show_entries_row(int r)
    {
        int n_entries=row_to_ind_col[r].size();

        cout << " ODE_solver_matrix::show_entries_row : row = " << r << endl;

        for(int i=0; i<n_entries; i++)
        {
            row_to_ind_col[r][i].show_entry();
            cout << scientific << " Jij = " << A[row_to_ind_col[r][i].index] << endl;
        }

        cout << endl;
    }

    //===================================================================================
    // display routines
    //===================================================================================
    void ODE_solver_matrix::show_whole_matrix(vector<double > &b)
    {
        //return;
        cout.precision(3);
        cout << endl;
        for(int r=0; r<(int)diags.size(); r++)
        {
            cout << r << " : ";
            for(int c=0; c<(int)diags.size(); c++)
            {
                cout.width(4);
                cout << left << scientific << Get_element(c, r) << "  ";
            }

            if(b.size()>0) cout << "  ||  " << b[r] << endl << endl;
            else cout << endl;
        }

        //show_entry(1, 2);
        //show_entries_col(1);
        //wait_f_r();
        return;
    }

    void ODE_solver_matrix::show_whole_matrix()
    {
        vector<double> b;
        show_whole_matrix(b);
        return;
    }

    //===================================================================================
    // Iterative biconjugate gradiant routine -- BiCGSTAB
    //
    // BiCGSTAB solves the unsymmetric linear system Ax = b
    // using the Preconditioned BiConjugate Gradient Stabilized method
    //
    // BiCGSTAB follows the algorithm described on p. 27 of the
    // SIAM Templates book.
    //
    // The return value indicates convergence within max_iter (input)
    // iterations (0), or no convergence within max_iter iterations (1).
    //
    // Upon successful return, output arguments have the following values:
    //
    //        x  --  approximate solution to Ax = b
    // max_iter  --  the number of iterations performed before the
    //               tolerance was reached
    //      tol  --  the residual after the final iteration
    //
    // This routine was adapted from the IML++ http://math.nist.gov/iml++/
    //===================================================================================
    int BiCGSTAB_JC(const ODE_solver_matrix &Jac, const vector<double> &b,
                    vector<double> &x, int &max_iter, double &tol)
    {
        int neq=b.size();
        double resid;
        double rho_1=1.0, rho_2=1.0, alpha=1.0, beta=1.0, omega=1.0;
        //
        vector<double> p(neq), phat(neq), s(neq), shat(neq);
        vector<double> t(neq), v(neq), r(neq), rtilde(neq);
        vector<double> invdiag(neq);

        //------------------------------------------------------
        // this is to precondition the Matrix (i.e. A=M*Atilde)
        // here simply the inverse diagonal elements are used
        //------------------------------------------------------
        Get_inverse_diags(Jac, invdiag);

        double normb = norm(b);
        //------------------------------------------------------
        // r = b - A*x
        //------------------------------------------------------
        A_times_b_is_c(Jac, x, r);
        for(int i=0; i<neq; i++) r[i]=b[i]-r[i];
        //
        rtilde = r;

        if (normb == 0.0) normb = 1;

        if ((resid = norm(r) / normb) <= tol)
        {
            tol = resid;
            max_iter = 0;
            return 0;
        }

        for (int k = 1; k <= max_iter; k++)
        {
            rho_1 = dot(rtilde, r);
            if (rho_1 == 0)
            {
                tol = norm(r) / normb;
                max_iter = k;
                return 2;
            }

            if (k == 1) p = r;
            else
            {
                beta = (rho_1/rho_2) * (alpha/omega);
                for(int i=0; i<neq; i++) p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
            //------------------------------------------------------
            // preconditioning phat
            //------------------------------------------------------
            for(int i=0; i<neq; i++) phat[i]=invdiag[i]*p[i];

            //------------------------------------------------------
            // v = A * phat;
            //------------------------------------------------------
            A_times_b_is_c(Jac, phat, v);

            alpha = rho_1 / dot(rtilde, v);
            for(int i=0; i<neq; i++) s[i] = r[i] - alpha * v[i];

            if ((resid = norm(s)/normb) < tol)
            {
                for(int i=0; i<neq; i++) x[i] += alpha * phat[i];
                tol = resid;
                max_iter = k;
                return 0;
            }
            //------------------------------------------------------
            // preconditioning shat
            //------------------------------------------------------
            for(int i=0; i<neq; i++) shat[i]=invdiag[i]*s[i];

            //------------------------------------------------------
            // t = A * shat;
            //------------------------------------------------------
            A_times_b_is_c(Jac, shat, t);

            omega = dot(t,s) / dot(t,t);
            //
            for(int i=0; i<neq; i++) x[i] += alpha * phat[i] + omega * shat[i];
            for(int i=0; i<neq; i++) r[i] = s[i] - omega * t[i];

            rho_2 = rho_1;
            if ((resid = norm(r) / normb) < tol)
            {
                tol = resid;
                max_iter = k;
                return 0;
            }

            if (omega == 0)
            {
                tol = norm(r) / normb;
                max_iter = k;
                return 3;
            }
        }

        tol = resid;
        return 1;
    }

    int check_convergence(double &resid, int neq, double tol,
                          vector<double> &x, vector<double> &xi)
    {
        int iterate_again=-1;
        resid=0.0;

        for(int i=0; i<neq; i++)
        {
            double dx=abs((xi[i]-x[i])/(1.0e-300+x[i]));
            resid+=dx*dx;

            //cout << i << " " << dx << " " << xi[i] << endl;
            if(dx>tol){ iterate_again=1; break; }
        }
        //wait_f_r("check_convergence");
        resid=sqrt(resid)/neq;

        return iterate_again;
    }

    //===================================================================================
    // improved version with better error check
    //===================================================================================
    int BiCGSTAB_JC_II(const ODE_solver_matrix &Jac, const vector<double> &b,
                       vector<double> &x, int &max_iter, double &tol)
    {
        int neq=b.size();
        double normb = norm(b);

        // return for trivial equation A * x = 0
        if(normb==0.0){ for(int i=0; i<neq; i++) x[i]=0.0; max_iter=0; return 0; }

        double rho_1=1.0, rho_2=1.0, alpha=1.0, beta=1.0, omega=1.0;
        vector<double> p(neq), phat(neq), s(neq), shat(neq);
        vector<double> t(neq), v(neq), r(neq), rtilde(neq);
        vector<double> invdiag(neq), xold;

        double resid=0.0;
        int iterate_again_max=4, iterate_again=iterate_again_max;

        //------------------------------------------------------
        // this is to precondition the Matrix (i.e. A=M*Atilde)
        // here simply the inverse diagonal elements are used
        //------------------------------------------------------
        Get_inverse_diags(Jac, invdiag);

        //for(unsigned int i=0; i<x.size(); i++) cout << b[i] << " " << x[i] << " norm " << normb << endl;
        //wait_f_r();

        //------------------------------------------------------
        // r = b - A*x
        //------------------------------------------------------
        A_times_b_is_c(Jac, x, r);
        for(int i=0; i<neq; i++) r[i]=b[i]-r[i];
        rtilde = r;
        xold = x;

//        if ((resid = norm(r) / normb) <= tol)
//        {
//            tol = resid;
//            max_iter = 0;
//            return 0;
//        }

        for (int k = 1; k <= max_iter; k++)
        {
            rho_1 = dot(rtilde, r);
            if (rho_1 == 0) { tol = norm(r) / normb; return 2; }

            if (k == 1) p = r;
            else
            {
                beta = (rho_1/rho_2) * (alpha/omega);
                for(int i=0; i<neq; i++) p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
            //------------------------------------------------------
            // preconditioning phat
            //------------------------------------------------------
            for(int i=0; i<neq; i++) phat[i]=invdiag[i]*p[i];

            //------------------------------------------------------
            // v = A * phat;
            //------------------------------------------------------
            A_times_b_is_c(Jac, phat, v);

            alpha = rho_1 / dot(rtilde, v);
            for(int i=0; i<neq; i++) s[i] = r[i] - alpha * v[i];

            //------------------------------------------------------
            // preconditioning shat
            //------------------------------------------------------
            for(int i=0; i<neq; i++) shat[i]=invdiag[i]*s[i];

            //------------------------------------------------------
            // t = A * shat;
            //------------------------------------------------------
            A_times_b_is_c(Jac, shat, t);

            omega = dot(t,s) / dot(t,t);
            //
            for(int i=0; i<neq; i++) x[i] += alpha * phat[i] + omega * shat[i];
            for(int i=0; i<neq; i++) r[i] = s[i] - omega * t[i];

            rho_2 = rho_1;

            //------------------------------------------------------
            // check result (and make sure it is the same n times)
            //------------------------------------------------------
            iterate_again+=check_convergence(resid, neq, tol, xold, x);
            iterate_again=min(iterate_again_max, iterate_again);

            //------------------------------------------------------
            // return result if converged
            //------------------------------------------------------
            if(iterate_again==0)
            { tol = resid; max_iter = k; return 0; }

            if (omega == 0) { tol = norm(r) / normb; return 3; }

            xold = x;
        }

        tol = resid;
        return 1;
    }

    //===================================================================================
    // reimplemented version of stabilized BiCG without preconditioning (for now)
    //===================================================================================
    int BiCGSTAB_JC_III(const ODE_solver_matrix &Jac, const vector<double> &b,
                        vector<double> &x, int &max_iter, double &tol)
    {
        int neq=b.size();
        double normb = norm(b);

        // return for trivial equation A * x = 0
        if(normb==0.0){ for(int i=0; i<neq; i++) x[i]=0.0; max_iter=0; return 0; }

        double rho_i=0.0, rho_im1=1.0, omega_i=0.0, omega_im1=1.0;
        double alpha=1.0, beta=0.0;
        vector<double> pi(neq, 0.0), nui(neq, 0.0), xi(neq, 0.0);
        vector<double> ri(neq, 0.0), rhat(neq, 0.0), h(neq, 0.0), s(neq, 0.0), t(neq, 0.0);

        double resid=0.0;
        int iterate_again_max=4, iterate_again=iterate_again_max;

        //------------------------------------------------------
        // r = b - A*x [i=0]
        //------------------------------------------------------
        A_times_b_is_c(Jac, x, ri);
        for(int i=0; i<neq; i++) ri[i]=b[i]-ri[i];
        rhat=ri;

        for (int k = 1; k <= max_iter; k++)
        {
            rho_i = dot(rhat, ri);
            beta = (rho_i / rho_im1) * ( alpha / omega_im1);

            // update pim1 --> pi
            for(int i=0; i<neq; i++) pi[i] = ri[i] + beta * (pi[i] - omega_im1 * nui[i]);

            // update nuim1 --> nui
            A_times_b_is_c(Jac, pi, nui);
            alpha=rho_i/dot(rhat, nui);

            for(int i=0; i<neq; i++) h[i] = x[i] + alpha * pi[i];

            //------------------------------------------------------
            // check result
            //------------------------------------------------------
            //iterate_again=check_convergence(resid, neq, tol, x, h);

            //------------------------------------------------------
            // return result if converged
            //------------------------------------------------------
            //if(!iterate_again)
            //{ x=h; tol = resid; max_iter = k; return 0; }

            //------------------------------------------------------
            // otherwise continue with other part
            //------------------------------------------------------
            for(int i=0; i<neq; i++) s[i] = ri[i] - alpha * nui[i];
            A_times_b_is_c(Jac, s, t);

            omega_i=dot(t, s)/dot(t, t);

            for(int i=0; i<neq; i++) xi[i] = h[i] + omega_i * s[i];

            //------------------------------------------------------
            // check result (and make sure it is the same n times)
            //------------------------------------------------------
            iterate_again+=check_convergence(resid, neq, tol, x, xi);
            iterate_again=min(iterate_again_max, iterate_again);

            //------------------------------------------------------
            // return result if converged
            //------------------------------------------------------
            if(iterate_again==0)
            { x=xi; tol = resid; max_iter = k; return 0; }

            // update rim1 --> ri
            for(int i=0; i<neq; i++) ri[i] = s[i] - omega_i * t[i];

            // switch variables
            omega_im1=omega_i;
            rho_im1=rho_i;
            x=xi;
        }

        return 1;
    }

    //===================================================================================
    // Solving the Matrix Equation A*x == b for x
    //===================================================================================
    int ODE_Solver_Solve_LA_system(ODE_solver_matrix &Jac,
                                   vector<double > &b,
                                   vector<double > &x,
                                   double tol,
                                   int verbose)
    {
        int ifail=0, maxit = 50000;                       // Maximum iterations

        //for(unsigned int i=0; i<x.size(); i++) x[i]=0.0; // reset x-vector

        //ifail = BiCGSTAB_JC(Jac, b, x, maxit, tol);      // initial version of Solve linear system
        ifail = BiCGSTAB_JC_II(Jac, b, x, maxit, tol);      // has better convergence
        //ifail = BiCGSTAB_JC_III(Jac, b, x, maxit, tol);      // Solve linear system

        if(verbose>0 && ifail!=0)
        {
            cout << endl << endl;
            cout << " BiCG flag = " << ifail << endl;
            cout << " iterations performed: " << maxit << endl;
            cout << " tolerance achieved  : " << tol << endl;
            cout << " Something bad happened to the ODE system. Restarting solver. "
                 << endl;

            if(isnan_JC(tol)) ifail=10;
        }

        if(verbose>0) cout << " average tolerance achieved  : " << tol
                           << " with " << maxit << " iterations" << endl;

        return ifail;
    }

    //===================================================================================
    // Gaussian elimination using matrix routines [ solves A * x == b ]
    //===================================================================================
    int ODE_Solver_Solve_LA_system_Gaussian(ODE_solver_matrix &M,
                                            vector<double > &b,
                                            vector<double > &x,
                                            int verbose)
    {
        if(verbose>1) M.show_whole_matrix(b);

        //===============================================================================
        // forward elimiation row by row
        //===============================================================================
        for(int r=0; r<M.dim; r++)
        {
            int n_row=M.col_to_ind_row[r].size();
            int n_col=M.row_to_ind_col[r].size();

            // pivot element of row
            double pivot=M.A[M.diags[r]];

            // divide pivot-row by pivot at c>=r
            for(int j=M.row_to_ind_col_diag[r]; j<n_col; j++)
            {
                long int index=M.row_to_ind_col[r][j].index;
                M.A[index]/=pivot;
            }

            // devide vector by pivot
            b[r]/=pivot;

            // eliminate sub-diagonals
            for(int subrow=M.col_to_ind_row_diag[r]+1; subrow<n_row; subrow++)
            {
                int piv_ind=M.col_to_ind_row[r][subrow].index;

                if(M.A[piv_ind]==0.0) continue; // save time

                int piv_row=M.col_to_ind_row[r][subrow].row;

                for(int j=M.row_to_ind_col_diag[r]+1; j<n_col; j++)
                {
                    // get col index
                    long int index=M.row_to_ind_col[r][j].index;

                    if(M.A[index]==0.0) continue; // save time

                    int c=M.row_to_ind_col[r][j].col;

                    // eliminate sub-diagonal elements
                    int n_eli_row=M.col_to_ind_row[c].size();
                    // conservative estimate of start index
                    // [TODO: may need checking for general cases... JC]
                    int eli_row_start=max(0, M.col_to_ind_row_diag[c]-(c-r));

                    for(int eli_row=eli_row_start; eli_row<n_eli_row; eli_row++)
                        if(M.col_to_ind_row[c][eli_row].row==piv_row)
                        {
                            long int eli_ind=M.col_to_ind_row[c][eli_row].index;
                            M.A[eli_ind]-=M.A[piv_ind]*M.A[index];
                            break; // save time
                        }
                }

                // subtract from target
                b[piv_row]-=M.A[piv_ind]*b[r];

                // eliminate sub-diagonal element of pivot column
                M.A[piv_ind]=0.0;
            }
        }

        if(verbose>1) M.show_whole_matrix(b);

        //===============================================================================
        // back-substitution
        //===============================================================================
        for(int k=0; k<M.dim; k++) x[k]=b[k];

        for(int c=M.dim-1; c>0; c--)
            for(int subrow=0; subrow<M.col_to_ind_row_diag[c]; subrow++)
            {
                int piv_ind=M.col_to_ind_row[c][subrow].index;
                int piv_row=M.col_to_ind_row[c][subrow].row;

                x[piv_row]-=M.A[piv_ind]*x[c];
            }

        if(verbose>2)
        {
            for(int k=0; k<M.dim; k++) cout << k << " " << x[k] << endl;
            wait_f_r(" Gauss-JC solver ");
        }

        return 0;
    }

    //===================================================================================
    // Solving the Matrix Equation A*x == b for x using GSL
    //===================================================================================
    int ODE_Solver_Solve_LA_system_GSL(ODE_solver_matrix &M,
                                       vector<double > &bi,
                                       vector<double > &x,
                                       int verbose)
    {
        //cout << " Entering GSL " << endl;

        int npxi=bi.size();

        if(npxi>10 && verbose>1) M.show_whole_matrix(bi);

        vector<double> a_data(npxi*npxi, 0.0);

        for(int r=0; r<(int)M.col.size(); r++)
        {
            int ind=M.col[r]+npxi*M.row[r];
            a_data[ind]=M.A[r];
        }

        gsl_matrix_view m=gsl_matrix_view_array (&a_data[0], npxi, npxi);
        gsl_vector_view b=gsl_vector_view_array (&bi[0], npxi);
        gsl_vector *x_GSL = gsl_vector_alloc (npxi);

        int s;
        gsl_permutation * p = gsl_permutation_alloc (npxi);
        //cout << " LU decomposition " << endl;
        gsl_linalg_LU_decomp (&m.matrix, p, &s);
        //cout << " LU solve " << endl;
        gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x_GSL);

        for(int k=0; k<npxi; k++) x[k]=x_GSL->data[k];

        gsl_permutation_free (p);
        gsl_vector_free (x_GSL);
        a_data.clear();

        if(npxi>10 && verbose>2)
        {
            for(int k=0; k<npxi; k++) cout << k << " " << x[k] << endl;
            wait_f_r();
        }

        return 0;
    }

    //===================================================================================
    void eliminate_n_lower_plus_phot(int nlow, double fac_el, int i_gamma, int col, int npxi,
                                     int n_elim, vector<int> &elim_ind,
                                     vector<double > &A, vector<double > &bi)
    {
        bi[col]/=fac_el;
        for(int k=1; k<=nlow; k++) bi[col+k]-=A[col+npxi*(col+k)]*bi[col];

        bi[i_gamma+0]-=A[col+npxi*(i_gamma+0)]*bi[col];
        bi[i_gamma+1]-=A[col+npxi*(i_gamma+1)]*bi[col];

        for(int i=n_elim-1; i>=0; i--)
        {
            A[elim_ind[i]+npxi*col]/=fac_el;
            for(int k=1; k<=nlow; k++)
                A[elim_ind[i]+npxi*(col+k)]-=A[col+npxi*(col+k)]*A[elim_ind[i]+npxi*col];

            A[elim_ind[i]+npxi*(i_gamma+0)]-=A[col+npxi*(i_gamma+0)]*A[elim_ind[i]+npxi*col];
            A[elim_ind[i]+npxi*(i_gamma+1)]-=A[col+npxi*(i_gamma+1)]*A[elim_ind[i]+npxi*col];
        }

        return;
    }

    //===================================================================================
    void eliminate_n_lower(int nlow, double fac_el, int col, int npxi,
                           int n_elim, vector<int> &elim_ind,
                           vector<double > &A, vector<double > &bi)
    {
        bi[col]/=fac_el;
        for(int k=1; k<=nlow; k++) bi[col+k]-=A[col+npxi*(col+k)]*bi[col];

        for(int i=n_elim-1; i>=0; i--)
        {
            A[elim_ind[i]+npxi*col]/=fac_el;
            for(int k=1; k<=nlow; k++)
                A[elim_ind[i]+npxi*(col+k)]-=A[col+npxi*(col+k)]*A[elim_ind[i]+npxi*col];
        }

        return;
    }

    //===================================================================================
    void eliminate_n_lower_plus_pol(int nlow, double fac_el, int i_gamma_P, int col, int npxi,
                                     int n_elim, vector<int> &elim_ind,
                                     vector<double > &A, vector<double > &bi)
    {
        bi[col]/=fac_el;
        for(int k=1; k<=nlow; k++) bi[col+k]-=A[col+npxi*(col+k)]*bi[col];

        bi[i_gamma_P+0]-=A[col+npxi*(i_gamma_P+0)]*bi[col];
        bi[i_gamma_P+2]-=A[col+npxi*(i_gamma_P+2)]*bi[col];

        for(int i=n_elim-1; i>=0; i--)
        {
            A[elim_ind[i]+npxi*col]/=fac_el;
            for(int k=1; k<=nlow; k++)
                A[elim_ind[i]+npxi*(col+k)]-=A[col+npxi*(col+k)]*A[elim_ind[i]+npxi*col];

            A[elim_ind[i]+npxi*(i_gamma_P+0)]-=A[col+npxi*(i_gamma_P+0)]*A[elim_ind[i]+npxi*col];
            A[elim_ind[i]+npxi*(i_gamma_P+2)]-=A[col+npxi*(i_gamma_P+2)]*A[elim_ind[i]+npxi*col];
        }

        return;
    }

    //===================================================================================
    // Solving the Matrix Equation A*x == b for x for the cosmological perturbation
    // equations. The sparesness and band structure are used for the high l-modes.
    //===================================================================================
    int ODE_Solver_Solve_LA_system_Anisotropies(ODE_solver_matrix &M,
                                                vector<double > &bi,
                                                vector<double > &x,
                                                int n_Nu,
                                                int verbose)
    {
        int npxi=bi.size();
        int i_gamma=n_Nu+5;
        int n_gamma=(npxi-i_gamma)/2;
        int i_gamma_P=i_gamma+n_gamma;
        int n_elim;
        vector<int> elim_ind(npxi);

        //===================================================================================
        // eliminate first column
        //===================================================================================
        int col=0;

        double fac_el=M.A[col+npxi*col];

//        A.show_whole_matrix(bi);

        elim_ind[0]=col;
        elim_ind[1]=col+1;
        elim_ind[2]=col+3;
        elim_ind[3]=col+5;
        elim_ind[4]=col+7;
        elim_ind[5]=i_gamma+0;
        elim_ind[6]=i_gamma+2;
        n_elim=7;

        eliminate_n_lower_plus_phot(6, fac_el, i_gamma, col, npxi, n_elim, elim_ind, M.A, bi);

        //===================================================================================
        // eliminate second column
        //===================================================================================
        col=1;

        fac_el=M.A[col+npxi*col];
        elim_ind[0]=col;
        elim_ind[1]=col+1;
        elim_ind[2]=col+2;
        elim_ind[3]=col+4;
        elim_ind[4]=col+6;
        elim_ind[5]=i_gamma+0;
        elim_ind[6]=i_gamma+2;
        n_elim=7;

        eliminate_n_lower_plus_phot(5, fac_el, i_gamma, col, npxi, n_elim, elim_ind, M.A, bi);

        //===================================================================================
        // eliminate third column
        //===================================================================================
        col=2;

        fac_el=M.A[col+npxi*col];
        elim_ind[0]=col;
        elim_ind[1]=col+1;
        elim_ind[2]=col+3;
        elim_ind[3]=col+5;
        elim_ind[4]=i_gamma+0;
        elim_ind[5]=i_gamma+2;
        n_elim=6;

        eliminate_n_lower_plus_phot(4, fac_el, i_gamma, col, npxi, n_elim, elim_ind, M.A, bi);

        //===================================================================================
        // eliminate fourth column
        //===================================================================================
        col=3;

        fac_el=M.A[col+npxi*col];
        elim_ind[0]=col;
        elim_ind[1]=col+1;
        elim_ind[2]=col+2;
        elim_ind[3]=col+4;
        elim_ind[4]=i_gamma+0;
        elim_ind[5]=i_gamma+2;
        n_elim=6;

        eliminate_n_lower_plus_phot(3, fac_el, i_gamma, col, npxi, n_elim, elim_ind, M.A, bi);

        //===================================================================================
        // eliminate fifth column
        //===================================================================================
        col=4;

        fac_el=M.A[col+npxi*col];
        elim_ind[0]=col;
        elim_ind[1]=col+1;
        elim_ind[2]=col+3;
        elim_ind[3]=i_gamma+0;
        elim_ind[4]=i_gamma+1;
        elim_ind[5]=i_gamma+2;
        n_elim=6;

        eliminate_n_lower_plus_phot(2, fac_el, i_gamma, col, npxi, n_elim, elim_ind, M.A, bi);

        //===================================================================================
        // eliminate first neutrino column
        //===================================================================================
        col=5;

        fac_el=M.A[col+npxi*col];
        elim_ind[0]=col;
        elim_ind[1]=col+1;
        elim_ind[2]=col+2;
        elim_ind[3]=i_gamma+0;
        elim_ind[4]=i_gamma+1;
        elim_ind[5]=i_gamma+2;
        n_elim=6;

        eliminate_n_lower_plus_phot(1, fac_el, i_gamma, col, npxi, n_elim, elim_ind, M.A, bi);

        //===================================================================================
        // eliminate remaining neutrino columns
        //===================================================================================
        for(col=6; col<i_gamma-1; col++)
        {
            fac_el=M.A[col+npxi*col];
            elim_ind[0]=col;
            elim_ind[1]=col+1;
            elim_ind[2]=i_gamma+0;
            elim_ind[3]=i_gamma+1;
            elim_ind[4]=i_gamma+2;
            n_elim=5;

            eliminate_n_lower_plus_phot(1, fac_el, i_gamma, col, npxi, n_elim, elim_ind, M.A, bi);
        }

        //===================================================================================
        // eliminate last neutrino column
        //===================================================================================
        col=i_gamma-1;

        fac_el=M.A[col+npxi*col];
        elim_ind[0]=col;
        elim_ind[1]=i_gamma+0;
        elim_ind[2]=i_gamma+1;
        elim_ind[3]=i_gamma+2;
        n_elim=4;

        eliminate_n_lower(2, fac_el, col, npxi, n_elim, elim_ind, M.A, bi);

        //===================================================================================
        // eliminate first photon column
        //===================================================================================
        col=i_gamma;

        fac_el=M.A[col+npxi*col];
        elim_ind[0]=col;
        elim_ind[1]=col+1;
        elim_ind[2]=col+2;
        n_elim=3;

        eliminate_n_lower(1, fac_el, col, npxi, n_elim, elim_ind, M.A, bi);

        //===================================================================================
        // eliminate second photon column
        //===================================================================================
        col=i_gamma+1;

        fac_el=M.A[col+npxi*col];
        elim_ind[0]=col;
        elim_ind[1]=col+1;
        n_elim=2;

        eliminate_n_lower(1, fac_el, col, npxi, n_elim, elim_ind, M.A, bi);

        //===================================================================================
        // eliminate remaining photon columns
        //===================================================================================
        for(col=i_gamma+2; col<i_gamma_P-1; col++)
        {
            fac_el=M.A[col+npxi*col];
            elim_ind[0]=col;
            elim_ind[1]=col+1;
            elim_ind[2]=i_gamma_P+0;
            elim_ind[3]=i_gamma_P+2;
            n_elim=4;

            eliminate_n_lower_plus_pol(1, fac_el, i_gamma_P, col, npxi, n_elim, elim_ind, M.A, bi);
        }

        //===================================================================================
        // eliminate last photon column
        //===================================================================================
        col=i_gamma_P-1;

        fac_el=M.A[col+npxi*col];
        elim_ind[0]=col;
        elim_ind[1]=i_gamma_P+0;
        elim_ind[2]=i_gamma_P+2;
        n_elim=3;

        eliminate_n_lower_plus_pol(0, fac_el, i_gamma_P, col, npxi, n_elim, elim_ind, M.A, bi);

        //===================================================================================
        // eliminate first polarization  column
        //===================================================================================
        col=i_gamma_P;

        fac_el=M.A[col+npxi*col];
        elim_ind[0]=col;
        elim_ind[1]=col+1;
        elim_ind[2]=col+2;
        n_elim=3;

        eliminate_n_lower(2, fac_el, col, npxi, n_elim, elim_ind, M.A, bi);

        //===================================================================================
        // eliminate remaining polarization columns
        //===================================================================================
        for(col=i_gamma_P+1; col<npxi-1; col++)
        {
            fac_el=M.A[col+npxi*col];
            elim_ind[0]=col;
            elim_ind[1]=col+1;
            n_elim=2;

            eliminate_n_lower(1, fac_el, col, npxi, n_elim, elim_ind, M.A, bi);
        }


        fac_el=M.A[col+npxi*col];
        bi[col]/=fac_el;
        M.A[col+npxi*col]/=fac_el;

//        A.show_whole_matrix(bi);

        //===================================================================================
        // solve
        //===================================================================================
        x[npxi-1]=bi[npxi-1];

        // polarization
        for(col=npxi-2; col>i_gamma_P; col--)
            x[col]=bi[col]-x[col+1]*M.A[col+1+npxi*col];

        // polarization monopole
        x[col]=bi[col]-x[col+1]*M.A[col+1+npxi*col]-x[col+2]*M.A[col+2+npxi*col];

        // photons (last col)
        col--;
        x[col]=bi[col];
        for(int i=0; i<3; i++) x[col]+=-x[i_gamma_P+i]*M.A[i_gamma_P+i+npxi*col];

        // photons (rest)
        for(col--; col>i_gamma; col--)
        {
            x[col]=bi[col]-x[col+1]*M.A[col+1+npxi*col];
            for(int i=0; i<3; i++)
                x[col]+=-x[i_gamma_P+i]*M.A[i_gamma_P+i+npxi*col];
        }

        // photon monopole
        x[col]=bi[col]-x[col+1]*M.A[col+1+npxi*col]-x[col+2]*M.A[col+2+npxi*col];

        // neutrinos
        col--;
        x[col]=bi[col];
        for(int i=0; i<3; i++) x[col]+=-x[i_gamma+i]*M.A[i_gamma+i+npxi*col];

        for(col--; col>5; col--)
        {
            x[col]=bi[col]-x[col+1]*M.A[col+1+npxi*col];
            for(int i=0; i<3; i++) x[col]+=-x[i_gamma+i]*M.A[i_gamma+i+npxi*col];
        }

        // neutrino monopole
        x[col]=bi[col]-x[col+1]*M.A[col+1+npxi*col]-x[col+2]*M.A[col+2+npxi*col];
        for(int i=0; i<3; i++) x[col]+=-x[i_gamma+i]*M.A[i_gamma+i+npxi*col];

        // baryon velocity
        col--;
        x[col]=bi[col]-x[col+1]*M.A[col+1+npxi*col]-x[col+3]*M.A[col+3+npxi*col];
        for(int i=0; i<3; i++) x[col]+=-x[i_gamma+i]*M.A[i_gamma+i+npxi*col];

        // baryon density
        col--;
        x[col] =bi[col]-x[col+1]*M.A[col+1+npxi*col]-x[col+2]*M.A[col+2+npxi*col];
        x[col]-=x[col+4]*M.A[col+4+npxi*col];
        x[col]+=-x[i_gamma+0]*M.A[i_gamma+0+npxi*col]-x[i_gamma+2]*M.A[i_gamma+2+npxi*col];

        // DM velocity
        col--;
        x[col] =bi[col]-x[col+1]*M.A[col+1+npxi*col]-x[col+3]*M.A[col+3+npxi*col];
        x[col]-=x[col+5]*M.A[col+5+npxi*col];
        x[col]+=-x[i_gamma+0]*M.A[i_gamma+0+npxi*col]-x[i_gamma+2]*M.A[i_gamma+2+npxi*col];

        // DM density
        col--;
        x[col] =bi[col]-x[col+1]*M.A[col+1+npxi*col]-x[col+2]*M.A[col+2+npxi*col];
        x[col]-=x[col+4]*M.A[col+4+npxi*col];
        x[col]-=x[col+6]*M.A[col+6+npxi*col];
        x[col]+=-x[i_gamma+0]*M.A[i_gamma+0+npxi*col]-x[i_gamma+2]*M.A[i_gamma+2+npxi*col];

        // phi
        col--;
        x[col] =bi[col]-x[col+1]*M.A[col+1+npxi*col]-x[col+3]*M.A[col+3+npxi*col];
        x[col]-=x[col+5]*M.A[col+5+npxi*col];
        x[col]-=x[col+7]*M.A[col+7+npxi*col];
        x[col]+=-x[i_gamma+0]*M.A[i_gamma+0+npxi*col]-x[i_gamma+2]*M.A[i_gamma+2+npxi*col];

        elim_ind.clear();

//        for(int k=0; k<(int)xc.size(); k++) cout << k << " " << x[k] << " " << xc[k] << endl;
//        wait_f_r();

        return 0;
    }
}

//==================================================================================================
//==================================================================================================
