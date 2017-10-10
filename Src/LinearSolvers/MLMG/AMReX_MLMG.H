#ifndef AMREX_ML_MG_H_
#define AMREX_ML_MG_H_

#include <AMReX_MLLinOp.H>

namespace amrex {

class MLMG
{
public:

    MLMG (MLLinOp& a_lp);
    ~MLMG ();

    void solve (const Vector<MultiFab*>& sol, const Vector<MultiFab const*>& a_rhs,
                Real a_tol_real, Real a_tol_abs);

    void setVerbose (int v) { verbose = v; }
    void setMaxIters (int i) { max_iters = i; }

private:

    int verbose = 1;
    int max_iters = 200;

    int nu1 = 1;
    

    MLLinOp& linop;
    int namrlevs;
    int finest_amr_lev;

                                      // Argument sol is the full solution
    Vector<MultiFab>          rhs_o;  // Original rhs

    // First Vector: Amr levels.  0 is the coarest level
    // Second Vector: MG levels.  0 is the finest level
    Vector<Vector<MultiFab> > rhs_c;  // = rhs - L(soln)    RHS of the correction form.
    Vector<Vector<MultiFab> > cor;    // L(cor) = rhs_c     The correction
    Vector<Vector<MultiFab> > res;    // = rhs_c - L(cor)   Residual of the correction form

    void oneIter (const Vector<MultiFab*>& sol);
    void miniCycle (int alev);
    void mgCycle ();
    void computeResidual (int alev);
};

}

#endif
