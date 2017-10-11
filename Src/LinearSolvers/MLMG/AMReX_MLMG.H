#ifndef AMREX_ML_MG_H_
#define AMREX_ML_MG_H_

#include <AMReX_MLLinOp.H>

namespace amrex {

class MLMG
{
public:

    using BCMode = MLLinOp::BCMode;

    MLMG (MLLinOp& a_lp);
    ~MLMG ();

    void solve (const Vector<MultiFab*>& a_sol, const Vector<MultiFab const*>& a_rhs,
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

    // To avoid confusion, terms like sol, cor, rhs, res, ... etc. are
    // in the frame of the original equation, not the correction form

    Vector<MultiFab*>         sol;  // alias to argument a_sol
    Vector<MultiFab>          rhs;  // Copy of original rhs

    // First Vector: Amr levels.  0 is the coarest level
    // Second Vector: MG levels.  0 is the finest level
    Vector<Vector<MultiFab> > res;    // = rhs - L(sol)
    Vector<Vector<MultiFab> > cor;    // L(cor) = res
    Vector<Vector<MultiFab> > rescor; // = res - L(cor)   Residual of the correction form

    Vector<MultiFab> cor_hold;

    void oneIter ();
    void miniCycle (int alev);
    void mgCycle ();
    void computeResidual (int alev, int mlev);
    void computeResWithCrseSolFineCor (int crse_amr_lev, int fine_amr_lev);
    void computeResWithCrseCorFineCor (int crse_amr_lev, int fine_amr_lev);
    void interpCorrection (int alev);
};

}

#endif
