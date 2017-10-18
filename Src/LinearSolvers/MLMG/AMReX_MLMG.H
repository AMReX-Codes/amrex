#ifndef AMREX_ML_MG_H_
#define AMREX_ML_MG_H_

#include <AMReX_MLLinOp.H>
#include <AMReX_iMultiFab.H>

namespace amrex {

class MLMG
{
public:

    using BCMode = MLLinOp::BCMode;

    enum class Cycle { Vcycle, Wcycle, Fcycle };

    enum class BottomSolver : int { smoother, bicgstab, cg };

    MLMG (MLLinOp& a_lp);
    ~MLMG ();

    void solve (const Vector<MultiFab*>& a_sol, const Vector<MultiFab const*>& a_rhs,
                Real a_tol_real, Real a_tol_abs);

    void setVerbose (int v) { verbose = v; }
    void setMaxIter (int n) { max_iters = n; }
    void setFixedIter (int nit) { do_fixed_number_of_iters = nit; }

    void setPreSmooth (int n) { nu1 = n; }
    void setPostSmooth (int n) { nu2 = n; }
    void setCtrRelax (int n) { nu0 = n; }
    void setFinalSmooth (int n) { nuf = n; }
    void setBottomSmooth (int n) { nub = n; }

    void setCycle (Cycle c) { cycle_type = c; }

    void setBottomSolver (BottomSolver s) { bottom_solver = s; }
    void setCGVerbose (int v) { cg_verbose = v; }
    void setCGMaxIter (int n) { cg_maxiter = n; }
    void setCGUnstableCriterion (Real x) { cg_unstable_criterion = x; }

    void setAlwaysUseBNorm (int flag) { always_use_bnorm = flag; }

private:

    int verbose = 1;
    int max_iters = 200;
    int do_fixed_number_of_iters = 0;

    int nu0 = 1;   // center
    int nu1 = 2;   // pre
    int nu2 = 2;   // post
    int nuf = 8;   // when smoother is used as bottom solver
    int nub = 0;   // aditional smoothing after bottom cg solver
    
    Cycle cycle_type = Cycle::Vcycle;

    BottomSolver bottom_solver = BottomSolver::bicgstab;
    int  cg_verbose            = 0;
    int  cg_maxiter            = 100;
    Real cg_unstable_criterion = 10.0;

    int always_use_bnorm = 0;

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

    Vector<iMultiFab> fine_mask;

    void oneIter (int iter);

    void miniCycle (int alev);
    void mgCycle ();

    void mgVcycle (int amrlev, int mglev);

    void bottomSolve ();

    void computeMLResidual (int amrlevmax, bool bndryregister_updated = false);
    void computeResidual (int alev);
    void computeResWithCrseSolFineCor (int crse_amr_lev, int fine_amr_lev);
    void computeResWithCrseCorFineCor (int fine_amr_lev);
    void interpCorrection (int alev);
    void interpCorrection (int alev, int mglev);

    void computeResOfCorrection (int amrlev, int mglev);

    Real ResNormInf (int amrlev, bool local = false);
    Real MLResNormInf (int alevmax, bool local = false);
    Real MLRhsNormInf (bool local = false);
    void buildFineMask ();
};

}

#endif
