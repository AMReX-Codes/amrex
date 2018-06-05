#ifndef AMREX_ML_MG_H_
#define AMREX_ML_MG_H_

#include <AMReX_MLLinOp.H>
#include <AMReX_iMultiFab.H>

#ifdef AMREX_USE_HYPRE
#include <AMReX_HypreABecLap2.H>
#endif

namespace amrex {

class MLMG
{
public:

    using BCMode = MLLinOp::BCMode;

    enum class BottomSolver : int { smoother, bicgstab, cg, hypre };

    MLMG (MLLinOp& a_lp);
    ~MLMG ();

    Real solve (const Vector<MultiFab*>& a_sol, const Vector<MultiFab const*>& a_rhs,
                Real a_tol_rel, Real a_tol_abs);

    void getGradSolution (const Vector<std::array<MultiFab*,AMREX_SPACEDIM> >& a_grad_sol);
    void getFluxes (const Vector<std::array<MultiFab*,AMREX_SPACEDIM> >& a_grad_sol);
    void compResidual (const Vector<MultiFab*>& a_res, const Vector<MultiFab*>& a_sol,
                       const Vector<MultiFab const*>& a_rhs);
    // out = L(in). Note that, if no acutal solve is needed, one could
    // turn off multigrid coarsening by constructing a MLLinOp object
    // with an appropriate LPInfo object (e.g., with LPInfo().setMaxCoarseningLevel(0)).
    void apply (const Vector<MultiFab*>& out, const Vector<MultiFab*>& in);

    void setVerbose (int v) { verbose = v; }
    void setMaxIter (int n) { max_iters = n; }
    void setMaxFmgIter (int n) { max_fmg_iters = n; }
    void setFixedIter (int nit) { do_fixed_number_of_iters = nit; }

    void setPreSmooth (int n) { nu1 = n; }
    void setPostSmooth (int n) { nu2 = n; }
    void setFinalSmooth (int n) { nuf = n; }
    void setBottomSmooth (int n) { nub = n; }

    void setBottomSolver (BottomSolver s) { bottom_solver = s; }
    void setBottomVerbose (int v) { bottom_verbose = v; }
    void setBottomMaxIter (int n) { bottom_maxiter = n; }
    void setCGVerbose (int v) { bottom_verbose = v; }
    void setCGMaxIter (int n) { bottom_maxiter = n; }

    void setAlwaysUseBNorm (int flag) { always_use_bnorm = flag; }

    void setFinalFillBC (int flag) { final_fill_bc = flag; }

    int numAMRLevels () const { return namrlevs; }

    void setNSolve (int flag) { do_nsolve = flag; }
    void setNSolveGridSize (int s) { nsolve_grid_size = s; }

private:

    int verbose = 1;
    int max_iters = 200;
    int do_fixed_number_of_iters = 0;

    int nu1 = 2;   // pre
    int nu2 = 2;   // post
    int nuf = 8;   // when smoother is used as bottom solver
    int nub = 0;   // aditional smoothing after bottom cg solver

    int max_fmg_iters = 0;

    BottomSolver bottom_solver = BottomSolver::bicgstab;
    int  bottom_verbose        = 0;
    int  bottom_maxiter        = 200;

    int always_use_bnorm = 0;

    int final_fill_bc = 0;

    MLLinOp& linop;
    int namrlevs;
    int finest_amr_lev;

    bool linop_prepared = false;

    // N Solve
    int do_nsolve = false;
    int nsolve_grid_size = 16;
    std::unique_ptr<MLLinOp> ns_linop;
    std::unique_ptr<MLMG> ns_mlmg;
    std::unique_ptr<MultiFab> ns_sol;
    std::unique_ptr<MultiFab> ns_rhs;

    // Hypre
#ifdef AMREX_USE_HYPRE
    std::unique_ptr<HypreABecLap2> hypre_solver;
    std::unique_ptr<MLMGBndry> hypre_bndry;
#endif

    // To avoid confusion, terms like sol, cor, rhs, res, ... etc. are
    // in the frame of the original equation, not the correction form
    Vector<std::unique_ptr<MultiFab> > sol_raii;
    Vector<MultiFab*>         sol;  // alias to argument a_sol
    Vector<MultiFab>          rhs;  // Copy of original rhs
                                    // L(sol) = rhs

    // First Vector: Amr levels.  0 is the coarest level
    // Second Vector: MG levels.  0 is the finest level
    Vector<Vector<MultiFab> >                   res; // = rhs - L(sol)
    Vector<Vector<std::unique_ptr<MultiFab> > > cor; // L(cor) = res
    Vector<Vector<std::unique_ptr<MultiFab> > > cor_hold;
    Vector<Vector<MultiFab> >                rescor; // = res - L(cor)  Residual of the correction form

    Vector<std::unique_ptr<iMultiFab> > fine_mask;

    enum timer_types { solve_time=0, iter_time, bottom_time, ntimers };
    Vector<Real> timer;

    void prepareForSolve (const Vector<MultiFab*>& a_sol, const Vector<MultiFab const*>& a_rhs);

    void prepareForNSolve ();

    void oneIter (int iter);

    void miniCycle (int alev);

    void mgVcycle (int amrlev, int mglev);
    void mgFcycle ();

    void bottomSolve ();
    void NSolve (MLMG& a_solver, MultiFab& a_sol, MultiFab& a_rhs);
    void actualBottomSolve ();

    void computeMLResidual (int amrlevmax);
    void computeResidual (int alev);
    void computeResWithCrseSolFineCor (int crse_amr_lev, int fine_amr_lev);
    void computeResWithCrseCorFineCor (int fine_amr_lev);
    void interpCorrection (int alev);
    void interpCorrection (int alev, int mglev);
    void addInterpCorrection (int alev, int mglev);

    void computeResOfCorrection (int amrlev, int mglev);

    Real ResNormInf (int amrlev, bool local = false);
    Real MLResNormInf (int alevmax, bool local = false);
    Real MLRhsNormInf (bool local = false);
    void buildFineMask ();

    void averageDownAndSync ();

    Real getNodalSum (int amrlve, int mglev, MultiFab& mf) const;

    void bottomSolveWithHypre (MultiFab& x, const MultiFab& b);
};

}

#endif
