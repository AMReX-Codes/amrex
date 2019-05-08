#ifndef AMREX_ML_MG_H_
#define AMREX_ML_MG_H_

#include <AMReX_MLLinOp.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MLCGSolver.H>

#ifdef AMREX_USE_HYPRE
#include <AMReX_Hypre.H>
#include <AMReX_HypreNodeLap.H>
#endif

namespace amrex {

#ifdef AMREX_USE_PETSC
class PETScABecLap;
#endif

class MLMG
{
public:

    friend class MLCGSolver;

    using BCMode = MLLinOp::BCMode;
    using Location = MLLinOp::Location;

    using BottomSolver = amrex::BottomSolver;

    MLMG (MLLinOp& a_lp);
    ~MLMG ();

    // Optional argument checkpoint_file is for debugging only.
    Real solve (const Vector<MultiFab*>& a_sol, const Vector<MultiFab const*>& a_rhs,
                Real a_tol_rel, Real a_tol_abs, const char* checkpoint_file = nullptr);

    void getGradSolution (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_grad_sol,
                          Location a_loc = Location::FaceCenter);

    /**
    * \brief For ``(alpha * a - beta * (del dot b grad)) phi = rhs``, flux means ``-b grad phi``
    *
    * \param a_flux
    * \param a_loc
    */
    void getFluxes (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_flux,
                    Location a_loc = Location::FaceCenter);
    void getFluxes (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_flux,
                    const Vector<MultiFab*> & a_sol,
                    Location a_loc = Location::FaceCenter);
    void getFluxes (const Vector<MultiFab*> & a_flux, Location a_loc = Location::CellCenter);
    void getFluxes (const Vector<MultiFab*> & a_flux, const Vector<MultiFab*> & a_sol, Location a_loc = Location::CellCenter);
    void compResidual (const Vector<MultiFab*>& a_res, const Vector<MultiFab*>& a_sol,
                       const Vector<MultiFab const*>& a_rhs);

    /**
    * \brief ``out = L(in)``. Note that, if no actual solve is needed, one could
    * turn off multigrid coarsening by constructing a MLLinOp object
    * with an appropriate LPInfo object (e.g., with LPInfo().setMaxCoarseningLevel(0)).
    *
    * \param out
    * \param in
    */
    void apply (const Vector<MultiFab*>& out, const Vector<MultiFab*>& in);

    void setVerbose (int v) noexcept { verbose = v; }
    void setMaxIter (int n) noexcept { max_iters = n; }
    void setMaxFmgIter (int n) noexcept { max_fmg_iters = n; }
    void setFixedIter (int nit) noexcept { do_fixed_number_of_iters = nit; }

    void setPreSmooth (int n) noexcept { nu1 = n; }
    void setPostSmooth (int n) noexcept { nu2 = n; }
    void setFinalSmooth (int n) noexcept { nuf = n; }
    void setBottomSmooth (int n) noexcept { nub = n; }

    void setBottomSolver (BottomSolver s) noexcept { bottom_solver = s; }
    void setBottomVerbose (int v) noexcept { bottom_verbose = v; }
    void setBottomMaxIter (int n) noexcept { bottom_maxiter = n; }
    void setBottomTolerance (Real t) noexcept { bottom_reltol = t; }
    void setCGVerbose (int v) noexcept { bottom_verbose = v; }
    void setCGMaxIter (int n) noexcept { bottom_maxiter = n; }
    void setCGTolerance (Real t) noexcept { bottom_reltol = t; }

    void setAlwaysUseBNorm (int flag) noexcept { always_use_bnorm = flag; }

    void setFinalFillBC (int flag) noexcept { final_fill_bc = flag; }

    int numAMRLevels () const noexcept { return namrlevs; }

    void setNSolve (int flag) noexcept { do_nsolve = flag; }
    void setNSolveGridSize (int s) noexcept { nsolve_grid_size = s; }

#ifdef AMREX_USE_HYPRE
    void setHypreInterface (Hypre::Interface f) noexcept {
        // must use ij interface for EB
#ifndef AMREX_USE_EB
        hypre_interface = f;
#endif
    }
#endif

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

    void computeVolInv ();
    void makeSolvable ();
    void makeSolvable (int amrlev, int mglev, MultiFab& mf);
    Real getNodalSum (int amrlev, int mglev, MultiFab& mf) const;

    void bottomSolveWithHypre (MultiFab& x, const MultiFab& b);

    void bottomSolveWithPETSc (MultiFab& x, const MultiFab& b);

    int bottomSolveWithCG (MultiFab& x, const MultiFab& b, MLCGSolver::Type type);

private:

    int verbose = 1;
    int max_iters = 200;
    int do_fixed_number_of_iters = 0;

    int nu1 = 2;       //!< pre
    int nu2 = 2;       //!< post
    int nuf = 8;       //!< when smoother is used as bottom solver
    int nub = 0;       //!< aditional smoothing after bottom cg solver

    int max_fmg_iters = 0;

    BottomSolver bottom_solver = BottomSolver::Default;
    int  bottom_verbose        = 0;
    int  bottom_maxiter        = 200;
    Real bottom_reltol         = 1.e-4;

    int always_use_bnorm = 0;

    int final_fill_bc = 0;

    MLLinOp& linop;
    int namrlevs;
    int finest_amr_lev;

    bool linop_prepared = false;
    long solve_called = 0;

    //! N Solve
    int do_nsolve = false;
    int nsolve_grid_size = 16;
    std::unique_ptr<MLLinOp> ns_linop;
    std::unique_ptr<MLMG> ns_mlmg;
    std::unique_ptr<MultiFab> ns_sol;
    std::unique_ptr<MultiFab> ns_rhs;

    //! Hypre
#ifdef AMREX_USE_HYPRE
#ifdef AMREX_USE_EB
    Hypre::Interface hypre_interface = Hypre::Interface::ij;
#else

    // Hypre::Interface hypre_interface = Hypre::Interface::structed;
    Hypre::Interface hypre_interface = Hypre::Interface::semi_structed;
//    Hypre::Interface hypre_interface = Hypre::Interface::ij;
#endif
    std::unique_ptr<Hypre> hypre_solver;
    std::unique_ptr<MLMGBndry> hypre_bndry;
    std::unique_ptr<HypreNodeLap> hypre_node_solver;
#endif

    //! PETSc
#ifdef AMREX_USE_PETSC
    std::unique_ptr<PETScABecLap> petsc_solver;
    std::unique_ptr<MLMGBndry> petsc_bndry;
#endif

    /**
    * \brief To avoid confusion, terms like sol, cor, rhs, res, ... etc. are
    * in the frame of the original equation, not the correction form
    */
    Vector<std::unique_ptr<MultiFab> > sol_raii;
    Vector<MultiFab*>         sol;      //!< alias to argument a_sol
    Vector<MultiFab>          rhs;      //!< Copy of original rhs
                                        //! L(sol) = rhs

    /**
    * \brief First Vector: Amr levels.  0 is the coarest level
    * Second Vector: MG levels.  0 is the finest level
    */
    Vector<Vector<MultiFab> >                   res;     //! = rhs - L(sol)
    Vector<Vector<std::unique_ptr<MultiFab> > > cor;     //!< L(cor) = res
    Vector<Vector<std::unique_ptr<MultiFab> > > cor_hold;
    Vector<Vector<MultiFab> >                   rescor;  //!< = res - L(cor)
                                                         //!  Residual of the correction form

    Vector<std::unique_ptr<iMultiFab> > fine_mask;

    Vector<Vector<Real> > volinv;      //!< used by makeSolvable

    Vector<std::unique_ptr<MultiFab> > scratch;

    enum timer_types { solve_time=0, iter_time, bottom_time, ntimers };
    Vector<Real> timer;

    void checkPoint (const Vector<MultiFab*>& a_sol, const Vector<MultiFab const*>& a_rhs,
                     Real a_tol_rel, Real a_tol_abs, const char* a_file_name) const;
};

}

#endif
