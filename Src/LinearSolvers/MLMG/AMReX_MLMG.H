#ifndef AMREX_ML_MG_H_
#define AMREX_ML_MG_H_
#include <AMReX_Config.H>

#include <AMReX_MLLinOp.H>
#include <AMReX_MLCGSolver.H>

namespace amrex {

template <typename MF>
class MLMGT
{
public:

    class error
        : public std::runtime_error
    {
    public :
        using std::runtime_error::runtime_error;
    };

    template <typename T> friend class MLCGSolverT;
    template <typename M> friend class GMRESMLMGT;

    using MFType = MF;
    using FAB = typename MLLinOpT<MF>::FAB;
    using RT  = typename MLLinOpT<MF>::RT;

    using BCMode   = typename MLLinOpT<MF>::BCMode;
    using Location = typename MLLinOpT<MF>::Location;

    using BottomSolver = amrex::BottomSolver;
    enum class CFStrategy : int {none,ghostnodes};

    MLMGT (MLLinOpT<MF>& a_lp);
    ~MLMGT ();

    MLMGT (MLMGT<MF> const&) = delete;
    MLMGT (MLMGT<MF> &&) = delete;
    MLMGT<MF>& operator= (MLMGT<MF> const&) = delete;
    MLMGT<MF>& operator= (MLMGT<MF> &&) = delete;

    // Optional argument checkpoint_file is for debugging only.
    template <typename AMF>
    RT solve (const Vector<AMF*>& a_sol, const Vector<AMF const*>& a_rhs,
              RT a_tol_rel, RT a_tol_abs, const char* checkpoint_file = nullptr);

    template <typename AMF>
    RT solve (std::initializer_list<AMF*> a_sol,
              std::initializer_list<AMF const*> a_rhs,
              RT a_tol_rel, RT a_tol_abs, const char* checkpoint_file = nullptr);

    template <typename AMF>
    void getGradSolution (const Vector<Array<AMF*,AMREX_SPACEDIM> >& a_grad_sol,
                          Location a_loc = Location::FaceCenter);

    template <typename AMF>
    void getGradSolution (std::initializer_list<Array<AMF*,AMREX_SPACEDIM>> a_grad_sol,
                          Location a_loc = Location::FaceCenter);

    /**
    * \brief For ``(alpha * a - beta * (del dot b grad)) phi = rhs``, flux means ``-b grad phi``
    */
    template <typename AMF>
    void getFluxes (const Vector<Array<AMF*,AMREX_SPACEDIM> >& a_flux,
                    Location a_loc = Location::FaceCenter);

    template <typename AMF>
    void getFluxes (std::initializer_list<Array<AMF*,AMREX_SPACEDIM>> a_flux,
                    Location a_loc = Location::FaceCenter);

    template <typename AMF>
    void getFluxes (const Vector<Array<AMF*,AMREX_SPACEDIM> >& a_flux,
                    const Vector<AMF*> & a_sol,
                    Location a_loc = Location::FaceCenter);

    template <typename AMF>
    void getFluxes (std::initializer_list<Array<AMF*,AMREX_SPACEDIM>> a_flux,
                    std::initializer_list<AMF*>  a_sol,
                    Location a_loc = Location::FaceCenter);

    template <typename AMF>
    void getFluxes (const Vector<AMF*> & a_flux,
                    Location a_loc = Location::CellCenter);

    template <typename AMF>
    void getFluxes (std::initializer_list<AMF*> a_flux,
                    Location a_loc = Location::CellCenter);

    template <typename AMF>
    void getFluxes (const Vector<AMF*> & a_flux,
                    const Vector<AMF*> & a_sol,
                    Location a_loc = Location::CellCenter);

    template <typename AMF>
    void getFluxes (std::initializer_list<AMF*> a_flux,
                    std::initializer_list<AMF*> a_sol,
                    Location a_loc = Location::CellCenter);

    void compResidual (const Vector<MF*>& a_res, const Vector<MF*>& a_sol,
                       const Vector<MF const*>& a_rhs);

#ifdef AMREX_USE_EB
    // Flux into the EB wall
    void getEBFluxes (const Vector<MF*>& a_eb_flux);
    void getEBFluxes (const Vector<MF*>& a_eb_flux, const Vector<MF*> & a_sol);
#endif

    /**
    * \brief ``out = L(in)``. Note that, if no actual solve is needed, one could
    * turn off multigrid coarsening by constructing a MLLinOp object
    * with an appropriate LPInfo object (e.g., with LPInfo().setMaxCoarseningLevel(0)).
    */
    void apply (const Vector<MF*>& out, const Vector<MF*>& in);

    void setThrowException (bool t) noexcept { throw_exception = t; }
    void setVerbose (int v) noexcept { verbose = v; }
    void setMaxIter (int n) noexcept { max_iters = n; }
    void setMaxFmgIter (int n) noexcept { max_fmg_iters = n; }
    void setFixedIter (int nit) noexcept { do_fixed_number_of_iters = nit; }

    void setPreSmooth (int n) noexcept { nu1 = n; }
    void setPostSmooth (int n) noexcept { nu2 = n; }
    void setFinalSmooth (int n) noexcept { nuf = n; }
    void setBottomSmooth (int n) noexcept { nub = n; }

    void setBottomSolver (BottomSolver s) noexcept { bottom_solver = s; }
    void setCFStrategy (CFStrategy a_cf_strategy) noexcept {cf_strategy = a_cf_strategy;}
    void setBottomVerbose (int v) noexcept { bottom_verbose = v; }
    void setBottomMaxIter (int n) noexcept { bottom_maxiter = n; }
    void setBottomTolerance (RT t) noexcept { bottom_reltol = t; }
    void setBottomToleranceAbs (RT t) noexcept { bottom_abstol = t;}
    RT getBottomToleranceAbs () noexcept{ return bottom_abstol; }

    void setAlwaysUseBNorm (int flag) noexcept { always_use_bnorm = flag; }

    void setFinalFillBC (int flag) noexcept { final_fill_bc = flag; }

    [[nodiscard]] int numAMRLevels () const noexcept { return namrlevs; }

    void setNSolve (int flag) noexcept { do_nsolve = flag; }
    void setNSolveGridSize (int s) noexcept { nsolve_grid_size = s; }

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    void setHypreInterface (Hypre::Interface f) noexcept {
        // must use ij interface for EB
#ifndef AMREX_USE_EB
        hypre_interface = f;
#else
        amrex::ignore_unused(f);
#endif
    }

    //! Set the namespace in input file for parsing HYPRE specific options
    void setHypreOptionsNamespace(const std::string& prefix) noexcept
    {
        hypre_options_namespace = prefix;
    }

    void setHypreOldDefault (bool l) noexcept {hypre_old_default = l;}
    void setHypreRelaxType (int n) noexcept {hypre_relax_type = n;}
    void setHypreRelaxOrder (int n) noexcept {hypre_relax_order = n;}
    void setHypreNumSweeps (int n) noexcept {hypre_num_sweeps = n;}
    void setHypreStrongThreshold (Real t) noexcept {hypre_strong_threshold = t;}
#endif

    template <typename AMF>
    void prepareForSolve (Vector<AMF*> const& a_sol, Vector<AMF const*> const& a_rhs);

    void prepareForNSolve ();

    void prepareLinOp ();

    void prepareMGcycle ();

    void prepareForGMRES ();

    void oneIter (int iter);

    void miniCycle (int amrlev);

    void mgVcycle (int amrlev, int mglev);
    void mgFcycle ();

    void bottomSolve ();
    void NSolve (MLMGT<MF>& a_solver, MF& a_sol, MF& a_rhs);
    void actualBottomSolve ();

    void computeMLResidual (int amrlevmax);
    void computeResidual (int alev);
    void computeResWithCrseSolFineCor (int calev, int falev);
    void computeResWithCrseCorFineCor (int falev);
    void interpCorrection (int alev);
    void interpCorrection (int alev, int mglev);
    void addInterpCorrection (int alev, int mglev);

    void computeResOfCorrection (int amrlev, int mglev);

    RT ResNormInf (int alev, bool local = false);
    RT MLResNormInf (int alevmax, bool local = false);
    RT MLRhsNormInf (bool local = false);

    void makeSolvable ();
    void makeSolvable (int amrlev, int mglev, MF& mf);

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    template <class TMF=MF,std::enable_if_t<std::is_same_v<TMF,MultiFab>,int> = 0>
    void bottomSolveWithHypre (MF& x, const MF& b);
#endif

#if defined(AMREX_USE_PETSC) && (AMREX_SPACEDIM > 1)
    template <class TMF=MF,std::enable_if_t<std::is_same_v<TMF,MultiFab>,int> = 0>
    void bottomSolveWithPETSc (MF& x, const MF& b);
#endif

    int bottomSolveWithCG (MF& x, const MF& b, typename MLCGSolverT<MF>::Type type);

    [[nodiscard]] RT getInitRHS () const noexcept { return m_rhsnorm0; }
    // Initial composite residual
    [[nodiscard]] RT getInitResidual () const noexcept { return m_init_resnorm0; }
    // Final composite residual
    [[nodiscard]] RT getFinalResidual () const noexcept { return m_final_resnorm0; }
    // Residuals on the *finest* AMR level after each iteration
    [[nodiscard]] Vector<RT> const& getResidualHistory () const noexcept { return m_iter_fine_resnorm0; }
    [[nodiscard]] int getNumIters () const noexcept { return m_iter_fine_resnorm0.size(); }
    [[nodiscard]] Vector<int> const& getNumCGIters () const noexcept { return m_niters_cg; }

    MLLinOpT<MF>& getLinOp () { return linop; }

private:

    bool throw_exception = false;
    int verbose = 1;

    int max_iters = 200;
    int do_fixed_number_of_iters = 0;

    int nu1 = 2;       //!< pre
    int nu2 = 2;       //!< post
    int nuf = 8;       //!< when smoother is used as bottom solver
    int nub = 0;       //!< additional smoothing after bottom cg solver

    int max_fmg_iters = 0;

    BottomSolver bottom_solver = BottomSolver::Default;
    CFStrategy cf_strategy     = CFStrategy::none;
    int  bottom_verbose        = 0;
    int  bottom_maxiter        = 200;
    RT bottom_reltol = std::is_same<RT,double>() ? RT(1.e-4) : RT(1.e-3);
    RT bottom_abstol = RT(-1.0);

    int always_use_bnorm = 0;

    int final_fill_bc = 0;

    MLLinOpT<MF>& linop;
    int ncomp;
    int namrlevs;
    int finest_amr_lev;

    bool linop_prepared = false;
    Long solve_called = 0;

    //! N Solve
    int do_nsolve = false;
    int nsolve_grid_size = 16;
    std::unique_ptr<MLLinOpT<MF>> ns_linop;
    std::unique_ptr<MLMGT<MF>> ns_mlmg;
    std::unique_ptr<MF> ns_sol;
    std::unique_ptr<MF> ns_rhs;

    //! Hypre
#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    // Hypre::Interface hypre_interface = Hypre::Interface::structed;
    // Hypre::Interface hypre_interface = Hypre::Interface::semi_structed;
    Hypre::Interface hypre_interface = Hypre::Interface::ij;

    std::unique_ptr<Hypre> hypre_solver;
    std::unique_ptr<MLMGBndryT<MF>> hypre_bndry;
    std::unique_ptr<HypreNodeLap> hypre_node_solver;

    std::string hypre_options_namespace = "hypre";
    bool hypre_old_default = true; // Falgout coarsening with modified classical interpolation
    int hypre_relax_type = 6;  // G-S/Jacobi hybrid relaxation
    int hypre_relax_order = 1; // uses C/F relaxation
    int hypre_num_sweeps = 2;  // Sweeps on each level
    Real hypre_strong_threshold = 0.25; // Hypre default is 0.25
#endif

    //! PETSc
#if defined(AMREX_USE_PETSC) && (AMREX_SPACEDIM > 1)
    std::unique_ptr<PETScABecLap> petsc_solver;
    std::unique_ptr<MLMGBndryT<MF>> petsc_bndry;
#endif

    /**
    * \brief To avoid confusion, terms like sol, cor, rhs, res, ... etc. are
    * in the frame of the original equation, not the correction form
    */
    Vector<MF> sol;      //!< Might be alias to argument a_sol
    Vector<MF> rhs;      //!< Copy of original rhs
                         //! L(sol) = rhs

    Vector<int> sol_is_alias;

    /**
    * \brief First Vector: Amr levels.  0 is the coarest level
    * Second Vector: MG levels.  0 is the finest level
    */
    Vector<Vector<MF> > res;     //! = rhs - L(sol)
    Vector<Vector<MF> > cor;     //!< L(cor) = res
    Vector<Vector<MF> > cor_hold;
    Vector<Vector<MF> > rescor;  //!< = res - L(cor)
                                 //!  Residual of the correction form

    enum timer_types { solve_time=0, iter_time, bottom_time, ntimers };
    Vector<double> timer;

    RT m_rhsnorm0 = RT(-1.0);
    RT m_init_resnorm0 = RT(-1.0);
    RT m_final_resnorm0 = RT(-1.0);
    Vector<int> m_niters_cg;
    Vector<RT> m_iter_fine_resnorm0; // Residual for each iteration at the finest level

    void checkPoint (const Vector<MultiFab*>& a_sol,
                     const Vector<MultiFab const*>& a_rhs,
                     RT a_tol_rel, RT a_tol_abs, const char* a_file_name) const;

};

template <typename MF>
MLMGT<MF>::MLMGT (MLLinOpT<MF>& a_lp)
    : linop(a_lp), ncomp(a_lp.getNComp()), namrlevs(a_lp.NAMRLevels()),
      finest_amr_lev(a_lp.NAMRLevels()-1)
{}

template <typename MF> MLMGT<MF>::~MLMGT () = default;

template <typename MF>
template <typename AMF>
auto
MLMGT<MF>::solve (std::initializer_list<AMF*> a_sol,
                  std::initializer_list<AMF const*> a_rhs,
                  RT a_tol_rel, RT a_tol_abs, const char* checkpoint_file) -> RT
{
    return solve(Vector<AMF*>(std::move(a_sol)),
                 Vector<AMF const*>(std::move(a_rhs)),
                 a_tol_rel, a_tol_abs, checkpoint_file);
}

template <typename MF>
template <typename AMF>
auto
MLMGT<MF>::solve (const Vector<AMF*>& a_sol, const Vector<AMF const*>& a_rhs,
                  RT a_tol_rel, RT a_tol_abs, const char* checkpoint_file) -> RT
{
    BL_PROFILE("MLMG::solve()");

    if constexpr (std::is_same<AMF,MultiFab>()) {
        if (checkpoint_file != nullptr) {
            checkPoint(a_sol, a_rhs, a_tol_rel, a_tol_abs, checkpoint_file);
        }
    }

    if (bottom_solver == BottomSolver::Default) {
        bottom_solver = linop.getDefaultBottomSolver();
    }

#if (defined(AMREX_USE_HYPRE) || defined(AMREX_USE_PETSC)) && (AMREX_SPACEDIM > 1)
    if (bottom_solver == BottomSolver::hypre || bottom_solver == BottomSolver::petsc) {
        int mo = linop.getMaxOrder();
        if (a_sol[0]->hasEBFabFactory()) {
            linop.setMaxOrder(2);
        } else {
            linop.setMaxOrder(std::min(3,mo));  // maxorder = 4 not supported
        }
    }
#endif

    bool is_nsolve = linop.m_parent;

    auto solve_start_time = amrex::second();

    RT& composite_norminf = m_final_resnorm0;

    m_niters_cg.clear();
    m_iter_fine_resnorm0.clear();

    prepareForSolve(a_sol, a_rhs);

    computeMLResidual(finest_amr_lev);

    bool local = true;
    RT resnorm0 = MLResNormInf(finest_amr_lev, local);
    RT rhsnorm0 = MLRhsNormInf(local);
    if (!is_nsolve) {
        ParallelAllReduce::Max<RT>({resnorm0, rhsnorm0}, ParallelContext::CommunicatorSub());

        if (verbose >= 1)
        {
            amrex::Print() << "MLMG: Initial rhs               = " << rhsnorm0 << "\n"
                           << "MLMG: Initial residual (resid0) = " << resnorm0 << "\n";
        }
    }

    m_init_resnorm0 = resnorm0;
    m_rhsnorm0 = rhsnorm0;

    RT max_norm;
    std::string norm_name;
    if (always_use_bnorm || rhsnorm0 >= resnorm0) {
        norm_name = "bnorm";
        max_norm = rhsnorm0;
    } else {
        norm_name = "resid0";
        max_norm = resnorm0;
    }
    const RT res_target = std::max(a_tol_abs, std::max(a_tol_rel,RT(1.e-16))*max_norm);

    if (!is_nsolve && resnorm0 <= res_target) {
        composite_norminf = resnorm0;
        if (verbose >= 1) {
            amrex::Print() << "MLMG: No iterations needed\n";
        }
    } else {
        auto iter_start_time = amrex::second();
        bool converged = false;

        const int niters = do_fixed_number_of_iters ? do_fixed_number_of_iters : max_iters;
        for (int iter = 0; iter < niters; ++iter)
        {
            oneIter(iter);

            converged = false;

            // Test convergence on the fine amr level
            computeResidual(finest_amr_lev);

            if (is_nsolve) { continue; }

            RT fine_norminf = ResNormInf(finest_amr_lev);
            m_iter_fine_resnorm0.push_back(fine_norminf);
            composite_norminf = fine_norminf;
            if (verbose >= 2) {
                amrex::Print() << "MLMG: Iteration " << std::setw(3) << iter+1 << " Fine resid/"
                               << norm_name << " = " << fine_norminf/max_norm << "\n";
            }
            bool fine_converged = (fine_norminf <= res_target);

            if (namrlevs == 1 && fine_converged) {
                converged = true;
            } else if (fine_converged) {
                // finest level is converged, but we still need to test the coarse levels
                computeMLResidual(finest_amr_lev-1);
                RT crse_norminf = MLResNormInf(finest_amr_lev-1);
                if (verbose >= 2) {
                    amrex::Print() << "MLMG: Iteration " << std::setw(3) << iter+1
                                   << " Crse resid/" << norm_name << " = "
                                   << crse_norminf/max_norm << "\n";
                }
                converged = (crse_norminf <= res_target);
                composite_norminf = std::max(fine_norminf, crse_norminf);
            } else {
                converged = false;
            }

            if (converged) {
                if (verbose >= 1) {
                    amrex::Print() << "MLMG: Final Iter. " << iter+1
                                   << " resid, resid/" << norm_name << " = "
                                   << composite_norminf << ", "
                                   << composite_norminf/max_norm << "\n";
                }
                break;
            } else {
                if (composite_norminf > RT(1.e20)*max_norm)
                {
                    if (verbose > 0) {
                        amrex::Print() << "MLMG: Failing to converge after " << iter+1 << " iterations."
                                       << " resid, resid/" << norm_name << " = "
                                       << composite_norminf << ", "
                                       << composite_norminf/max_norm << "\n";
                    }

                    if ( throw_exception ) {
                        throw error("MLMG blew up.");
                    } else {
                        amrex::Abort("MLMG failing so lets stop here");
                    }
                }
            }
        }

        if (!converged && do_fixed_number_of_iters == 0) {
            if (verbose > 0) {
                amrex::Print() << "MLMG: Failed to converge after " << max_iters << " iterations."
                               << " resid, resid/" << norm_name << " = "
                               << composite_norminf << ", "
                               << composite_norminf/max_norm << "\n";
            }

            if ( throw_exception ) {
                throw error("MLMG failed to converge.");
            } else {
                amrex::Abort("MLMG failed.");
            }
        }
        timer[iter_time] = amrex::second() - iter_start_time;
    }

    linop.postSolve(sol);

    IntVect ng_back = final_fill_bc ? IntVect(1) : IntVect(0);
    if (linop.hasHiddenDimension()) {
        ng_back[linop.hiddenDirection()] = 0;
    }
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (!sol_is_alias[alev]) {
            LocalCopy(*a_sol[alev], sol[alev], 0, 0, ncomp, ng_back);
        }
    }

    timer[solve_time] = amrex::second() - solve_start_time;
    if (verbose >= 1) {
        ParallelReduce::Max<double>(timer.data(), timer.size(), 0,
                                    ParallelContext::CommunicatorSub());
        if (ParallelContext::MyProcSub() == 0)
        {
            amrex::AllPrint() << "MLMG: Timers: Solve = " << timer[solve_time]
                              << " Iter = " << timer[iter_time]
                              << " Bottom = " << timer[bottom_time] << "\n";
        }
    }

    ++solve_called;

    return composite_norminf;
}

template <typename MF>
template <typename AMF>
void
MLMGT<MF>::getGradSolution (const Vector<Array<AMF*,AMREX_SPACEDIM> >& a_grad_sol, Location a_loc)
{
    BL_PROFILE("MLMG::getGradSolution()");
    for (int alev = 0; alev <= finest_amr_lev; ++alev) {
        if constexpr (std::is_same<AMF,MF>()) {
            linop.compGrad(alev, a_grad_sol[alev], sol[alev], a_loc);
        } else {
            Array<MF,AMREX_SPACEDIM> grad_sol;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                auto const& amf = *(a_grad_sol[alev][idim]);
                grad_sol[idim].define(boxArray(amf), DistributionMap(amf), ncomp, 0);
            }
            linop.compGrad(alev, GetArrOfPtrs(grad_sol), sol[alev], a_loc);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                LocalCopy(*a_grad_sol[alev][idim], grad_sol[idim], 0, 0, ncomp, IntVect(0));
            }
        }
    }
}

template <typename MF>
template <typename AMF>
void
MLMGT<MF>::getGradSolution (std::initializer_list<Array<AMF*,AMREX_SPACEDIM>> a_grad_sol, Location a_loc)
{
    getGradSolution(Vector<Array<AMF*,AMREX_SPACEDIM>>(std::move(a_grad_sol)), a_loc);
}

template <typename MF>
template <typename AMF>
void
MLMGT<MF>::getFluxes (const Vector<Array<AMF*,AMREX_SPACEDIM> >& a_flux,
                      Location a_loc)
{
    if (!linop.isCellCentered()) {
        amrex::Abort("Calling wrong getFluxes for nodal solver");
    }

    AMREX_ASSERT(sol.size() == a_flux.size());

    if constexpr (std::is_same<AMF,MF>()) {
        getFluxes(a_flux, GetVecOfPtrs(sol), a_loc);
    } else {
        Vector<Array<MF,AMREX_SPACEDIM>> fluxes(namrlevs);
        for (int ilev = 0; ilev < namrlevs; ++ilev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                auto const& amf = *(a_flux[ilev][idim]);
                fluxes[ilev][idim].define(boxArray(amf), DistributionMap(amf), ncomp, 0);
            }
        }
        getFluxes(GetVecOfArrOfPtrs(fluxes), GetVecOfPtrs(sol), a_loc);
        for (int ilev = 0; ilev < namrlevs; ++ilev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                LocalCopy(*a_flux[ilev][idim], fluxes[ilev][idim], 0, 0, ncomp, IntVect(0));
            }
        }
    }
}

template <typename MF>
template <typename AMF>
void
MLMGT<MF>::getFluxes (std::initializer_list<Array<AMF*,AMREX_SPACEDIM>> a_flux,
                      Location a_loc)
{
    getFluxes(Vector<Array<AMF*,AMREX_SPACEDIM>>(std::move(a_flux)), a_loc);
}

template <typename MF>
template <typename AMF>
void
MLMGT<MF>::getFluxes (const Vector<Array<AMF*,AMREX_SPACEDIM> >& a_flux,
                      const Vector<AMF*>& a_sol, Location a_loc)
{
    BL_PROFILE("MLMG::getFluxes()");

    if (!linop.isCellCentered()) {
       amrex::Abort("Calling wrong getFluxes for nodal solver");
    }

    if constexpr (std::is_same<AMF,MF>()) {
        linop.getFluxes(a_flux, a_sol, a_loc);
    } else {
        Vector<Array<MF,AMREX_SPACEDIM>> fluxes(namrlevs);
        for (int ilev = 0; ilev < namrlevs; ++ilev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                auto const& amf = *(a_flux[ilev][idim]);
                fluxes[ilev][idim].define(boxArray(amf), DistributionMap(amf), ncomp, 0);
            }
            LocalCopy(sol[ilev], *a_sol[ilev], 0, 0, ncomp, nGrowVect(sol[ilev]));
        }
        linop.getFluxes(GetVecOfArrOfPtrs(fluxes), GetVecOfPtrs(sol), a_loc);
        for (int ilev = 0; ilev < namrlevs; ++ilev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                LocalCopy(*a_flux[ilev][idim], fluxes[ilev][idim], 0, 0, ncomp, IntVect(0));
            }
        }
    }
}

template <typename MF>
template <typename AMF>
void
MLMGT<MF>::getFluxes (std::initializer_list<Array<AMF*,AMREX_SPACEDIM>> a_flux,
                      std::initializer_list<AMF*> a_sol, Location a_loc)
{
    getFluxes(Vector<Array<AMF*,AMREX_SPACEDIM>>(std::move(a_flux)),
              Vector<AMF*>(std::move(a_sol)), a_loc);
}

template <typename MF>
template <typename AMF>
void
MLMGT<MF>::getFluxes (const Vector<AMF*> & a_flux, Location a_loc)
{
    AMREX_ASSERT(sol.size() == a_flux.size());
    if constexpr (std::is_same<AMF,MF>()) {
        getFluxes(a_flux, GetVecOfPtrs(sol), a_loc);
    } else {
        Vector<MF> fluxes(namrlevs);
        for (int ilev = 0; ilev < namrlevs; ++ilev) {
            auto const& amf = *a_flux[ilev];
            fluxes[ilev].define(boxArray(amf), DistributionMap(amf), ncomp, 0);
        }
        getFluxes(GetVecOfPtrs(fluxes), GetVecOfPtrs(sol), a_loc);
        for (int ilev = 0; ilev < namrlevs; ++ilev) {
            LocalCopy(*a_flux[ilev], fluxes[ilev], 0, 0, ncomp, IntVect(0));
        }
    }
}

template <typename MF>
template <typename AMF>
void
MLMGT<MF>::getFluxes (std::initializer_list<AMF*> a_flux, Location a_loc)
{
    getFluxes(Vector<AMF*>(std::move(a_flux)), a_loc);
}

template <typename MF>
template <typename AMF>
void
MLMGT<MF>::getFluxes (const Vector<AMF*> & a_flux,
                      const Vector<AMF*>& a_sol, Location /*a_loc*/)
{
    AMREX_ASSERT(nComp(*a_flux[0]) >= AMREX_SPACEDIM);

    if constexpr (! std::is_same<AMF,MF>()) {
        for (int alev = 0; alev < namrlevs; ++alev) {
            LocalCopy(sol[alev], *a_sol[alev], 0, 0, ncomp, nGrowVect(sol[alev]));
        }
    }

    if (linop.isCellCentered())
    {
        Vector<Array<MF,AMREX_SPACEDIM> > ffluxes(namrlevs);
        for (int alev = 0; alev < namrlevs; ++alev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const int mglev = 0;
                int nghost = 0;
                if (cf_strategy == CFStrategy::ghostnodes) { nghost = linop.getNGrow(alev); }
                ffluxes[alev][idim].define(amrex::convert(linop.m_grids[alev][mglev],
                                                          IntVect::TheDimensionVector(idim)),
                                           linop.m_dmap[alev][mglev], ncomp, nghost, MFInfo(),
                                           *linop.m_factory[alev][mglev]);
            }
        }
        if constexpr (std::is_same<AMF,MF>()) {
            getFluxes(amrex::GetVecOfArrOfPtrs(ffluxes), a_sol, Location::FaceCenter);
        } else {
            getFluxes(amrex::GetVecOfArrOfPtrs(ffluxes), GetVecOfPtrs(sol), Location::FaceCenter);
        }
        for (int alev = 0; alev < namrlevs; ++alev) {
#ifdef AMREX_USE_EB
            EB_average_face_to_cellcenter(*a_flux[alev], 0, amrex::GetArrOfConstPtrs(ffluxes[alev]));
#else
            average_face_to_cellcenter(*a_flux[alev], 0, amrex::GetArrOfConstPtrs(ffluxes[alev]));
#endif
        }

    } else {
        if constexpr (std::is_same<AMF,MF>()) {
            linop.getFluxes(a_flux, a_sol);
        } else {
            Vector<MF> fluxes(namrlevs);
            for (int ilev = 0; ilev < namrlevs; ++ilev) {
                auto const& amf = *a_flux[ilev];
                fluxes[ilev].define(boxArray(amf), DistributionMap(amf), ncomp, 0);
            }
            linop.getFluxes(GetVecOfPtrs(fluxes), GetVecOfPtrs(sol));
            for (int ilev = 0; ilev < namrlevs; ++ilev) {
                LocalCopy(*a_flux[ilev], fluxes[ilev], 0, 0, ncomp, IntVect(0));
            }
        }
    }
}

template <typename MF>
template <typename AMF>
void
MLMGT<MF>::getFluxes (std::initializer_list<AMF*> a_flux,
                      std::initializer_list<AMF*> a_sol, Location a_loc)
{
    getFluxes(Vector<AMF*>(std::move(a_flux)),
              Vector<AMF*>(std::move(a_sol)), a_loc);
}

#ifdef AMREX_USE_EB
template <typename MF>
void
MLMGT<MF>::getEBFluxes (const Vector<MF*>& a_eb_flux)
{
    if (!linop.isCellCentered()) {
       amrex::Abort("getEBFluxes is for cell-centered only");
    }

    AMREX_ASSERT(sol.size() == a_eb_flux.size());
    getEBFluxes(a_eb_flux, GetVecOfPtrs(sol));
}

template <typename MF>
void
MLMGT<MF>::getEBFluxes (const Vector<MF*>& a_eb_flux, const Vector<MF*>& a_sol)
{
    BL_PROFILE("MLMG::getEBFluxes()");

    if (!linop.isCellCentered()) {
       amrex::Abort("getEBFluxes is for cell-centered only");
    }

    linop.getEBFluxes(a_eb_flux, a_sol);
}
#endif

template <typename MF>
void
MLMGT<MF>::compResidual (const Vector<MF*>& a_res, const Vector<MF*>& a_sol,
                         const Vector<MF const*>& a_rhs)
{
    BL_PROFILE("MLMG::compResidual()");

    IntVect ng_sol(1);
    if (linop.hasHiddenDimension()) { ng_sol[linop.hiddenDirection()] = 0; }

    sol.resize(namrlevs);
    sol_is_alias.resize(namrlevs,true);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (cf_strategy == CFStrategy::ghostnodes || nGrowVect(*a_sol[alev]) == ng_sol)
        {
            sol[alev] = linop.makeAlias(*a_sol[alev]);
            sol_is_alias[alev] = true;
        }
        else
        {
            if (sol_is_alias[alev])
            {
                sol[alev] = linop.make(alev, 0, ng_sol);
            }
            LocalCopy(sol[alev], *a_sol[alev], 0, 0, ncomp, IntVect(0));
        }
    }

    prepareLinOp();

    const auto& amrrr = linop.AMRRefRatio();

    for (int alev = finest_amr_lev; alev >= 0; --alev) {
        const MF* crse_bcdata = (alev > 0) ? &(sol[alev-1]) : nullptr;
        const MF* prhs = a_rhs[alev];
#if (AMREX_SPACEDIM != 3)
        int nghost = (cf_strategy == CFStrategy::ghostnodes) ? linop.getNGrow(alev) : 0;
        MF rhstmp(boxArray(*prhs), DistributionMap(*prhs), ncomp, nghost,
                  MFInfo(), *linop.Factory(alev));
        LocalCopy(rhstmp, *prhs, 0, 0, ncomp, IntVect(nghost));
        linop.applyMetricTerm(alev, 0, rhstmp);
        linop.unimposeNeumannBC(alev, rhstmp);
        linop.applyInhomogNeumannTerm(alev, rhstmp);
        prhs = &rhstmp;
#endif
        linop.solutionResidual(alev, *a_res[alev], sol[alev], *prhs, crse_bcdata);
        if (alev < finest_amr_lev) {
            linop.reflux(alev, *a_res[alev], sol[alev], *prhs,
                         *a_res[alev+1], sol[alev+1], *a_rhs[alev+1]);
            if (linop.isCellCentered()) {
#ifdef AMREX_USE_EB
                EB_average_down(*a_res[alev+1], *a_res[alev], 0, ncomp, amrrr[alev]);
#else
                average_down(*a_res[alev+1], *a_res[alev], 0, ncomp, amrrr[alev]);
#endif
            }
        }
    }


#if (AMREX_SPACEDIM != 3)
    for (int alev = 0; alev <= finest_amr_lev; ++alev) {
        linop.unapplyMetricTerm(alev, 0, *a_res[alev]);
    }
#endif
}

template <typename MF>
void
MLMGT<MF>::apply (const Vector<MF*>& out, const Vector<MF*>& a_in)
{
    BL_PROFILE("MLMG::apply()");

    Vector<MF*> in(namrlevs);
    Vector<MF> in_raii(namrlevs);
    Vector<MF> rh(namrlevs);
    int nghost = 0;
    IntVect ng_sol(1);
    if (linop.hasHiddenDimension()) { ng_sol[linop.hiddenDirection()] = 0; }

    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (cf_strategy == CFStrategy::ghostnodes)
        {
            nghost = linop.getNGrow(alev);
            in[alev] = a_in[alev];
        }
        else if (nGrowVect(*a_in[alev]) == ng_sol)
        {
            in[alev] = a_in[alev];
        }
        else
        {
            IntVect ng = ng_sol;
            if (cf_strategy == CFStrategy::ghostnodes) { ng = IntVect(nghost); }
            in_raii[alev] = linop.make(alev, 0, ng);
            LocalCopy(in_raii[alev], *a_in[alev], 0, 0, ncomp, IntVect(nghost));
            in[alev] = &(in_raii[alev]);
        }
        rh[alev] = linop.make(alev, 0, IntVect(nghost));
        setVal(rh[alev], RT(0.0));
    }

    prepareLinOp();

    for (int alev = 0; alev < namrlevs; ++alev) {
        linop.applyInhomogNeumannTerm(alev, rh[alev]);
    }

    const auto& amrrr = linop.AMRRefRatio();

    for (int alev = finest_amr_lev; alev >= 0; --alev) {
        const MF* crse_bcdata = (alev > 0) ? in[alev-1] : nullptr;
        linop.solutionResidual(alev, *out[alev], *in[alev], rh[alev], crse_bcdata);
        if (alev < finest_amr_lev) {
            linop.reflux(alev, *out[alev], *in[alev], rh[alev],
                         *out[alev+1], *in[alev+1], rh[alev+1]);
            if (linop.isCellCentered()) {
                if constexpr (IsMultiFabLike_v<MF>) {
#ifdef AMREX_USE_EB
                    EB_average_down(*out[alev+1], *out[alev], 0, nComp(*out[alev]), amrrr[alev]);
#else
                    average_down(*out[alev+1], *out[alev], 0, nComp(*out[alev]), amrrr[alev]);
#endif
                } else {
                    amrex::Abort("MLMG: TODO average_down for non-MultiFab");
                }
            }
        }
    }

#if (AMREX_SPACEDIM != 3)
    for (int alev = 0; alev <= finest_amr_lev; ++alev) {
        linop.unapplyMetricTerm(alev, 0, *out[alev]);
    }
#endif

    for (int alev = 0; alev <= finest_amr_lev; ++alev) {
        if (cf_strategy == CFStrategy::ghostnodes) { nghost = linop.getNGrow(alev); }
        Scale(*out[alev], RT(-1), 0, nComp(*out[alev]), nghost);
    }
}

template <typename MF>
template <typename AMF>
void
MLMGT<MF>::prepareForSolve (Vector<AMF*> const& a_sol, Vector<AMF const*> const& a_rhs)
{
    BL_PROFILE("MLMG::prepareForSolve()");

    AMREX_ASSERT(namrlevs <= a_sol.size());
    AMREX_ASSERT(namrlevs <= a_rhs.size());

    timer.assign(ntimers, 0.0);

    IntVect ng_rhs(0);
    IntVect ng_sol(1);
    if (linop.hasHiddenDimension()) { ng_sol[linop.hiddenDirection()] = 0; }

    if (!linop_prepared) {
        linop.prepareForSolve();
        linop_prepared = true;
    } else if (linop.needsUpdate()) {
        linop.update();

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
        hypre_solver.reset();
        hypre_bndry.reset();
        hypre_node_solver.reset();
#endif

#ifdef AMREX_USE_PETSC
        petsc_solver.reset();
        petsc_bndry.reset();
#endif
    }

    sol.resize(namrlevs);
    sol_is_alias.resize(namrlevs,false);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (cf_strategy == CFStrategy::ghostnodes)
        {
            if constexpr (std::is_same<AMF,MF>()) {
                sol[alev] = linop.makeAlias(*a_sol[alev]);
                sol_is_alias[alev] = true;
            } else {
                amrex::Abort("Type conversion not supported for CFStrategy::ghostnodes");
            }
        }
        else
        {
            if (nGrowVect(*a_sol[alev]) == ng_sol) {
                if constexpr (std::is_same<AMF,MF>()) {
                    sol[alev] = linop.makeAlias(*a_sol[alev]);
                    setBndry(sol[alev], RT(0.0), 0, ncomp);
                    sol_is_alias[alev] = true;
                }
            }
            if (!sol_is_alias[alev]) {
                if (!solve_called) {
                    sol[alev] = linop.make(alev, 0, ng_sol);
                }
                LocalCopy(sol[alev], *a_sol[alev], 0, 0, ncomp, IntVect(0));
                setBndry(sol[alev], RT(0.0), 0, ncomp);
            }
        }
    }

    rhs.resize(namrlevs);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (cf_strategy == CFStrategy::ghostnodes) { ng_rhs = IntVect(linop.getNGrow(alev)); }
        if (!solve_called) {
            rhs[alev] = linop.make(alev, 0, ng_rhs);
        }
        LocalCopy(rhs[alev], *a_rhs[alev], 0, 0, ncomp, ng_rhs);
        linop.applyMetricTerm(alev, 0, rhs[alev]);
        linop.unimposeNeumannBC(alev, rhs[alev]);
        linop.applyInhomogNeumannTerm(alev, rhs[alev]);
        linop.applyOverset(alev, rhs[alev]);
        linop.scaleRHS(alev, rhs[alev]);

#ifdef AMREX_USE_EB
        const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(linop.Factory(alev));
        if (factory && !factory->isAllRegular()) {
            if constexpr (std::is_same<MF,MultiFab>()) {
                EB_set_covered(rhs[alev], 0, ncomp, 0, RT(0.0));
                EB_set_covered(sol[alev], 0, ncomp, 0, RT(0.0));
            } else {
                amrex::Abort("TODO: MLMG with EB only works with MultiFab");
            }
        }
#endif
    }

    for (int falev = finest_amr_lev; falev > 0; --falev)
    {
        linop.averageDownSolutionRHS(falev-1, sol[falev-1], rhs[falev-1], sol[falev], rhs[falev]);
    }

    // enforce solvability if appropriate
    if (linop.isSingular(0) && linop.getEnforceSingularSolvable())
    {
        makeSolvable();
    }

    IntVect ng = linop.getNGrowVectRestriction();
    if (cf_strategy == CFStrategy::ghostnodes) { ng = ng_rhs; }
    if (!solve_called) {
        linop.make(res, ng);
        linop.make(rescor, ng);
    }
    for (int alev = 0; alev <= finest_amr_lev; ++alev)
    {
        const int nmglevs = linop.NMGLevels(alev);
        for (int mglev = 0; mglev < nmglevs; ++mglev)
        {
            setVal(res   [alev][mglev], RT(0.0));
            setVal(rescor[alev][mglev], RT(0.0));
        }
    }

    if (cf_strategy != CFStrategy::ghostnodes) { ng = ng_sol; }
    cor.resize(namrlevs);
    for (int alev = 0; alev <= finest_amr_lev; ++alev)
    {
        const int nmglevs = linop.NMGLevels(alev);
        cor[alev].resize(nmglevs);
        for (int mglev = 0; mglev < nmglevs; ++mglev)
        {
            if (!solve_called) {
                IntVect _ng = ng;
                if (cf_strategy == CFStrategy::ghostnodes) { _ng=IntVect(linop.getNGrow(alev,mglev)); }
                cor[alev][mglev] = linop.make(alev, mglev, _ng);
            }
            setVal(cor[alev][mglev], RT(0.0));
        }
    }

    cor_hold.resize(std::max(namrlevs-1,1));
    {
        const int alev = 0;
        const int nmglevs = linop.NMGLevels(alev);
        cor_hold[alev].resize(nmglevs);
        for (int mglev = 0; mglev < nmglevs-1; ++mglev)
        {
            if (!solve_called) {
                IntVect _ng = ng;
                if (cf_strategy == CFStrategy::ghostnodes) { _ng=IntVect(linop.getNGrow(alev,mglev)); }
                cor_hold[alev][mglev] = linop.make(alev, mglev, _ng);
            }
            setVal(cor_hold[alev][mglev], RT(0.0));
        }
    }
    for (int alev = 1; alev < finest_amr_lev; ++alev)
    {
        cor_hold[alev].resize(1);
        if (!solve_called) {
            IntVect _ng = ng;
            if (cf_strategy == CFStrategy::ghostnodes) { _ng=IntVect(linop.getNGrow(alev)); }
            cor_hold[alev][0] = linop.make(alev, 0, _ng);
        }
        setVal(cor_hold[alev][0], RT(0.0));
    }

    if (linop.m_parent // no embedded N-Solve
        || !linop.supportNSolve())
    {
        do_nsolve = false;
    }

    if (do_nsolve && ns_linop == nullptr)
    {
        prepareForNSolve();
    }

    if (verbose >= 2) {
        amrex::Print() << "MLMG: # of AMR levels: " << namrlevs << "\n"
                       << "      # of MG levels on the coarsest AMR level: " << linop.NMGLevels(0)
                       << "\n";
        if (ns_linop) {
            amrex::Print() << "      # of MG levels in N-Solve: " << ns_linop->NMGLevels(0) << "\n"
                           << "      # of grids in N-Solve: " << ns_linop->m_grids[0][0].size() << "\n";
        }
    }
}

template <typename MF>
void
MLMGT<MF>::prepareLinOp ()
{
    if (!linop_prepared) {
        linop.prepareForSolve();
        linop_prepared = true;
    } else if (linop.needsUpdate()) {
        linop.update();
    }
}

template <typename MF>
void
MLMGT<MF>::prepareForGMRES ()
{
    prepareLinOp();
    linop.prepareForGMRES();
}

template <typename MF>
void
MLMGT<MF>::prepareMGcycle ()
{
    if (res.empty()) {
        IntVect ng = linop.getNGrowVectRestriction();
        linop.make(res, ng);
        linop.make(rescor, ng);

        for (int alev = 0; alev <= finest_amr_lev; ++alev)
        {
            const int nmglevs = linop.NMGLevels(alev);
            for (int mglev = 0; mglev < nmglevs; ++mglev)
            {
                setVal(res   [alev][mglev], RT(0.0));
                setVal(rescor[alev][mglev], RT(0.0));
            }
        }

        IntVect ng_sol(1);
        if (linop.hasHiddenDimension()) { ng_sol[linop.hiddenDirection()] = 0; }
        ng = ng_sol;

        cor.resize(namrlevs);
        for (int alev = 0; alev <= finest_amr_lev; ++alev)
        {
            const int nmglevs = linop.NMGLevels(alev);
            cor[alev].resize(nmglevs);
            for (int mglev = 0; mglev < nmglevs; ++mglev)
            {
                cor[alev][mglev] = linop.make(alev, mglev, ng);
                setVal(cor[alev][mglev], RT(0.0));
            }
        }

        cor_hold.resize(std::max(namrlevs-1,1));
        {
            const int alev = 0;
            const int nmglevs = linop.NMGLevels(alev);
            cor_hold[alev].resize(nmglevs);
            for (int mglev = 0; mglev < nmglevs-1; ++mglev)
            {
                cor_hold[alev][mglev] = linop.make(alev, mglev, ng);
                setVal(cor_hold[alev][mglev], RT(0.0));
            }
        }
        for (int alev = 1; alev < finest_amr_lev; ++alev)
        {
            cor_hold[alev].resize(1);
            cor_hold[alev][0] = linop.make(alev, 0, ng);
            setVal(cor_hold[alev][0], RT(0.0));
        }
    }
}

template <typename MF>
void
MLMGT<MF>::prepareForNSolve ()
{
    if constexpr (IsMultiFabLike_v<MF>) {
        ns_linop = linop.makeNLinOp(nsolve_grid_size);

        int nghost = 0;
        if (cf_strategy == CFStrategy::ghostnodes) { nghost = linop.getNGrow(); }

        const BoxArray& ba = (*ns_linop).m_grids[0][0];
        const DistributionMapping& dm =(*ns_linop).m_dmap[0][0];

        int ng = 1;
        if (cf_strategy == CFStrategy::ghostnodes) { ng = nghost; }
        ns_sol = std::make_unique<MF>(ba, dm, ncomp, ng, MFInfo(), *(ns_linop->Factory(0,0)));
        ng = 0;
        if (cf_strategy == CFStrategy::ghostnodes) { ng = nghost; }
        ns_rhs = std::make_unique<MF>(ba, dm, ncomp, ng, MFInfo(), *(ns_linop->Factory(0,0)));
        setVal(*ns_sol, RT(0.0));
        setVal(*ns_rhs, RT(0.0));

        ns_linop->setLevelBC(0, ns_sol.get());

        ns_mlmg = std::make_unique<MLMGT<MF>>(*ns_linop);
        ns_mlmg->setVerbose(0);
        ns_mlmg->setFixedIter(1);
        ns_mlmg->setMaxFmgIter(20);
        ns_mlmg->setBottomSolver(BottomSolver::smoother);
    }
}

// in  : Residual (res) on the finest AMR level
// out : sol on all AMR levels
template <typename MF>
void MLMGT<MF>::oneIter (int iter)
{
    BL_PROFILE("MLMG::oneIter()");

    for (int alev = finest_amr_lev; alev > 0; --alev)
    {
        miniCycle(alev);

        IntVect nghost(0);
        if (cf_strategy == CFStrategy::ghostnodes) { nghost = IntVect(linop.getNGrow(alev)); }
        LocalAdd(sol[alev], cor[alev][0], 0, 0, ncomp, nghost);

        // compute residual for the coarse AMR level
        computeResWithCrseSolFineCor(alev-1,alev);

        if (alev != finest_amr_lev) {
            std::swap(cor_hold[alev][0], cor[alev][0]); // save it for the up cycle
        }
    }

    // coarsest amr level
    {
        // enforce solvability if appropriate
        if (linop.isSingular(0) && linop.getEnforceSingularSolvable())
        {
            makeSolvable(0,0,res[0][0]);
        }

        if (iter < max_fmg_iters) {
            mgFcycle();
        } else {
            mgVcycle(0, 0);
        }

        IntVect nghost(0);
        if (cf_strategy == CFStrategy::ghostnodes) { nghost = IntVect(linop.getNGrow(0)); }
        LocalAdd(sol[0], cor[0][0], 0, 0, ncomp, nghost);
    }

    for (int alev = 1; alev <= finest_amr_lev; ++alev)
    {
        // (Fine AMR correction) = I(Coarse AMR correction)
        interpCorrection(alev);

        IntVect nghost(0);
        if (cf_strategy == CFStrategy::ghostnodes) { nghost = IntVect(linop.getNGrow(alev)); }
        LocalAdd(sol[alev], cor[alev][0], 0, 0, ncomp, nghost);

        if (alev != finest_amr_lev) {
            LocalAdd(cor_hold[alev][0], cor[alev][0], 0, 0, ncomp, nghost);
        }

        // Update fine AMR level correction
        computeResWithCrseCorFineCor(alev);

        miniCycle(alev);

        LocalAdd(sol[alev], cor[alev][0], 0, 0, ncomp, nghost);

        if (alev != finest_amr_lev) {
            LocalAdd(cor[alev][0], cor_hold[alev][0], 0, 0, ncomp, nghost);
        }
    }

    linop.averageDownAndSync(sol);
}

template <typename MF>
void
MLMGT<MF>::miniCycle (int amrlev)
{
    BL_PROFILE("MLMG::miniCycle()");
    const int mglev = 0;
    mgVcycle(amrlev, mglev);
}

// in   : Residual (res)
// out  : Correction (cor) from bottom to this function's local top
template <typename MF>
void
MLMGT<MF>::mgVcycle (int amrlev, int mglev_top)
{
    BL_PROFILE("MLMG::mgVcycle()");

    const int mglev_bottom = linop.NMGLevels(amrlev) - 1;

    for (int mglev = mglev_top; mglev < mglev_bottom; ++mglev)
    {
        BL_PROFILE_VAR("MLMG::mgVcycle_down::"+std::to_string(mglev), blp_mgv_down_lev);

        if (verbose >= 4)
        {
            RT norm = norminf(res[amrlev][mglev],0,ncomp,IntVect(0));
            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev
                           << "   DN: Norm before smooth " << norm << "\n";
        }

        setVal(cor[amrlev][mglev], RT(0.0));
        bool skip_fillboundary = true;
        for (int i = 0; i < nu1; ++i) {
            linop.smooth(amrlev, mglev, cor[amrlev][mglev], res[amrlev][mglev], skip_fillboundary);
            skip_fillboundary = false;
        }

        // rescor = res - L(cor)
        computeResOfCorrection(amrlev, mglev);

        if (verbose >= 4)
        {
            RT norm = norminf(rescor[amrlev][mglev],0,ncomp,IntVect(0));
            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev
                           << "   DN: Norm after  smooth " << norm << "\n";
        }

        // res_crse = R(rescor_fine); this provides res/b to the level below
        linop.restriction(amrlev, mglev+1, res[amrlev][mglev+1], rescor[amrlev][mglev]);
    }

    BL_PROFILE_VAR("MLMG::mgVcycle_bottom", blp_bottom);
    if (amrlev == 0)
    {
        if (verbose >= 4)
        {
            RT norm = norminf(res[amrlev][mglev_bottom],0,ncomp,IntVect(0));
            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev_bottom
                           << "   DN: Norm before bottom " << norm << "\n";
        }
        bottomSolve();
        if (verbose >= 4)
        {
            computeResOfCorrection(amrlev, mglev_bottom);
            RT norm = norminf(rescor[amrlev][mglev_bottom],0,ncomp,IntVect(0));
            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev_bottom
                           << "   UP: Norm after  bottom " << norm << "\n";
        }
    }
    else
    {
        if (verbose >= 4)
        {
            RT norm = norminf(res[amrlev][mglev_bottom],0,ncomp,IntVect(0));
            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev_bottom
                           << "       Norm before smooth " << norm << "\n";
        }
        setVal(cor[amrlev][mglev_bottom], RT(0.0));
        bool skip_fillboundary = true;
        for (int i = 0; i < nu1; ++i) {
            linop.smooth(amrlev, mglev_bottom, cor[amrlev][mglev_bottom],
                         res[amrlev][mglev_bottom], skip_fillboundary);
            skip_fillboundary = false;
        }
        if (verbose >= 4)
        {
            computeResOfCorrection(amrlev, mglev_bottom);
            RT norm = norminf(rescor[amrlev][mglev_bottom],0,ncomp,IntVect(0));
            amrex::Print() << "AT LEVEL "  << amrlev  << " " << mglev_bottom
                           << "       Norm after  smooth " << norm << "\n";
        }
    }
    BL_PROFILE_VAR_STOP(blp_bottom);

    for (int mglev = mglev_bottom-1; mglev >= mglev_top; --mglev)
    {
        BL_PROFILE_VAR("MLMG::mgVcycle_up::"+std::to_string(mglev), blp_mgv_up_lev);
        // cor_fine += I(cor_crse)
        addInterpCorrection(amrlev, mglev);
        if (verbose >= 4)
        {
            computeResOfCorrection(amrlev, mglev);
            RT norm = norminf(rescor[amrlev][mglev],0,ncomp,IntVect(0));
            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev
                           << "   UP: Norm before smooth " << norm << "\n";
        }
        for (int i = 0; i < nu2; ++i) {
            linop.smooth(amrlev, mglev, cor[amrlev][mglev], res[amrlev][mglev]);
        }

        if (cf_strategy == CFStrategy::ghostnodes) { computeResOfCorrection(amrlev, mglev); }

        if (verbose >= 4)
        {
            computeResOfCorrection(amrlev, mglev);
            RT norm = norminf(rescor[amrlev][mglev],0,ncomp,IntVect(0));
            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev
                           << "   UP: Norm after  smooth " << norm << "\n";
        }
    }
}

// FMG cycle on the coarsest AMR level.
// in:  Residual on the top MG level (i.e., 0)
// out: Correction (cor) on all MG levels
template <typename MF>
void
MLMGT<MF>::mgFcycle ()
{
    BL_PROFILE("MLMG::mgFcycle()");

#ifdef AMREX_USE_EB
    AMREX_ASSERT(linop.isCellCentered());
#endif

    const int amrlev = 0;
    const int mg_bottom_lev = linop.NMGLevels(amrlev) - 1;
    IntVect nghost(0);
    if (cf_strategy == CFStrategy::ghostnodes) { nghost = IntVect(linop.getNGrow(amrlev)); }

    for (int mglev = 1; mglev <= mg_bottom_lev; ++mglev)
    {
        linop.avgDownResMG(mglev, res[amrlev][mglev], res[amrlev][mglev-1]);
    }

    bottomSolve();

    for (int mglev = mg_bottom_lev-1; mglev >= 0; --mglev)
    {
        // cor_fine = I(cor_crse)
        interpCorrection(amrlev, mglev);

        // rescor = res - L(cor)
        computeResOfCorrection(amrlev, mglev);
        // res = rescor; this provides b to the vcycle below
        LocalCopy(res[amrlev][mglev], rescor[amrlev][mglev], 0, 0, ncomp, nghost);

        // save cor; do v-cycle; add the saved to cor
        std::swap(cor[amrlev][mglev], cor_hold[amrlev][mglev]);
        mgVcycle(amrlev, mglev);
        LocalAdd(cor[amrlev][mglev], cor_hold[amrlev][mglev], 0, 0, ncomp, nghost);
    }
}

// At the true bottom of the coarsest AMR level.
// in  : Residual (res) as b
// out : Correction (cor) as x
template <typename MF>
void
MLMGT<MF>::bottomSolve ()
{
    if (do_nsolve)
    {
        NSolve(*ns_mlmg, *ns_sol, *ns_rhs);
    }
    else
    {
        actualBottomSolve();
    }
}

template <typename MF>
void
MLMGT<MF>::NSolve (MLMGT<MF>& a_solver, MF& a_sol, MF& a_rhs)
{
    BL_PROFILE("MLMG::NSolve()");

    setVal(a_sol, RT(0.0));

    MF const& res_bottom = res[0].back();
    if (BoxArray::SameRefs(boxArray(a_rhs),boxArray(res_bottom)) &&
        DistributionMapping::SameRefs(DistributionMap(a_rhs),DistributionMap(res_bottom)))
    {
        LocalCopy(a_rhs, res_bottom, 0, 0, ncomp, IntVect(0));
    } else {
        setVal(a_rhs, RT(0.0));
        ParallelCopy(a_rhs, res_bottom, 0, 0, ncomp);
    }

    a_solver.solve(Vector<MF*>{&a_sol}, Vector<MF const*>{&a_rhs},
                   RT(-1.0), RT(-1.0));

    linop.copyNSolveSolution(cor[0].back(), a_sol);
}

template <typename MF>
void
MLMGT<MF>::actualBottomSolve ()
{
    BL_PROFILE("MLMG::actualBottomSolve()");

    if (!linop.isBottomActive()) { return; }

    auto bottom_start_time = amrex::second();

    ParallelContext::push(linop.BottomCommunicator());

    const int amrlev = 0;
    const int mglev = linop.NMGLevels(amrlev) - 1;
    auto& x = cor[amrlev][mglev];
    auto& b = res[amrlev][mglev];

    setVal(x, RT(0.0));

    if (bottom_solver == BottomSolver::smoother)
    {
        bool skip_fillboundary = true;
        for (int i = 0; i < nuf; ++i) {
            linop.smooth(amrlev, mglev, x, b, skip_fillboundary);
            skip_fillboundary = false;
        }
    }
    else
    {
        MF* bottom_b = &b;
        MF raii_b;
        if (linop.isBottomSingular() && linop.getEnforceSingularSolvable())
        {
            const IntVect ng = nGrowVect(b);
            raii_b = linop.make(amrlev, mglev, ng);
            LocalCopy(raii_b, b, 0, 0, ncomp, ng);
            bottom_b = &raii_b;

            makeSolvable(amrlev,mglev,*bottom_b);
        }

        if (bottom_solver == BottomSolver::hypre)
        {
#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
            if constexpr (std::is_same<MF,MultiFab>()) {
                bottomSolveWithHypre(x, *bottom_b);
            } else
#endif
            {
                amrex::Abort("Using Hypre as bottom solver not supported in this case");
            }
        }
        else if (bottom_solver == BottomSolver::petsc)
        {
#if defined(AMREX_USE_PETSC) && (AMREX_SPACEDIM > 1)
            if constexpr (std::is_same<MF,MultiFab>()) {
                bottomSolveWithPETSc(x, *bottom_b);
            } else
#endif
            {
                amrex::Abort("Using PETSc as bottom solver not supported in this case");
            }
        }
        else
        {
            typename MLCGSolverT<MF>::Type cg_type;
            if (bottom_solver == BottomSolver::cg ||
                bottom_solver == BottomSolver::cgbicg) {
                cg_type = MLCGSolverT<MF>::Type::CG;
            } else {
                cg_type = MLCGSolverT<MF>::Type::BiCGStab;
            }
            int ret = bottomSolveWithCG(x, *bottom_b, cg_type);
            // If the MLMG solve failed then set the correction to zero
            if (ret != 0) {
                setVal(cor[amrlev][mglev], RT(0.0));
                if (bottom_solver == BottomSolver::cgbicg ||
                    bottom_solver == BottomSolver::bicgcg) {
                    if (bottom_solver == BottomSolver::cgbicg) {
                        cg_type = MLCGSolverT<MF>::Type::BiCGStab; // switch to bicg
                    } else {
                        cg_type = MLCGSolverT<MF>::Type::CG; // switch to cg
                    }
                    ret = bottomSolveWithCG(x, *bottom_b, cg_type);
                    if (ret != 0) {
                        setVal(cor[amrlev][mglev], RT(0.0));
                    } else { // switch permanently
                        if (cg_type == MLCGSolverT<MF>::Type::CG) {
                            bottom_solver = BottomSolver::cg;
                        } else {
                            bottom_solver = BottomSolver::bicgstab;
                        }
                    }
                }
            }
            const int n = (ret==0) ? nub : nuf;
            for (int i = 0; i < n; ++i) {
                linop.smooth(amrlev, mglev, x, b);
            }
        }
    }

    ParallelContext::pop();

    if (! timer.empty()) {
        timer[bottom_time] += amrex::second() - bottom_start_time;
    }
}

template <typename MF>
int
MLMGT<MF>::bottomSolveWithCG (MF& x, const MF& b, typename MLCGSolverT<MF>::Type type)
{
    MLCGSolverT<MF> cg_solver(linop);
    cg_solver.setSolver(type);
    cg_solver.setVerbose(bottom_verbose);
    cg_solver.setMaxIter(bottom_maxiter);
    cg_solver.setInitSolnZeroed(true);
    if (cf_strategy == CFStrategy::ghostnodes) { cg_solver.setNGhost(linop.getNGrow()); }

    int ret = cg_solver.solve(x, b, bottom_reltol, bottom_abstol);
    if (ret != 0 && verbose > 1) {
        amrex::Print() << "MLMG: Bottom solve failed.\n";
    }
    m_niters_cg.push_back(cg_solver.getNumIters());
    return ret;
}

// Compute multi-level Residual (res) up to amrlevmax.
template <typename MF>
void
MLMGT<MF>::computeMLResidual (int amrlevmax)
{
    BL_PROFILE("MLMG::computeMLResidual()");

    const int mglev = 0;
    for (int alev = amrlevmax; alev >= 0; --alev) {
        const MF* crse_bcdata = (alev > 0) ? &(sol[alev-1]) : nullptr;
        linop.solutionResidual(alev, res[alev][mglev], sol[alev], rhs[alev], crse_bcdata);
        if (alev < finest_amr_lev) {
            linop.reflux(alev, res[alev][mglev], sol[alev], rhs[alev],
                         res[alev+1][mglev], sol[alev+1], rhs[alev+1]);
        }
    }
}

// Compute single AMR level residual without masking.
template <typename MF>
void
MLMGT<MF>::computeResidual (int alev)
{
    BL_PROFILE("MLMG::computeResidual()");
    const MF* crse_bcdata = (alev > 0) ? &(sol[alev-1]) : nullptr;
    linop.solutionResidual(alev, res[alev][0], sol[alev], rhs[alev], crse_bcdata);
}

// Compute coarse AMR level composite residual with coarse solution and fine correction
template <typename MF>
void
MLMGT<MF>::computeResWithCrseSolFineCor (int calev, int falev)
{
    BL_PROFILE("MLMG::computeResWithCrseSolFineCor()");

    IntVect nghost(0);
    if (cf_strategy == CFStrategy::ghostnodes) {
        nghost = IntVect(std::min(linop.getNGrow(falev),linop.getNGrow(calev)));
    }

    MF&       crse_sol = sol[calev];
    const MF& crse_rhs = rhs[calev];
    MF&       crse_res = res[calev][0];

    MF&       fine_sol = sol[falev];
    const MF& fine_rhs = rhs[falev];
    MF&       fine_cor = cor[falev][0];
    MF&       fine_res = res[falev][0];
    MF&    fine_rescor = rescor[falev][0];

    const MF* crse_bcdata = (calev > 0) ? &(sol[calev-1]) : nullptr;
    linop.solutionResidual(calev, crse_res, crse_sol, crse_rhs, crse_bcdata);

    linop.correctionResidual(falev, 0, fine_rescor, fine_cor, fine_res, BCMode::Homogeneous);
    LocalCopy(fine_res, fine_rescor, 0, 0, ncomp, nghost);

    linop.reflux(calev, crse_res, crse_sol, crse_rhs, fine_res, fine_sol, fine_rhs);

    linop.avgDownResAmr(calev, crse_res, fine_res);
}

// Compute fine AMR level residual fine_res = fine_res - L(fine_cor) with coarse providing BC.
template <typename MF>
void
MLMGT<MF>::computeResWithCrseCorFineCor (int falev)
{
    BL_PROFILE("MLMG::computeResWithCrseCorFineCor()");

    IntVect nghost(0);
    if (cf_strategy == CFStrategy::ghostnodes) {
        nghost = IntVect(linop.getNGrow(falev));
    }

    const MF& crse_cor = cor[falev-1][0];

    MF& fine_cor    = cor   [falev][0];
    MF& fine_res    = res   [falev][0];
    MF& fine_rescor = rescor[falev][0];

    // fine_rescor = fine_res - L(fine_cor)
    linop.correctionResidual(falev, 0, fine_rescor, fine_cor, fine_res,
                             BCMode::Inhomogeneous, &crse_cor);
    LocalCopy(fine_res, fine_rescor, 0, 0, ncomp, nghost);
}

// Interpolate correction from coarse to fine AMR level.
template <typename MF>
void
MLMGT<MF>::interpCorrection (int alev)
{
    BL_PROFILE("MLMG::interpCorrection_1");

    IntVect nghost(0);
    if (cf_strategy == CFStrategy::ghostnodes) {
        nghost = IntVect(linop.getNGrow(alev));
    }

    MF const& crse_cor = cor[alev-1][0];
    MF      & fine_cor = cor[alev  ][0];

    const Geometry& crse_geom = linop.Geom(alev-1,0);

    int ng_src = 0;
    int ng_dst = linop.isCellCentered() ? 1 : 0;
    if (cf_strategy == CFStrategy::ghostnodes)
    {
        ng_src = linop.getNGrow(alev-1);
        ng_dst = linop.getNGrow(alev-1);
    }

    MF cfine = linop.makeCoarseAmr(alev, IntVect(ng_dst));
    setVal(cfine, RT(0.0));
    ParallelCopy(cfine, crse_cor, 0, 0, ncomp, IntVect(ng_src), IntVect(ng_dst),
                 crse_geom.periodicity());

    linop.interpolationAmr(alev, fine_cor, cfine, nghost); // NOLINT(readability-suspicious-call-argument)
}

// Interpolate correction between MG levels
// inout: Correction (cor) on coarse MG lev.  (out due to FillBoundary)
// out  : Correction (cor) on fine MG lev.
template <typename MF>
void
MLMGT<MF>::interpCorrection (int alev, int mglev)
{
    BL_PROFILE("MLMG::interpCorrection_2");

    MF& crse_cor = cor[alev][mglev+1];
    MF& fine_cor = cor[alev][mglev  ];
    linop.interpAssign(alev, mglev, fine_cor, crse_cor);
}

// (Fine MG level correction) += I(Coarse MG level correction)
template <typename MF>
void
MLMGT<MF>::addInterpCorrection (int alev, int mglev)
{
    BL_PROFILE("MLMG::addInterpCorrection()");

    const MF& crse_cor = cor[alev][mglev+1];
    MF&       fine_cor = cor[alev][mglev  ];

    MF cfine;
    const MF* cmf;

    if (linop.isMFIterSafe(alev, mglev, mglev+1))
    {
        cmf = &crse_cor;
    }
    else
    {
        cfine = linop.makeCoarseMG(alev, mglev, IntVect(0));
        ParallelCopy(cfine, crse_cor, 0, 0, ncomp);
        cmf = &cfine;
    }

    linop.interpolation(alev, mglev, fine_cor, *cmf);
}

// Compute rescor = res - L(cor)
// in   : res
// inout: cor (out due to FillBoundary in linop.correctionResidual)
// out  : rescor
template <typename MF>
void
MLMGT<MF>::computeResOfCorrection (int amrlev, int mglev)
{
    BL_PROFILE("MLMG:computeResOfCorrection()");
    MF      & x =    cor[amrlev][mglev];
    const MF& b =    res[amrlev][mglev];
    MF      & r = rescor[amrlev][mglev];
    linop.correctionResidual(amrlev, mglev, r, x, b, BCMode::Homogeneous);
}

// Compute single-level masked inf-norm of Residual (res).
template <typename MF>
auto
MLMGT<MF>::ResNormInf (int alev, bool local) -> RT
{
    BL_PROFILE("MLMG::ResNormInf()");
    return linop.normInf(alev, res[alev][0], local);
}

// Computes multi-level masked inf-norm of Residual (res).
template <typename MF>
auto
MLMGT<MF>::MLResNormInf (int alevmax, bool local) -> RT
{
    BL_PROFILE("MLMG::MLResNormInf()");
    RT r = RT(0.0);
    for (int alev = 0; alev <= alevmax; ++alev)
    {
        r = std::max(r, ResNormInf(alev,true));
    }
    if (!local) { ParallelAllReduce::Max(r, ParallelContext::CommunicatorSub()); }
    return r;
}

// Compute multi-level masked inf-norm of RHS (rhs).
template <typename MF>
auto
MLMGT<MF>::MLRhsNormInf (bool local) -> RT
{
    BL_PROFILE("MLMG::MLRhsNormInf()");
    RT r = RT(0.0);
    for (int alev = 0; alev <= finest_amr_lev; ++alev) {
        auto t = linop.normInf(alev, rhs[alev], true);
        r = std::max(r, t);
    }
    if (!local) { ParallelAllReduce::Max(r, ParallelContext::CommunicatorSub()); }
    return r;
}

template <typename MF>
void
MLMGT<MF>::makeSolvable ()
{
    auto const& offset = linop.getSolvabilityOffset(0, 0, rhs[0]);
    if (verbose >= 4) {
        for (int c = 0; c < ncomp; ++c) {
            amrex::Print() << "MLMG: Subtracting " << offset[c] << " from rhs component "
                           << c << "\n";
        }
    }
    for (int alev = 0; alev < namrlevs; ++alev) {
        linop.fixSolvabilityByOffset(alev, 0, rhs[alev], offset);
    }
}

template <typename MF>
void
MLMGT<MF>::makeSolvable (int amrlev, int mglev, MF& mf)
{
    auto const& offset = linop.getSolvabilityOffset(amrlev, mglev, mf);
    if (verbose >= 4) {
        for (int c = 0; c < ncomp; ++c) {
            amrex::Print() << "MLMG: Subtracting " << offset[c]
                           << " from mf component c = " << c
                           << " on level (" << amrlev << ", " << mglev << ")\n";
        }
    }
    linop.fixSolvabilityByOffset(amrlev, mglev, mf, offset);
}

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
template <typename MF>
template <class TMF,std::enable_if_t<std::is_same_v<TMF,MultiFab>,int>>
void
MLMGT<MF>::bottomSolveWithHypre (MF& x, const MF& b)
{
    const int amrlev = 0;
    const int mglev  = linop.NMGLevels(amrlev) - 1;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncomp == 1, "bottomSolveWithHypre doesn't work with ncomp > 1");

    if (linop.isCellCentered())
    {
        if (hypre_solver == nullptr)  // We should reuse the setup
        {
            hypre_solver = linop.makeHypre(hypre_interface);

            hypre_solver->setVerbose(bottom_verbose);
            if (hypre_interface == amrex::Hypre::Interface::ij) {
                hypre_solver->setHypreOptionsNamespace(hypre_options_namespace);
            } else {
                hypre_solver->setHypreOldDefault(hypre_old_default);
                hypre_solver->setHypreRelaxType(hypre_relax_type);
                hypre_solver->setHypreRelaxOrder(hypre_relax_order);
                hypre_solver->setHypreNumSweeps(hypre_num_sweeps);
                hypre_solver->setHypreStrongThreshold(hypre_strong_threshold);
            }

            const BoxArray& ba = linop.m_grids[amrlev].back();
            const DistributionMapping& dm = linop.m_dmap[amrlev].back();
            const Geometry& geom = linop.m_geom[amrlev].back();

            hypre_bndry = std::make_unique<MLMGBndryT<MF>>(ba, dm, ncomp, geom);
            hypre_bndry->setHomogValues();
            const Real* dx = linop.m_geom[0][0].CellSize();
            IntVect crse_ratio = linop.m_coarse_data_crse_ratio.allGT(0) ? linop.m_coarse_data_crse_ratio : IntVect(1);
            RealVect bclocation(AMREX_D_DECL(0.5*dx[0]*crse_ratio[0],
                                             0.5*dx[1]*crse_ratio[1],
                                             0.5*dx[2]*crse_ratio[2]));
            hypre_bndry->setLOBndryConds(linop.m_lobc, linop.m_hibc, IntVect(-1), bclocation,
                                         linop.m_coarse_fine_bc_type);
        }

        // IJ interface understands absolute tolerance API of hypre
        amrex::Real hypre_abstol =
            (hypre_interface == amrex::Hypre::Interface::ij)
            ? bottom_abstol : Real(-1.0);
        hypre_solver->solve(
            x, b, bottom_reltol, hypre_abstol, bottom_maxiter, *hypre_bndry,
            linop.getMaxOrder());
    }
    else
    {
        if (hypre_node_solver == nullptr)
        {
            hypre_node_solver =
                linop.makeHypreNodeLap(bottom_verbose, hypre_options_namespace);
        }
        hypre_node_solver->solve(x, b, bottom_reltol, bottom_abstol, bottom_maxiter);
    }

    // For singular problems there may be a large constant added to all values of the solution
    // For precision reasons we enforce that the average of the correction from hypre is 0
    if (linop.isSingular(amrlev) && linop.getEnforceSingularSolvable())
    {
        makeSolvable(amrlev, mglev, x);
    }
}
#endif

#if defined(AMREX_USE_PETSC) && (AMREX_SPACEDIM > 1)
template <typename MF>
template <class TMF,std::enable_if_t<std::is_same_v<TMF,MultiFab>,int>>
void
MLMGT<MF>::bottomSolveWithPETSc (MF& x, const MF& b)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncomp == 1, "bottomSolveWithPETSc doesn't work with ncomp > 1");

    if(petsc_solver == nullptr)
    {
        petsc_solver = linop.makePETSc();
        petsc_solver->setVerbose(bottom_verbose);

        const BoxArray& ba = linop.m_grids[0].back();
        const DistributionMapping& dm = linop.m_dmap[0].back();
        const Geometry& geom = linop.m_geom[0].back();

        petsc_bndry = std::make_unique<MLMGBndryT<MF>>(ba, dm, ncomp, geom);
        petsc_bndry->setHomogValues();
        const Real* dx = linop.m_geom[0][0].CellSize();
        auto crse_ratio = linop.m_coarse_data_crse_ratio.allGT(0) ? linop.m_coarse_data_crse_ratio : IntVect(1);
        RealVect bclocation(AMREX_D_DECL(0.5*dx[0]*crse_ratio[0],
                                         0.5*dx[1]*crse_ratio[1],
                                         0.5*dx[2]*crse_ratio[2]));
        petsc_bndry->setLOBndryConds(linop.m_lobc, linop.m_hibc, IntVect(-1), bclocation,
                                     linop.m_coarse_fine_bc_type);
    }
    petsc_solver->solve(x, b, bottom_reltol, Real(-1.), bottom_maxiter, *petsc_bndry,
                        linop.getMaxOrder());
}
#endif

template <typename MF>
void
MLMGT<MF>::checkPoint (const Vector<MultiFab*>& a_sol,
                       const Vector<MultiFab const*>& a_rhs,
                       RT a_tol_rel, RT a_tol_abs, const char* a_file_name) const
{
    std::string file_name(a_file_name);
    UtilCreateCleanDirectory(file_name, false);

    if (ParallelContext::IOProcessorSub())
    {
        std::string HeaderFileName(std::string(a_file_name)+"/Header");
        std::ofstream HeaderFile;
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if( ! HeaderFile.good()) {
            FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);

        HeaderFile << linop.name() << "\n"
                   << "a_tol_rel = " << a_tol_rel << "\n"
                   << "a_tol_abs = " << a_tol_abs << "\n"
                   << "verbose = " << verbose << "\n"
                   << "max_iters = " << max_iters << "\n"
                   << "nu1 = " << nu1 << "\n"
                   << "nu2 = " << nu2 << "\n"
                   << "nuf = " << nuf << "\n"
                   << "nub = " << nub << "\n"
                   << "max_fmg_iters = " << max_fmg_iters << "\n"
                   << "bottom_solver = " << static_cast<int>(bottom_solver) << "\n"
                   << "bottom_verbose = " << bottom_verbose << "\n"
                   << "bottom_maxiter = " << bottom_maxiter << "\n"
                   << "bottom_reltol = " << bottom_reltol << "\n"
                   << "always_use_bnorm = " << always_use_bnorm << "\n"
                   << "namrlevs = " << namrlevs << "\n"
                   << "finest_amr_lev = " << finest_amr_lev << "\n"
                   << "linop_prepared = " << linop_prepared << "\n"
                   << "solve_called = " << solve_called << "\n";

        for (int ilev = 0; ilev <= finest_amr_lev; ++ilev) {
            UtilCreateCleanDirectory(file_name+"/Level_"+std::to_string(ilev), false);
        }
    }

    ParallelContext::BarrierSub();

    for (int ilev = 0; ilev <= finest_amr_lev; ++ilev) {
        VisMF::Write(*a_sol[ilev], file_name+"/Level_"+std::to_string(ilev)+"/sol");
        VisMF::Write(*a_rhs[ilev], file_name+"/Level_"+std::to_string(ilev)+"/rhs");
    }

    linop.checkPoint(file_name+"/linop");
}

extern template class MLMGT<MultiFab>;

using MLMG = MLMGT<MultiFab>;

}

#endif
