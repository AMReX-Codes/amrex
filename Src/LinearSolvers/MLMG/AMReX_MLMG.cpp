#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_MLABecLaplacian.H>

#ifdef AMREX_USE_PETSC
#include <petscksp.h>
#include <AMReX_PETSc.H>
#endif

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MLEBABecLap.H>
#endif

// sol: full solution
// rhs: rhs of the original equation L(sol) = rhs
// res: rhs of the residual equation L(cor) = res
// usually res is the result of rhs-L(sol), but is the result of res-L(cor).
// cor/cor_hold: correction
// rescor: res - L(cor)

// x and b as in L(x) = b.  Here x could be either sol or cor.

// BCMode: Inhomogeneous for original equation, Homogeneous for residual equation

// LinOp functions:
//     solutionresidual()  : rhs - L(sol), sol.FillBoundary() will be called.
//                           BC data can be optionally provided.
//     correctionResidual(): res - L(cor), cor.FillBoundary() will be called.
//                           There are BC modes: Homogeneous and Inhomogeneous.
//                           For Inhomogeneous, BC data can be optionally provided.
//     reflux()            : Given sol on crse and fine AMR levels, reflux coarse res at crse/fine.
//     smooth()            : L(cor) = res. cor.FillBoundary() will be called.

namespace amrex {

MLMG::MLMG (MLLinOp& a_lp)
    : linop(a_lp),
      namrlevs(a_lp.NAMRLevels()),
      finest_amr_lev(a_lp.NAMRLevels()-1)
{}

MLMG::~MLMG ()
{}

Real
MLMG::solve (const Vector<MultiFab*>& a_sol, const Vector<MultiFab const*>& a_rhs,
             Real a_tol_rel, Real a_tol_abs, const char* checkpoint_file)
{
    Vector<Any> any_sol(namrlevs);
    Vector<Any> any_rhs(namrlevs);
    for (int lev = 0; lev < namrlevs; ++lev) {
        any_sol[lev] = MultiFab(*a_sol[lev], amrex::make_alias, 0, a_sol[lev]->nComp());
        any_rhs[lev] = MultiFab(*a_rhs[lev], amrex::make_alias, 0, a_rhs[lev]->nComp());
    }
    return solve(any_sol, any_rhs, a_tol_rel, a_tol_abs, checkpoint_file);
}

Real
MLMG::solve (Vector<Any>& a_sol, const Vector<Any>& a_rhs,
             Real a_tol_rel, Real a_tol_abs, const char* checkpoint_file)
{
    BL_PROFILE("MLMG::solve()");

    if (checkpoint_file != nullptr) {
        if (a_sol[0].is<MultiFab>()) {
            Vector<MultiFab*> mf_sol(namrlevs);
            Vector<MultiFab const*> mf_rhs(namrlevs);
            for (int lev = 0; lev < namrlevs; ++lev) {
                mf_sol[lev] = &(a_sol[lev].get<MultiFab>());
                mf_rhs[lev] = &(a_rhs[lev].get<MultiFab>());
            }
            checkPoint(mf_sol, mf_rhs, a_tol_rel, a_tol_abs, checkpoint_file);
        } else {
            amrex::Abort("MLMG::solve: checkpoint not supported for non-MultiFab type");
        }
    }

    if (bottom_solver == BottomSolver::Default) {
        bottom_solver = linop.getDefaultBottomSolver();
    }

#if defined(AMREX_USE_HYPRE) || defined(AMREX_USE_PETSC)
    if (bottom_solver == BottomSolver::hypre || bottom_solver == BottomSolver::petsc) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_sol[0].is<MultiFab>(),
                                         "Non-MultiFab type not supported for hypre and petsc");
        int mo = linop.getMaxOrder();
        if (a_sol[0].get<MultiFab>().hasEBFabFactory()) {
            linop.setMaxOrder(2);
        } else {
            linop.setMaxOrder(std::min(3,mo));  // maxorder = 4 not supported
        }
    }
#endif

    bool is_nsolve = linop.m_parent;

    auto solve_start_time = amrex::second();

    Real& composite_norminf = m_final_resnorm0;

    m_niters_cg.clear();
    m_iter_fine_resnorm0.clear();

    prepareForSolve(a_sol, a_rhs);

    computeMLResidual(finest_amr_lev);

    bool local = true;
    Real resnorm0 = MLResNormInf(finest_amr_lev, local);
    Real rhsnorm0 = MLRhsNormInf(local);
    if (!is_nsolve) {
        ParallelAllReduce::Max<Real>({resnorm0, rhsnorm0}, ParallelContext::CommunicatorSub());

        if (verbose >= 1)
        {
            amrex::Print() << "MLMG: Initial rhs               = " << rhsnorm0 << "\n"
                           << "MLMG: Initial residual (resid0) = " << resnorm0 << "\n";
        }
    }

    m_init_resnorm0 = resnorm0;
    m_rhsnorm0 = rhsnorm0;

    Real max_norm;
    std::string norm_name;
    if (always_use_bnorm || rhsnorm0 >= resnorm0) {
        norm_name = "bnorm";
        max_norm = rhsnorm0;
    } else {
        norm_name = "resid0";
        max_norm = resnorm0;
    }
    const Real res_target = std::max(a_tol_abs, std::max(a_tol_rel,Real(1.e-16))*max_norm);

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

            if (is_nsolve) continue;

            Real fine_norminf = ResNormInf(finest_amr_lev);
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
                Real crse_norminf = MLResNormInf(finest_amr_lev-1);
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
              if (composite_norminf > Real(1.e20)*max_norm)
              {
                  if (verbose > 0) {
                      amrex::Print() << "MLMG: Failing to converge after " << iter+1 << " iterations."
                                     << " resid, resid/" << norm_name << " = "
                                     << composite_norminf << ", "
                                     << composite_norminf/max_norm << "\n";
                  }
                  amrex::Abort("MLMG failing so lets stop here");
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
            amrex::Abort("MLMG failed");
        }
        timer[iter_time] = amrex::second() - iter_start_time;
    }

    IntVect ng_back = final_fill_bc ? IntVect(1) : IntVect(0);
    if (linop.hasHiddenDimension()) {
        ng_back[linop.hiddenDirection()] = 0;
    }
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (!sol_is_alias[alev]) {
            linop.AnyCopy(a_sol[alev], sol[alev], ng_back);
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

// in  : Residual (res) on the finest AMR level
// out : sol on all AMR levels
void MLMG::oneIter (int iter)
{
    BL_PROFILE("MLMG::oneIter()");

    for (int alev = finest_amr_lev; alev > 0; --alev)
    {
        miniCycle(alev);

        IntVect nghost(0);
        if (cf_strategy == CFStrategy::ghostnodes) nghost = IntVect(linop.getNGrow(alev));
        linop.AnyAdd(sol[alev], cor[alev][0], nghost);

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
        if (cf_strategy == CFStrategy::ghostnodes) nghost = IntVect(linop.getNGrow(0));
        linop.AnyAdd(sol[0], cor[0][0], nghost);
    }

    for (int alev = 1; alev <= finest_amr_lev; ++alev)
    {
        // (Fine AMR correction) = I(Coarse AMR correction)
        interpCorrection(alev);

        IntVect nghost(0);
        if (cf_strategy == CFStrategy::ghostnodes) nghost = IntVect(linop.getNGrow(alev));
        linop.AnyAdd(sol[alev], cor[alev][0], nghost);

        if (alev != finest_amr_lev) {
            linop.AnyAdd(cor_hold[alev][0], cor[alev][0], nghost);
        }

        // Update fine AMR level correction
        computeResWithCrseCorFineCor(alev);

        miniCycle(alev);

        linop.AnyAdd(sol[alev], cor[alev][0], nghost);

        if (alev != finest_amr_lev) {
            linop.AnyAdd(cor[alev][0], cor_hold[alev][0], nghost);
        }
    }

    linop.AnyAverageDownAndSync(sol);
}

// Compute multi-level Residual (res) up to amrlevmax.
void
MLMG::computeMLResidual (int amrlevmax)
{
    BL_PROFILE("MLMG::computeMLResidual()");

    const int mglev = 0;
    for (int alev = amrlevmax; alev >= 0; --alev) {
        const Any* crse_bcdata = (alev > 0) ? &(sol[alev-1]) : nullptr;
        linop.AnySolutionResidual(alev, res[alev][mglev], sol[alev], rhs[alev], crse_bcdata);
        if (alev < finest_amr_lev) {
            linop.AnyReflux(alev, res[alev][mglev], sol[alev], rhs[alev],
                            res[alev+1][mglev], sol[alev+1], rhs[alev+1]);
        }
    }
}

// Compute single AMR level residual without masking.
void
MLMG::computeResidual (int alev)
{
    BL_PROFILE("MLMG::computeResidual()");
    const Any* crse_bcdata = (alev > 0) ? &(sol[alev-1]) : nullptr;
    linop.AnySolutionResidual(alev, res[alev][0], sol[alev], rhs[alev], crse_bcdata);
}

// Compute coarse AMR level composite residual with coarse solution and fine correction
void
MLMG::computeResWithCrseSolFineCor (int calev, int falev)
{
    BL_PROFILE("MLMG::computeResWithCrseSolFineCor()");

    IntVect nghost(0);
    if (cf_strategy == CFStrategy::ghostnodes) nghost = IntVect(std::min(linop.getNGrow(falev),linop.getNGrow(calev)));

    Any&       crse_sol = sol[calev];
    const Any& crse_rhs = rhs[calev];
    Any&       crse_res = res[calev][0];

    Any&       fine_sol = sol[falev];
    const Any& fine_rhs = rhs[falev];
    Any&       fine_cor = cor[falev][0];
    Any&       fine_res = res[falev][0];
    Any&    fine_rescor = rescor[falev][0];

    const Any* crse_bcdata = (calev > 0) ? &(sol[calev-1]) : nullptr;
    linop.AnySolutionResidual(calev, crse_res, crse_sol, crse_rhs, crse_bcdata);

    linop.AnyCorrectionResidual(falev, 0, fine_rescor, fine_cor, fine_res, BCMode::Homogeneous);
    linop.AnyCopy(fine_res, fine_rescor, nghost);

    linop.AnyReflux(calev, crse_res, crse_sol, crse_rhs, fine_res, fine_sol, fine_rhs);

    linop.AnyAvgDownResAmr(calev, crse_res, fine_res);
}

// Compute fine AMR level residual fine_res = fine_res - L(fine_cor) with coarse providing BC.
void
MLMG::computeResWithCrseCorFineCor (int falev)
{
    BL_PROFILE("MLMG::computeResWithCrseCorFineCor()");

    IntVect nghost(0);
    if (cf_strategy == CFStrategy::ghostnodes) nghost = IntVect(linop.getNGrow(falev));

    const Any& crse_cor = cor[falev-1][0];

    Any& fine_cor    = cor   [falev][0];
    Any& fine_res    = res   [falev][0];
    Any& fine_rescor = rescor[falev][0];

    // fine_rescor = fine_res - L(fine_cor)
    linop.AnyCorrectionResidual(falev, 0, fine_rescor, fine_cor, fine_res,
                                BCMode::Inhomogeneous, &crse_cor);
    linop.AnyCopy(fine_res, fine_rescor, nghost);
}

void
MLMG::miniCycle (int amrlev)
{
    BL_PROFILE("MLMG::miniCycle()");
    const int mglev = 0;
    mgVcycle(amrlev, mglev);
}

// in   : Residual (res)
// out  : Correction (cor) from bottom to this function's local top
void
MLMG::mgVcycle (int amrlev, int mglev_top)
{
    BL_PROFILE("MLMG::mgVcycle()");

    const int mglev_bottom = linop.NMGLevels(amrlev) - 1;

    for (int mglev = mglev_top; mglev < mglev_bottom; ++mglev)
    {
        BL_PROFILE_VAR("MLMG::mgVcycle_down::"+std::to_string(mglev), blp_mgv_down_lev);

        if (verbose >= 4)
        {
            Real norm = linop.AnyNormInf(res[amrlev][mglev]);
            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev
                           << "   DN: Norm before smooth " << norm << "\n";
        }

        linop.AnySetToZero(cor[amrlev][mglev]);
        bool skip_fillboundary = true;
        for (int i = 0; i < nu1; ++i) {
            linop.AnySmooth(amrlev, mglev, cor[amrlev][mglev], res[amrlev][mglev],
                            skip_fillboundary);
            skip_fillboundary = false;
        }

        // rescor = res - L(cor)
        computeResOfCorrection(amrlev, mglev);

        if (verbose >= 4)
        {
            Real norm = linop.AnyNormInf(rescor[amrlev][mglev]);
            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev
                           << "   DN: Norm after  smooth " << norm << "\n";
        }

        // res_crse = R(rescor_fine); this provides res/b to the level below
        linop.AnyRestriction(amrlev, mglev+1, res[amrlev][mglev+1], rescor[amrlev][mglev]);
    }

    BL_PROFILE_VAR("MLMG::mgVcycle_bottom", blp_bottom);
    if (amrlev == 0)
    {
        if (verbose >= 4)
        {
            Real norm = linop.AnyNormInf(res[amrlev][mglev_bottom]);
            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev_bottom
                           << "   DN: Norm before bottom " << norm << "\n";
        }
        bottomSolve();
        if (verbose >= 4)
        {
            computeResOfCorrection(amrlev, mglev_bottom);
            Real norm = linop.AnyNormInf(rescor[amrlev][mglev_bottom]);

            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev_bottom
                           << "   UP: Norm after  bottom " << norm << "\n";
        }
    }
    else
    {
        if (verbose >= 4)
        {
            Real norm = linop.AnyNormInf(res[amrlev][mglev_bottom]);
            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev_bottom
                           << "       Norm before smooth " << norm << "\n";
        }
        linop.AnySetToZero(cor[amrlev][mglev_bottom]);
        bool skip_fillboundary = true;
        for (int i = 0; i < nu1; ++i) {
            linop.AnySmooth(amrlev, mglev_bottom, cor[amrlev][mglev_bottom],
                            res[amrlev][mglev_bottom], skip_fillboundary);
            skip_fillboundary = false;
        }
        if (verbose >= 4)
        {
            computeResOfCorrection(amrlev, mglev_bottom);
            Real norm = linop.AnyNormInf(rescor[amrlev][mglev_bottom]);
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
            Real norm = linop.AnyNormInf(rescor[amrlev][mglev]);
            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev
                           << "   UP: Norm before smooth " << norm << "\n";
        }
        for (int i = 0; i < nu2; ++i) {
            linop.AnySmooth(amrlev, mglev, cor[amrlev][mglev], res[amrlev][mglev]);
        }

        if (cf_strategy == CFStrategy::ghostnodes) computeResOfCorrection(amrlev, mglev);

        if (verbose >= 4)
        {
            computeResOfCorrection(amrlev, mglev);
            Real norm = linop.AnyNormInf(rescor[amrlev][mglev]);
            amrex::Print() << "AT LEVEL "  << amrlev << " " << mglev
                           << "   UP: Norm after  smooth " << norm << "\n";
        }
    }
}

// FMG cycle on the coarsest AMR level.
// in:  Residual on the top MG level (i.e., 0)
// out: Correction (cor) on all MG levels
void
MLMG::mgFcycle ()
{
    BL_PROFILE("MLMG::mgFcycle()");

#ifdef AMREX_USE_EB
    AMREX_ASSERT(linop.isCellCentered());
#endif

    const int amrlev = 0;
    const int mg_bottom_lev = linop.NMGLevels(amrlev) - 1;
    IntVect nghost(0);
    if (cf_strategy == CFStrategy::ghostnodes) nghost = IntVect(linop.getNGrow(amrlev));

    for (int mglev = 1; mglev <= mg_bottom_lev; ++mglev)
    {
        linop.AnyAvgDownResMG(mglev, res[amrlev][mglev], res[amrlev][mglev-1]);
    }

    bottomSolve();

    for (int mglev = mg_bottom_lev-1; mglev >= 0; --mglev)
    {
        // cor_fine = I(cor_crse)
        interpCorrection(amrlev, mglev);

        // rescor = res - L(cor)
        computeResOfCorrection(amrlev, mglev);
        // res = rescor; this provides b to the vcycle below
        linop.AnyCopy(res[amrlev][mglev], rescor[amrlev][mglev], nghost);

        // save cor; do v-cycle; add the saved to cor
        std::swap(cor[amrlev][mglev], cor_hold[amrlev][mglev]);
        mgVcycle(amrlev, mglev);
        linop.AnyAdd(cor[amrlev][mglev], cor_hold[amrlev][mglev], nghost);
    }
}

// Interpolate correction from coarse to fine AMR level.
void
MLMG::interpCorrection (int alev)
{
    BL_PROFILE("MLMG::interpCorrection_1");

    IntVect nghost(0);
    if (cf_strategy == CFStrategy::ghostnodes) nghost = IntVect(linop.getNGrow(alev));

    Any const& crse_cor = cor[alev-1][0];
    Any      & fine_cor = cor[alev  ][0];

    const Geometry& crse_geom = linop.Geom(alev-1,0);

    int ng_src = 0;
    int ng_dst = linop.isCellCentered() ? 1 : 0;
    if (cf_strategy == CFStrategy::ghostnodes)
    {
        ng_src = linop.getNGrow(alev-1);
        ng_dst = linop.getNGrow(alev-1);
    }

    Any cfine = linop.AnyMakeCoarseAmr(alev, IntVect(ng_dst));
    linop.AnySetToZero(cfine);
    linop.AnyParallelCopy(cfine, crse_cor, IntVect(ng_src), IntVect(ng_dst), crse_geom.periodicity());

    linop.AnyInterpolationAmr(alev, fine_cor, cfine, nghost);
}

// Interpolate correction between MG levels
// inout: Correction (cor) on coarse MG lev.  (out due to FillBoundary)
// out  : Correction (cor) on fine MG lev.
void
MLMG::interpCorrection (int alev, int mglev)
{
    BL_PROFILE("MLMG::interpCorrection_2");

    Any& crse_cor = cor[alev][mglev+1];
    Any& fine_cor = cor[alev][mglev  ];
    linop.AnyInterpAssignMG(alev, mglev, fine_cor, crse_cor);
}

// (Fine MG level correction) += I(Coarse MG level correction)
void
MLMG::addInterpCorrection (int alev, int mglev)
{
    BL_PROFILE("MLMG::addInterpCorrection()");

    const Any& crse_cor = cor[alev][mglev+1];
    Any&       fine_cor = cor[alev][mglev  ];

    Any cfine;
    const Any* cany;

    if (linop.isMFIterSafe(alev, mglev, mglev+1))
    {
        cany = &crse_cor;
    }
    else
    {
        cfine = linop.AnyMakeCoarseMG(alev, mglev, IntVect(0));
        linop.AnyParallelCopy(cfine,crse_cor,IntVect(0),IntVect(0));
        cany = &cfine;
    }

    linop.AnyInterpolationMG(alev, mglev, fine_cor, *cany);
}

// Compute rescor = res - L(cor)
// in   : res
// inout: cor (out due to FillBoundary in linop.correctionResidual)
// out  : rescor
void
MLMG::computeResOfCorrection (int amrlev, int mglev)
{
    BL_PROFILE("MLMG:computeResOfCorrection()");
    Any      & x =    cor[amrlev][mglev];
    const Any& b =    res[amrlev][mglev];
    Any      & r = rescor[amrlev][mglev];
    linop.AnyCorrectionResidual(amrlev, mglev, r, x, b, BCMode::Homogeneous);
}

// At the true bottom of the coarset AMR level.
// in  : Residual (res) as b
// out : Correction (cor) as x
void
MLMG::bottomSolve ()
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

void
MLMG::NSolve (MLMG& a_solver, MultiFab& a_sol, MultiFab& a_rhs)
{
    BL_PROFILE("MLMG::NSolve()");

    a_sol.setVal(0.0);

    MultiFab const& res_bottom = res[0].back().get<MultiFab>();
    if (BoxArray::SameRefs(a_rhs.boxArray(),res_bottom.boxArray()) &&
        DistributionMapping::SameRefs(a_rhs.DistributionMap(),res_bottom.DistributionMap()))
    {
        MultiFab::Copy(a_rhs, res_bottom, 0, 0, a_rhs.nComp(), 0);
    } else {
        a_rhs.setVal(0.0);
        a_rhs.ParallelCopy(res_bottom);
    }

    a_solver.solve({&a_sol}, {&a_rhs}, Real(-1.0), Real(-1.0));

    linop.copyNSolveSolution(cor[0].back().get<MultiFab>(), a_sol);
}

void
MLMG::actualBottomSolve ()
{
    BL_PROFILE("MLMG::actualBottomSolve()");

    if (!linop.isBottomActive()) return;

    auto bottom_start_time = amrex::second();

    ParallelContext::push(linop.BottomCommunicator());

    const int amrlev = 0;
    const int mglev = linop.NMGLevels(amrlev) - 1;
    auto& x = cor[amrlev][mglev];
    auto& b = res[amrlev][mglev];

    linop.AnySetToZero(x);

    if (bottom_solver == BottomSolver::smoother)
    {
        bool skip_fillboundary = true;
        for (int i = 0; i < nuf; ++i) {
            linop.AnySmooth(amrlev, mglev, x, b, skip_fillboundary);
            skip_fillboundary = false;
        }
    }
    else
    {
        Any* bottom_b = &b;
        Any raii_b;
        if (linop.isBottomSingular() && linop.getEnforceSingularSolvable())
        {
            const IntVect ng = linop.AnyGrowVect(b);
            raii_b = linop.AnyMake(amrlev, mglev, ng);
            linop.AnyCopy(raii_b, b, ng);
            bottom_b = &raii_b;

            makeSolvable(amrlev,mglev,*bottom_b);
        }

        if (bottom_solver == BottomSolver::hypre)
        {
#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
            bottomSolveWithHypre(x, *bottom_b);
#endif
        }
        else if (bottom_solver == BottomSolver::petsc)
        {
            bottomSolveWithPETSc(x, *bottom_b);
        }
        else
        {
            MLCGSolver::Type cg_type;
            if (bottom_solver == BottomSolver::cg ||
                bottom_solver == BottomSolver::cgbicg) {
                cg_type = MLCGSolver::Type::CG;
            } else {
                cg_type = MLCGSolver::Type::BiCGStab;
            }
            int ret = bottomSolveWithCG(x, *bottom_b, cg_type);
            // If the MLMG solve failed then set the correction to zero
            if (ret != 0) {
                linop.AnySetToZero(cor[amrlev][mglev]);
                if (bottom_solver == BottomSolver::cgbicg ||
                    bottom_solver == BottomSolver::bicgcg) {
                    if (bottom_solver == BottomSolver::cgbicg) {
                        cg_type = MLCGSolver::Type::BiCGStab; // switch to bicg
                    } else {
                        cg_type = MLCGSolver::Type::CG; // switch to cg
                    }
                    ret = bottomSolveWithCG(x, *bottom_b, cg_type);
                    if (ret != 0) {
                        linop.AnySetToZero(cor[amrlev][mglev]);
                    } else { // switch permanently
                        if (cg_type == MLCGSolver::Type::CG) {
                            bottom_solver = BottomSolver::cg;
                        } else {
                            bottom_solver = BottomSolver::bicgstab;
                        }
                    }
                }
            }
            const int n = (ret==0) ? nub : nuf;
            for (int i = 0; i < n; ++i) {
                linop.AnySmooth(amrlev, mglev, x, b);
            }
        }
    }

    ParallelContext::pop();

    timer[bottom_time] += amrex::second() - bottom_start_time;
}

int
MLMG::bottomSolveWithCG (Any& x, const Any& b, MLCGSolver::Type type)
{
    MLCGSolver cg_solver(this, linop);
    cg_solver.setSolver(type);
    cg_solver.setVerbose(bottom_verbose);
    cg_solver.setMaxIter(bottom_maxiter);
    if (cf_strategy == CFStrategy::ghostnodes) cg_solver.setNGhost(linop.getNGrow());

    int ret = cg_solver.solve(x, b, bottom_reltol, bottom_abstol);
    if (ret != 0 && verbose > 1) {
        amrex::Print() << "MLMG: Bottom solve failed.\n";
    }
    m_niters_cg.push_back(cg_solver.getNumIters());
    return ret;
}

// Compute single-level masked inf-norm of Residual (res).
Real
MLMG::ResNormInf (int alev, bool local)
{
    BL_PROFILE("MLMG::ResNormInf()");
    return linop.AnyNormInfMask(alev, res[alev][0], local);
}

// Computes multi-level masked inf-norm of Residual (res).
Real
MLMG::MLResNormInf (int alevmax, bool local)
{
    BL_PROFILE("MLMG::MLResNormInf()");
    Real r = 0.0;
    for (int alev = 0; alev <= alevmax; ++alev)
    {
        r = std::max(r, ResNormInf(alev,true));
    }
    if (!local) ParallelAllReduce::Max(r, ParallelContext::CommunicatorSub());
    return r;
}

// Compute multi-level masked inf-norm of RHS (rhs).
Real
MLMG::MLRhsNormInf (bool local)
{
    BL_PROFILE("MLMG::MLRhsNormInf()");
    Real r = 0.0_rt;
    for (int alev = 0; alev <= finest_amr_lev; ++alev) {
        auto t = linop.AnyNormInfMask(alev, rhs[alev], true);
        r = std::max(r, t);
    }
    if (!local) ParallelAllReduce::Max(r, ParallelContext::CommunicatorSub());
    return r;
}

void
MLMG::prepareForSolve (Vector<Any>& a_sol, const Vector<Any>& a_rhs)
{
    BL_PROFILE("MLMG::prepareForSolve()");

    AMREX_ASSERT(namrlevs <= a_sol.size());
    AMREX_ASSERT(namrlevs <= a_rhs.size());

    timer.assign(ntimers, 0.0);

    IntVect ng_rhs(0);
    IntVect ng_sol(1);
    if (linop.hasHiddenDimension()) ng_sol[linop.hiddenDirection()] = 0;

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
    sol_is_alias.resize(namrlevs);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (cf_strategy == CFStrategy::ghostnodes)
        {
            sol[alev] = linop.AnyMakeAlias(a_sol[alev]);
            sol_is_alias[alev] = true;
        }
        else if (linop.AnyGrowVect(a_sol[alev]) == ng_sol)
        {
            sol[alev] = linop.AnyMakeAlias(a_sol[alev]);
            linop.AnySetBndryToZero(sol[alev]);
            sol_is_alias[alev] = true;
        }
        else
        {
            if (!solve_called) {
                sol[alev] = linop.AnyMake(alev, 0, ng_sol);
            }
            linop.AnyCopy(sol[alev], a_sol[alev], IntVect(0));
            linop.AnySetBndryToZero(sol[alev]);
            sol_is_alias[alev] = false;
        }
    }

    rhs.resize(namrlevs);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (cf_strategy == CFStrategy::ghostnodes) ng_rhs = IntVect(linop.getNGrow(alev));
        if (!solve_called) {
            rhs[alev] = linop.AnyMake(alev, 0, ng_rhs);
        }
        linop.AnyCopy(rhs[alev], a_rhs[alev], ng_rhs);
        linop.applyMetricTerm(alev, 0, rhs[alev]);
        linop.unimposeNeumannBC(alev, rhs[alev]);
        linop.applyInhomogNeumannTerm(alev, rhs[alev]);
        linop.applyOverset(alev, rhs[alev]);
        linop.scaleRHS(alev, rhs[alev]);

#ifdef AMREX_USE_EB
        auto factory = dynamic_cast<EBFArrayBoxFactory const*>(linop.Factory(alev));
        if (factory) {
            linop.AnySetCoveredToZero(rhs[alev]);
            linop.AnySetCoveredToZero(sol[alev]);
        }
#endif
    }

    for (int falev = finest_amr_lev; falev > 0; --falev)
    {
        linop.AnyAverageDownSolutionRHS(falev-1, sol[falev-1], rhs[falev-1],
                                        sol[falev], rhs[falev]);
    }

    // enforce solvability if appropriate
    if (linop.isSingular(0) && linop.getEnforceSingularSolvable())
    {
        makeSolvable();
    }

    IntVect ng = linop.isCellCentered() ? IntVect(0) : IntVect(1);
    if (cf_strategy == CFStrategy::ghostnodes) ng = ng_rhs;
    if (!solve_called) {
        linop.make(res, ng);
        linop.make(rescor, ng);
    }
    for (int alev = 0; alev <= finest_amr_lev; ++alev)
    {
        const int nmglevs = linop.NMGLevels(alev);
        for (int mglev = 0; mglev < nmglevs; ++mglev)
        {
            linop.AnySetToZero(res   [alev][mglev]);
            linop.AnySetToZero(rescor[alev][mglev]);
        }
    }

    if (cf_strategy != CFStrategy::ghostnodes) ng = ng_sol;
    cor.resize(namrlevs);
    for (int alev = 0; alev <= finest_amr_lev; ++alev)
    {
        const int nmglevs = linop.NMGLevels(alev);
        cor[alev].resize(nmglevs);
        for (int mglev = 0; mglev < nmglevs; ++mglev)
        {
            if (!solve_called) {
                IntVect _ng = ng;
                if (cf_strategy == CFStrategy::ghostnodes) _ng=IntVect(linop.getNGrow(alev,mglev));
                cor[alev][mglev] = linop.AnyMake(alev, mglev, _ng);
            }
            linop.AnySetToZero(cor[alev][mglev]);
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
                if (cf_strategy == CFStrategy::ghostnodes) _ng=IntVect(linop.getNGrow(alev,mglev));
                cor_hold[alev][mglev] = linop.AnyMake(alev, mglev, _ng);
            }
            linop.AnySetToZero(cor_hold[alev][mglev]);
        }
    }
    for (int alev = 1; alev < finest_amr_lev; ++alev)
    {
        cor_hold[alev].resize(1);
        if (!solve_called) {
            IntVect _ng = ng;
            if (cf_strategy == CFStrategy::ghostnodes) _ng=IntVect(linop.getNGrow(alev));
            cor_hold[alev][0] = linop.AnyMake(alev, 0, _ng);
        }
        linop.AnySetToZero(cor_hold[alev][0]);
    }

    if (linop.m_parent) {
        do_nsolve = false;  // no embedded N-Solve
    } else if (!linop.supportNSolve()) {
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

void
MLMG::prepareForNSolve ()
{
    ns_linop = linop.makeNLinOp(nsolve_grid_size);

    const int ncomp = linop.getNComp();
    int nghost = 0;
    if (cf_strategy == CFStrategy::ghostnodes) nghost = linop.getNGrow();

    const BoxArray& ba = (*ns_linop).m_grids[0][0];
    const DistributionMapping& dm =(*ns_linop).m_dmap[0][0];

    int ng = 1;
    if (cf_strategy == CFStrategy::ghostnodes) ng = nghost;
    ns_sol = std::make_unique<MultiFab>(ba, dm, ncomp, ng, MFInfo(), *(ns_linop->Factory(0,0)));
    ng = 0;
    if (cf_strategy == CFStrategy::ghostnodes) ng = nghost;
    ns_rhs = std::make_unique<MultiFab>(ba, dm, ncomp, ng, MFInfo(), *(ns_linop->Factory(0,0)));
    ns_sol->setVal(0.0);
    ns_rhs->setVal(0.0);

    ns_linop->setLevelBC(0, ns_sol.get());

    ns_mlmg = std::make_unique<MLMG>(*ns_linop);
    ns_mlmg->setVerbose(0);
    ns_mlmg->setFixedIter(1);
    ns_mlmg->setMaxFmgIter(20);
    ns_mlmg->setBottomSolver(BottomSolver::smoother);
}

void
MLMG::getGradSolution (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_grad_sol,
                       Location a_loc)
{
    BL_PROFILE("MLMG::getGradSolution()");
    for (int alev = 0; alev <= finest_amr_lev; ++alev) {
        linop.compGrad(alev, a_grad_sol[alev], sol[alev].get<MultiFab>(), a_loc);
    }
}

void
MLMG::getFluxes (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_flux,
                 Location a_loc)
{
    if (!linop.isCellCentered()) {
       amrex::Abort("Calling wrong getFluxes for nodal solver");
    }

    AMREX_ASSERT(sol.size() == a_flux.size());
    Vector<MultiFab*> solmf;
    for (auto & s : sol) {
        solmf.push_back(&(s.get<MultiFab>()));
    }
    getFluxes(a_flux, solmf, a_loc);
}

void
MLMG::getFluxes (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_flux,
                 const Vector<MultiFab*>& a_sol,
                 Location a_loc)
{
    BL_PROFILE("MLMG::getFluxes()");

    if (!linop.isCellCentered()) {
       amrex::Abort("Calling wrong getFluxes for nodal solver");
    }

    linop.getFluxes(a_flux, a_sol, a_loc);
}

void
MLMG::getFluxes (const Vector<MultiFab*> & a_flux, Location a_loc)
{
    AMREX_ASSERT(sol.size() == a_flux.size());
    Vector<MultiFab*> solmf;
    for (auto & s : sol) {
        solmf.push_back(&(s.get<MultiFab>()));
    }
    getFluxes(a_flux, solmf, a_loc);
}

void
MLMG::getFluxes (const Vector<MultiFab*> & a_flux, const Vector<MultiFab*>& a_sol, Location /*a_loc*/)
{
    AMREX_ASSERT(a_flux[0]->nComp() >= AMREX_SPACEDIM);

    if (linop.isCellCentered())
    {
        Vector<Array<MultiFab,AMREX_SPACEDIM> > ffluxes(namrlevs);
        for (int alev = 0; alev < namrlevs; ++alev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const int mglev = 0;
                const int ncomp = linop.getNComp();
                int nghost = 0;
                if (cf_strategy == CFStrategy::ghostnodes) nghost = linop.getNGrow(alev);
                ffluxes[alev][idim].define(amrex::convert(linop.m_grids[alev][mglev],
                                                          IntVect::TheDimensionVector(idim)),
                                           linop.m_dmap[alev][mglev], ncomp, nghost, MFInfo(),
                                           *linop.m_factory[alev][mglev]);
            }
        }
        getFluxes(amrex::GetVecOfArrOfPtrs(ffluxes), a_sol, Location::FaceCenter);
        for (int alev = 0; alev < namrlevs; ++alev) {
#ifdef AMREX_USE_EB
            EB_average_face_to_cellcenter(*a_flux[alev], 0, amrex::GetArrOfConstPtrs(ffluxes[alev]));
#else
            average_face_to_cellcenter(*a_flux[alev], 0, amrex::GetArrOfConstPtrs(ffluxes[alev]));
#endif
        }

    } else {
        linop.getFluxes(a_flux, a_sol);
    }
}

#ifdef AMREX_USE_EB
void
MLMG::getEBFluxes (const Vector<MultiFab*>& a_eb_flux)
{
    if (!linop.isCellCentered()) {
       amrex::Abort("getEBFluxes is for cell-centered only");
    }

    AMREX_ASSERT(sol.size() == a_eb_flux.size());
    Vector<MultiFab*> solmf;
    for (auto & s : sol) {
        solmf.push_back(&(s.get<MultiFab>()));
    }
    getEBFluxes(a_eb_flux, solmf);
}

void
MLMG::getEBFluxes (const Vector<MultiFab*>& a_eb_flux, const Vector<MultiFab*>& a_sol)
{
    BL_PROFILE("MLMG::getEBFluxes()");

    if (!linop.isCellCentered()) {
       amrex::Abort("getEBFluxes is for cell-centered only");
    }

    linop.getEBFluxes(a_eb_flux, a_sol);
}
#endif

void
MLMG::compResidual (const Vector<MultiFab*>& a_res, const Vector<MultiFab*>& a_sol,
                    const Vector<MultiFab const*>& a_rhs)
{
    BL_PROFILE("MLMG::compResidual()");

    const int ncomp = linop.getNComp();
    IntVect ng_sol(1);
    if (linop.hasHiddenDimension()) ng_sol[linop.hiddenDirection()] = 0;

    sol.resize(namrlevs);
    sol_is_alias.resize(namrlevs,true);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (cf_strategy == CFStrategy::ghostnodes || a_sol[alev]->nGrowVect() == ng_sol)
        {
            sol[alev] = linop.AnyMakeAlias(*a_sol[alev]);
            sol_is_alias[alev] = true;
        }
        else
        {
            if (sol_is_alias[alev])
            {
                sol[alev] = linop.AnyMake(alev, 0, ng_sol);
            }
            MultiFab::Copy(sol[alev].get<MultiFab>(), *a_sol[alev], 0, 0, ncomp, 0);
        }
    }

    if (!linop_prepared) {
        linop.prepareForSolve();
        linop_prepared = true;
    } else if (linop.needsUpdate()) {
        linop.update();
    }

    const auto& amrrr = linop.AMRRefRatio();

    for (int alev = finest_amr_lev; alev >= 0; --alev) {
        const MultiFab* crse_bcdata = (alev > 0) ? &(sol[alev-1].get<MultiFab>()) : nullptr;
        const MultiFab* prhs = a_rhs[alev];
#if (AMREX_SPACEDIM != 3)
        int nghost = (cf_strategy == CFStrategy::ghostnodes) ? linop.getNGrow(alev) : 0;
        Any rhstmp_a(MultiFab(prhs->boxArray(), prhs->DistributionMap(), ncomp, nghost,
                              MFInfo(), *linop.Factory(alev)));
        MultiFab& rhstmp = rhstmp_a.get<MultiFab>();
        MultiFab::Copy(rhstmp, *prhs, 0, 0, ncomp, nghost);
        linop.applyMetricTerm(alev, 0, rhstmp_a);
        linop.unimposeNeumannBC(alev, rhstmp_a);
        linop.applyInhomogNeumannTerm(alev, rhstmp_a);
        prhs = &rhstmp;
#endif
        linop.solutionResidual(alev, *a_res[alev], sol[alev].get<MultiFab>(), *prhs, crse_bcdata);
        if (alev < finest_amr_lev) {
            linop.reflux(alev, *a_res[alev], sol[alev].get<MultiFab>(), *prhs,
                         *a_res[alev+1], sol[alev+1].get<MultiFab>(), *a_rhs[alev+1]);
            if (linop.isCellCentered()) {
#ifdef AMREX_USE_EB
                amrex::EB_average_down(*a_res[alev+1], *a_res[alev], 0, ncomp, amrrr[alev]);
#else
                amrex::average_down(*a_res[alev+1], *a_res[alev], 0, ncomp, amrrr[alev]);
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

void
MLMG::apply (const Vector<MultiFab*>& out, const Vector<MultiFab*>& a_in)
{
    BL_PROFILE("MLMG::apply()");

    Vector<MultiFab*> in(namrlevs);
    Vector<MultiFab> in_raii(namrlevs);
    Vector<MultiFab> rh(namrlevs);
    int nghost = 0;
    IntVect ng_sol(1);
    if (linop.hasHiddenDimension()) ng_sol[linop.hiddenDirection()] = 0;

    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (cf_strategy == CFStrategy::ghostnodes)
        {
            nghost = linop.getNGrow(alev);
            in[alev] = a_in[alev];
        }
        else if (a_in[alev]->nGrowVect() == ng_sol)
        {
            in[alev] = a_in[alev];
        }
        else
        {
            IntVect ng = ng_sol;
            if (cf_strategy == CFStrategy::ghostnodes) ng = IntVect(nghost);
            in_raii[alev].define(a_in[alev]->boxArray(),
                                 a_in[alev]->DistributionMap(),
                                 a_in[alev]->nComp(), ng,
                                 MFInfo(), *linop.Factory(alev));
            MultiFab::Copy(in_raii[alev], *a_in[alev], 0, 0, a_in[alev]->nComp(), nghost);
            in[alev] = &(in_raii[alev]);
        }
        rh[alev].define(a_in[alev]->boxArray(),
                        a_in[alev]->DistributionMap(),
                        a_in[alev]->nComp(), nghost, MFInfo(),
                        *linop.Factory(alev));
        rh[alev].setVal(0.0);
    }

    if (!linop_prepared) {
        linop.prepareForSolve();
        linop_prepared = true;
    } else if (linop.needsUpdate()) {
        linop.update();
    }

    for (int alev = 0; alev < namrlevs; ++alev) {
        Any a(MultiFab(rh[alev], amrex::make_alias, 0, rh[alev].nComp()));
        linop.applyInhomogNeumannTerm(alev, a);
    }

    const auto& amrrr = linop.AMRRefRatio();

    for (int alev = finest_amr_lev; alev >= 0; --alev) {
        const MultiFab* crse_bcdata = (alev > 0) ? in[alev-1] : nullptr;
        linop.solutionResidual(alev, *out[alev], *in[alev], rh[alev], crse_bcdata);
        if (alev < finest_amr_lev) {
            linop.reflux(alev, *out[alev], *in[alev], rh[alev],
                         *out[alev+1], *in[alev+1], rh[alev+1]);
            if (linop.isCellCentered()) {
#ifdef AMREX_USE_EB
                amrex::EB_average_down(*out[alev+1], *out[alev], 0, out[alev]->nComp(), amrrr[alev]);
#else
                amrex::average_down(*out[alev+1], *out[alev], 0, out[alev]->nComp(), amrrr[alev]);
#endif
            }
        }
    }

#if (AMREX_SPACEDIM != 3)
    for (int alev = 0; alev <= finest_amr_lev; ++alev) {
        linop.unapplyMetricTerm(alev, 0, *out[alev]);
    }
#endif

    for (int alev = 0; alev <= finest_amr_lev; ++alev) {
        if (cf_strategy == CFStrategy::ghostnodes)  nghost = linop.getNGrow(alev);
        out[alev]->negate(nghost);
    }
}

void
MLMG::makeSolvable ()
{
    auto const& offset = linop.getSolvabilityOffset(0, 0, rhs[0]);
    if (verbose >= 4) {
        const int ncomp = offset.size();
        for (int c = 0; c < ncomp; ++c) {
            amrex::Print() << "MLMG: Subtracting " << offset[c] << " from rhs component "
                           << c << "\n";
        }
    }
    for (int alev = 0; alev < namrlevs; ++alev) {
        linop.fixSolvabilityByOffset(alev, 0, rhs[alev], offset);
    }
}

void
MLMG::makeSolvable (int amrlev, int mglev, Any& mf)
{
    auto const& offset = linop.getSolvabilityOffset(amrlev, mglev, mf);
    if (verbose >= 4) {
        const int ncomp = offset.size();
        for (int c = 0; c < ncomp; ++c) {
            amrex::Print() << "MLMG: Subtracting " << offset[c]
                           << " from mf component c = " << c
                           << " on level (" << amrlev << ", " << mglev << ")\n";
        }
    }
    linop.fixSolvabilityByOffset(amrlev, mglev, mf, offset);
}

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
void
MLMG::bottomSolveWithHypre (Any& a_x, const Any& a_b)
{
    AMREX_ASSERT(a_x.is<MultiFab>());
    MultiFab& x = a_x.get<MultiFab>();
    MultiFab const& b = a_b.get<MultiFab>();

    const int amrlev = 0;
    const int mglev  = linop.NMGLevels(amrlev) - 1;

    const int ncomp = linop.getNComp();
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

            hypre_bndry = std::make_unique<MLMGBndry>(ba, dm, ncomp, geom);
            hypre_bndry->setHomogValues();
            const Real* dx = linop.m_geom[0][0].CellSize();
            int crse_ratio = linop.m_coarse_data_crse_ratio > 0 ? linop.m_coarse_data_crse_ratio : 1;
            RealVect bclocation(AMREX_D_DECL(0.5*dx[0]*crse_ratio,
                                             0.5*dx[1]*crse_ratio,
                                             0.5*dx[2]*crse_ratio));
            hypre_bndry->setLOBndryConds(linop.m_lobc, linop.m_hibc, -1, bclocation);
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
        makeSolvable(amrlev, mglev, a_x);
    }
}
#endif

void
MLMG::bottomSolveWithPETSc (Any& a_x, const Any& a_b)
{
#if !defined(AMREX_USE_PETSC)
    amrex::ignore_unused(a_x,a_b);
    amrex::Abort("bottomSolveWithPETSc is called without building with PETSc");
#else
    AMREX_ASSERT(a_x.is<MultiFab>());
    MultiFab& x = a_x.get<MultiFab>();
    MultiFab const& b = a_b.get<MultiFab>();

    const int ncomp = linop.getNComp();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncomp == 1, "bottomSolveWithPETSc doesn't work with ncomp > 1");

    if(petsc_solver == nullptr)
    {
        petsc_solver = linop.makePETSc();
        petsc_solver->setVerbose(bottom_verbose);

        const BoxArray& ba = linop.m_grids[0].back();
        const DistributionMapping& dm = linop.m_dmap[0].back();
        const Geometry& geom = linop.m_geom[0].back();

        petsc_bndry = std::make_unique<MLMGBndry>(ba, dm, ncomp, geom);
        petsc_bndry->setHomogValues();
        const Real* dx = linop.m_geom[0][0].CellSize();
        int crse_ratio = linop.m_coarse_data_crse_ratio > 0 ? linop.m_coarse_data_crse_ratio : 1;
        RealVect bclocation(AMREX_D_DECL(0.5*dx[0]*crse_ratio,
                                         0.5*dx[1]*crse_ratio,
                                         0.5*dx[2]*crse_ratio));
        petsc_bndry->setLOBndryConds(linop.m_lobc, linop.m_hibc, -1, bclocation);
    }
    petsc_solver->solve(x, b, bottom_reltol, Real(-1.), bottom_maxiter, *petsc_bndry, linop.getMaxOrder());
#endif
}

void
MLMG::checkPoint (const Vector<MultiFab*>& a_sol, const Vector<MultiFab const*>& a_rhs,
                  Real a_tol_rel, Real a_tol_abs, const char* a_file_name) const
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

}
