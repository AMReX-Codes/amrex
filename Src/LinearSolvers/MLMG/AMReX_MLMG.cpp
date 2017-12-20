
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_MLMG_F.H>

// sol: full solution
// rhs: rhs of the original equation L(sol) = rhs
// res: rhs of the residual equation L(cor) = res
//      usually res is the result of rhs-L(sol), but is the result of res-L(cor).
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
             Real a_tol_rel, Real a_tol_abs)
{
    BL_PROFILE_REGION("MLMG::solve()");

    Real solve_start_time = amrex::second();

    Real composite_norminf;

    prepareForSolve(a_sol, a_rhs);

    computeMLResidual(finest_amr_lev);

    bool local = true;
    Real resnorm0 = MLResNormInf(finest_amr_lev, local); 
    Real rhsnorm0 = MLRhsNormInf(local); 
    ParallelDescriptor::ReduceRealMax({resnorm0, rhsnorm0}, rhs[0].color());

    if (verbose >= 1)
    {
        amrex::Print() << "MLMG: Initial rhs               = " << rhsnorm0 << "\n"
                       << "MLMG: Initial residual (resid0) = " << resnorm0 << "\n";
    }

    Real max_norm;
    std::string norm_name;
    if (always_use_bnorm or rhsnorm0 >= resnorm0) {
        norm_name = "bnorm";
        max_norm = rhsnorm0;
    } else {
        norm_name = "resid0";
        max_norm = resnorm0;
    }
    const Real res_target = std::max(a_tol_abs, std::max(a_tol_rel,1.e-13)*max_norm);

    if (resnorm0 <= res_target)
    {
        composite_norminf = resnorm0;
        if (verbose >= 1) {
            amrex::Print() << "MLMG: No iterations needed\n";
        }
    }
    else
    {
        Real iter_start_time = amrex::second();
        bool converged = false;
        for (int iter = 0; iter < max_iters; ++iter)
        {
            oneIter(iter);

            converged = false;

            // Test convergence on the fine amr level
            computeResidual(finest_amr_lev);
            Real fine_norminf = res[finest_amr_lev][0].norm0();
            composite_norminf = fine_norminf;
            if (verbose >= 2) {
                amrex::Print() << "MLMG: Iteration " << std::setw(3) << iter+1 << " Fine resid/"
                               << norm_name << " = " << fine_norminf/max_norm << "\n";
            }
            bool fine_converged = (fine_norminf <= res_target);

            if (namrlevs == 1 and fine_converged)
            {
                converged = true;
            }
            else if (fine_converged)
            {
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
            }
            else
            {
                converged = false;
            }

            if (converged)
            {
                if (verbose >= 1) {
                    amrex::Print() << "MLMG: Final Iter. " << iter+1
                                   << " resid, resid/" << norm_name << " = "
                                   << composite_norminf << ", "
                                   << composite_norminf/max_norm << "\n";
                }
                break;
            }
        }
        if (!converged) {
            amrex::Print() << "MLMG: Failed to converge after " << max_iters << " iterations."
                           << " resid, resid/" << norm_name << " = "
                           << composite_norminf << ", "
                           << composite_norminf/max_norm << "\n";
            amrex::Abort("MLMG failed");
        }
        timer[iter_time] = amrex::second() - iter_start_time;
    }

    if (final_fill_bc) fillSolutionBC();

    timer[solve_time] = amrex::second() - solve_start_time;
    if (verbose >= 1) {
        ParallelDescriptor::ReduceRealMax(timer.data(), timer.size());
        amrex::Print() << "MLMG: Timers: Solve = " << timer[solve_time]
                       << " Iter = " << timer[iter_time]
                       << " Bottom = " << timer[bottom_time] << "\n";
    }

    return composite_norminf;
}

// in  : Residual (res) on the finest AMR level
// out : sol on all AMR levels
void
MLMG::oneIter (int iter)
{
    BL_PROFILE("MLMG::oneIter()");

    for (int alev = finest_amr_lev; alev > 0; --alev)
    {
        miniCycle(alev);

        MultiFab::Add(*sol[alev], *cor[alev][0], 0, 0, 1, 0);

        // compute residual for the coarse AMR level
        computeResWithCrseSolFineCor(alev-1,alev);

        if (alev != finest_amr_lev) {
            std::swap(cor_hold[alev][0], cor[alev][0]); // save it for the up cycle
        }
    }

    // coarest amr level
    {
        // enforce solvability if appropriate
        if (linop.isSingular(0)) {
            Real offset = res[0][0].sum() / linop.Geom(0,0).Domain().d_numPts();
            res[0][0].plus(-offset, 0, 1);
        }

        if (iter < max_fmg_iters) {
            mgFcycle ();
        } else {
            mgVcycle (0, 0);
        }

        MultiFab::Add(*sol[0], *cor[0][0], 0, 0, 1, 0);
    }

    for (int alev = 1; alev <= finest_amr_lev; ++alev)
    {
        // (Fine AMR correction) = I(Coarse AMR correction)
        interpCorrection(alev);

        MultiFab::Add(*sol[alev], *cor[alev][0], 0, 0, 1, 0);

        if (alev != finest_amr_lev) {
            MultiFab::Add(*cor_hold[alev][0], *cor[alev][0], 0, 0, 1, 0);
        }

        // Update fine AMR level correction
        computeResWithCrseCorFineCor(alev);

        miniCycle(alev);

        MultiFab::Add(*sol[alev], *cor[alev][0], 0, 0, 1, 0);

        if (alev != finest_amr_lev) {
            MultiFab::Add(*cor[alev][0], *cor_hold[alev][0], 0, 0, 1, 0);
        }
    }

    if (linop.isCellCentered())
    {
        const auto& amrrr = linop.AMRRefRatio();
        for (int falev = finest_amr_lev; falev > 0; --falev)
        {
            amrex::average_down(*sol[falev], *sol[falev-1], 0, 1, amrrr[falev-1]);
        }
    }
}

// Compute multi-level Residual (res) up to amrlevmax.
void
MLMG::computeMLResidual (int amrlevmax)
{
    BL_PROFILE("MLMG::computeMLResidual()");

    const int mglev = 0;
    for (int alev = amrlevmax; alev >= 0; --alev) {
        const MultiFab* crse_bcdata = (alev > 0) ? sol[alev-1] : nullptr;
        linop.solutionResidual(alev, res[alev][mglev], *sol[alev], rhs[alev], crse_bcdata);
        if (alev < finest_amr_lev) {
            linop.reflux(alev, res[alev][mglev], *sol[alev], rhs[alev],
                         res[alev+1][mglev], *sol[alev+1]);
        }
    }
}

// Compute single AMR level residual without masking.
void
MLMG::computeResidual (int alev)
{
    BL_PROFILE("MLMG::computeResidual()");

    MultiFab& x = *sol[alev];
    const MultiFab& b = rhs[alev];
    MultiFab& r = res[alev][0];

    const MultiFab* crse_bcdata = nullptr;
    if (alev > 0) {
        crse_bcdata = sol[alev-1];
    }
    linop.solutionResidual(alev, r, x, b, crse_bcdata);
}

// Compute coarse AMR level composite residual with coarse solution and fine correction
void
MLMG::computeResWithCrseSolFineCor (int calev, int falev)
{
    BL_PROFILE("MLMG::computeResWithCrseSolFineCor()");

    MultiFab& crse_sol = *sol[calev];
    const MultiFab& crse_rhs = rhs[calev];
    MultiFab& crse_res = res[calev][0];

    MultiFab& fine_sol = *sol[falev];
    MultiFab& fine_cor = *cor[falev][0];
    MultiFab& fine_res = res[falev][0];
    MultiFab& fine_rescor = rescor[falev][0];
    
    const MultiFab* crse_bcdata = nullptr;
    if (calev > 0) {
        crse_bcdata = sol[calev-1];
    }
    linop.solutionResidual(calev, crse_res, crse_sol, crse_rhs, crse_bcdata);

    linop.correctionResidual(falev, 0, fine_rescor, fine_cor, fine_res, BCMode::Homogeneous);
    MultiFab::Copy(fine_res, fine_rescor, 0, 0, 1, 0);

    linop.reflux(calev, crse_res, crse_sol, crse_rhs, fine_res, fine_sol);

    if (linop.isCellCentered()) {
        const int amrrr = linop.AMRRefRatio(calev);
        amrex::average_down(fine_res, crse_res, 0, 1, amrrr);
    }
}

// Compute fine AMR level residual fine_res = fine_res - L(fine_cor) with coarse providing BC.
void
MLMG::computeResWithCrseCorFineCor (int falev)
{
    BL_PROFILE("MLMG::computeResWithCrseCorFineCor()");

    const MultiFab& crse_cor = *cor[falev-1][0];

    MultiFab& fine_cor = *cor[falev][0];
    MultiFab& fine_res = res[falev][0];
    MultiFab& fine_rescor = rescor[falev][0];

    // fine_rescor = fine_res - L(fine_cor)
    linop.correctionResidual(falev, 0, fine_rescor, fine_cor, fine_res,
                             BCMode::Inhomogeneous, &crse_cor);
    MultiFab::Copy(fine_res, fine_rescor, 0, 0, 1, 0);
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
    BL_PROFILE_VAR_NS("MLMG::mgVcycle_up", blp_up);
    BL_PROFILE_VAR_NS("MLMG::mgVcycle_bottom", blp_bottom);
    BL_PROFILE_VAR_NS("MLMG::mgVcycle_down", blp_down);

    const int mglev_bottom = linop.NMGLevels(amrlev) - 1;

    BL_PROFILE_VAR_START(blp_down);
    for (int mglev = mglev_top; mglev < mglev_bottom; ++mglev)
    {
        cor[amrlev][mglev]->setVal(0.0);
        bool skip_fillboundary = true;
        for (int i = 0; i < nu1; ++i) {
            linop.smooth(amrlev, mglev, *cor[amrlev][mglev], res[amrlev][mglev],
                         skip_fillboundary);
            skip_fillboundary = false;
        }

        // rescor = res - L(cor)
        computeResOfCorrection(amrlev, mglev);

        // res_crse = R(rescor_fine); this provides res/b to the level below
        linop.restriction(amrlev, mglev+1, res[amrlev][mglev+1], rescor[amrlev][mglev]);
    }
    BL_PROFILE_VAR_STOP(blp_down);

    BL_PROFILE_VAR_START(blp_bottom);
    if (amrlev == 0)
    {
        bottomSolve();
    }
    else
    {
        cor[amrlev][mglev_bottom]->setVal(0.0);
        bool skip_fillboundary = true;
        for (int i = 0; i < nu1; ++i) {
            linop.smooth(amrlev, mglev_bottom, *cor[amrlev][mglev_bottom], res[amrlev][mglev_bottom],
                         skip_fillboundary);
            skip_fillboundary = false;
        }
    }
    BL_PROFILE_VAR_STOP(blp_bottom);

    BL_PROFILE_VAR_START(blp_up);
    for (int mglev = mglev_bottom-1; mglev >= mglev_top; --mglev)
    {
        // cor_fine += I(cor_crse)
        addInterpCorrection(amrlev, mglev);
        for (int i = 0; i < nu2; ++i) {
            linop.smooth(amrlev, mglev, *cor[amrlev][mglev], res[amrlev][mglev]);
        }
    }
    BL_PROFILE_VAR_STOP(blp_up);
}

// FMG cycle on the coarest AMR level.
// in:  Residual on the top MG level (i.e., 0)
// out: Correction (cor) on all MG levels
void
MLMG::mgFcycle ()
{
    BL_PROFILE("MLMG::mgFcycle()");

    const int amrlev = 0;
    const int ratio = 2;
    const int mg_bottom_lev = linop.NMGLevels(amrlev) - 1;

    for (int mglev = 1; mglev <= mg_bottom_lev; ++mglev)
    {
        amrex::average_down(res[amrlev][mglev-1], res[amrlev][mglev], 0, 1, ratio);
    }

    bottomSolve();

    for (int mglev = mg_bottom_lev-1; mglev >= 0; --mglev)
    {
        // cor_fine = I(cor_crse)
        interpCorrection (amrlev, mglev);

        // rescor = res - L(cor)
        computeResOfCorrection(amrlev, mglev);
        // res = rescor; this provides b to the vcycle below
        MultiFab::Copy(res[amrlev][mglev], rescor[amrlev][mglev], 0,0,1,0);

        // save cor; do v-cycle; add the saved to cor
        std::swap(cor[amrlev][mglev], cor_hold[amrlev][mglev]);
        mgVcycle(amrlev, mglev);
        MultiFab::Add(*cor[amrlev][mglev], *cor_hold[amrlev][mglev], 0, 0, 1, 0);
    }
}

// Interpolate correction from coarse to fine AMR level.
void
MLMG::interpCorrection (int alev)
{
    BL_PROFILE("MLMG::interpCorrection_1");

    const MultiFab& crse_cor = *cor[alev-1][0];
    MultiFab& fine_cor = *cor[alev][0];

    BoxArray ba = fine_cor.boxArray();
    const int amrrr = linop.AMRRefRatio(alev-1);
    IntVect refratio{amrrr};
    ba.coarsen(refratio);

    const Geometry& crse_geom = linop.Geom(alev-1,0);

    const int ng = linop.isCellCentered() ? 1 : 0;
    MultiFab cfine(ba, fine_cor.DistributionMap(), 1, ng);
    cfine.setVal(0.0);
    cfine.ParallelCopy(crse_cor, 0, 0, 1, 0, ng, crse_geom.periodicity());

    if (linop.isCellCentered())
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(fine_cor, MFItInfo().EnableTiling().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            amrex_mlmg_lin_cc_interp(BL_TO_FORTRAN_BOX(bx),
                                     BL_TO_FORTRAN_ANYD(fine_cor[mfi]),
                                     BL_TO_FORTRAN_ANYD(cfine[mfi]),
                                     &refratio[0]);
        }
    }
    else
    {
        AMREX_ALWAYS_ASSERT(amrrr == 2);
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            FArrayBox tmpfab;
            for (MFIter mfi(fine_cor, true); mfi.isValid(); ++mfi)
            {
                const Box& fbx = mfi.tilebox();
                const Box& cbx = amrex::coarsen(fbx,2);
                const Box& tmpbx = amrex::refine(cbx,2);
                tmpfab.resize(tmpbx);
                amrex_mlmg_lin_nd_interp(BL_TO_FORTRAN_BOX(cbx),
                                         BL_TO_FORTRAN_BOX(tmpbx),
                                         BL_TO_FORTRAN_ANYD(tmpfab),
                                         BL_TO_FORTRAN_ANYD(cfine[mfi]));
                fine_cor[mfi].copy(tmpfab, fbx, 0, fbx, 0, 1);
            }
        }
    }
}

// Interpolate correction between MG levels
// inout: Correction (cor) on coarse MG lev.  (out due to FillBoundary)
// out  : Correction (cor) on fine MG lev.
void
MLMG::interpCorrection (int alev, int mglev)
{
    BL_PROFILE("MLMG::interpCorrection_2");

    MultiFab& crse_cor = *cor[alev][mglev+1];
    MultiFab& fine_cor = *cor[alev][mglev  ];

    const Geometry& crse_geom = linop.Geom(alev,mglev+1);
    const int refratio = 2;

    MultiFab cfine;
    const MultiFab* cmf;
    
    if (amrex::isMFIterSafe(crse_cor, fine_cor))
    {
        crse_cor.FillBoundary(crse_geom.periodicity());
        cmf = &crse_cor;
    }
    else
    {
        BoxArray cba = fine_cor.boxArray();
        cba.coarsen(refratio);
        const int nc = crse_cor.nComp();
        const int ng = crse_cor.nGrow();
        cfine.define(cba, fine_cor.DistributionMap(), nc, ng);
        cfine.setVal(0.0);
        cfine.ParallelCopy(crse_cor, 0, 0, nc, 0, ng, crse_geom.periodicity());
        cmf = & cfine;
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(fine_cor, MFItInfo().EnableTiling().SetDynamic(true)); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        amrex_mlmg_lin_cc_interp(BL_TO_FORTRAN_BOX(bx),
                                 BL_TO_FORTRAN_ANYD(fine_cor[mfi]),
                                 BL_TO_FORTRAN_ANYD(  (*cmf)[mfi]),
                                 &refratio);
    }
}

// (Fine MG level correction) += I(Coarse MG level correction)
void
MLMG::addInterpCorrection (int alev, int mglev)
{
    BL_PROFILE("MLMG::addInterpCorrection()");

    const MultiFab& crse_cor = *cor[alev][mglev+1];
    MultiFab&       fine_cor = *cor[alev][mglev  ];

    const int refratio = 2;
    MultiFab cfine;
    const MultiFab* cmf;

    if (amrex::isMFIterSafe(crse_cor, fine_cor))
    {
        cmf = &crse_cor;
    }
    else
    {
        BoxArray cba = fine_cor.boxArray();
        cba.coarsen(refratio);
        const int nc = crse_cor.nComp();
        const int ng = 0;
        cfine.define(cba, fine_cor.DistributionMap(), nc, ng);
        cfine.ParallelCopy(crse_cor);
        cmf = &cfine;
    }

    linop.interpolation(alev, mglev, fine_cor, *cmf);
}

// Compute rescor = res - L(cor)
// in   : res
// inout: cor (out due to FillBoundary in linop.correctionResidual)
// out  : rescor
void
MLMG::computeResOfCorrection (int amrlev, int mglev)
{
    BL_PROFILE("MLMG:computeResOfCorrection()");
    MultiFab& x = *cor[amrlev][mglev];
    const MultiFab& b = res[amrlev][mglev];
    MultiFab& r = rescor[amrlev][mglev];
    linop.correctionResidual(amrlev, mglev, r, x, b, BCMode::Homogeneous);
}

// At the true bottom of the coarset AMR level.
// in  : Residual (res) as b
// out : Correction (cor) as x
void
MLMG::bottomSolve ()
{
    BL_PROFILE("MLMG::bottomSolve()");

    if (!linop.isBottomActive()) return;

    Real bottom_start_time = amrex::second();

    int old_sn = ParallelDescriptor::SeqNum(3);
    
    const int amrlev = 0;
    const int mglev = linop.NMGLevels(amrlev) - 1;
    MultiFab& x = *cor[amrlev][mglev];
    MultiFab& b = res[amrlev][mglev];

    if (linop.isSingular(amrlev)) {
        const bool local = true;
        Real offset = b.sum(0,local) / linop.Geom(amrlev,mglev).Domain().d_numPts();
        ParallelAllReduce::Sum(offset, linop.BottomCommunicator());
        b.plus(-offset, 0, 1);
    }

    x.setVal(0.0);

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
        MLCGSolver::Solver solver_type;
        if (bottom_solver == BottomSolver::bicgstab)
        {
            solver_type = MLCGSolver::Solver::BiCGStab;
        }
        else if (bottom_solver == BottomSolver::cg)
        {
            solver_type = MLCGSolver::Solver::CG;         
        }

        MLCGSolver cg_solver(linop, solver_type);
        cg_solver.setVerbose(cg_verbose);
        cg_solver.setMaxIter(cg_maxiter);
        cg_solver.setUnstableCriterion(cg_unstable_criterion);

        const Real cg_rtol = 1.e-4;
        const Real cg_atol = -1.0;
        int ret = cg_solver.solve(x, b, cg_rtol, cg_atol);
        if (ret != 0 && verbose >= 1) {
            amrex::Print() << "MLMG: Bottom solve failed.\n";
        }
        const int n = ret==0 ? nub : nuf;
        for (int i = 0; i < n; ++i) {
            linop.smooth(amrlev, mglev, x, b);
        }
    }

    ParallelDescriptor::SeqNum(2, old_sn);

    timer[bottom_time] += amrex::second() - bottom_start_time;
}

// Compute single-level masked inf-norm of Residual (res).
Real
MLMG::ResNormInf (int alev, bool local)
{
    BL_PROFILE("MLMG::ResNormInf()");
    const int mglev = 0;
    if (alev < finest_amr_lev) {
        return res[alev][mglev].norm0(fine_mask[alev],0,0,local);
    } else {
        return res[alev][mglev].norm0(0,0,local);
    }
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
    if (!local) ParallelDescriptor::ReduceRealMax(r, rhs[0].color());
    return r;
}

// Compute multi-level masked inf-norm of RHS (rhs).
Real
MLMG::MLRhsNormInf (bool local)
{
    BL_PROFILE("MLMG::MLRhsNormInf()");
    Real r = 0.0;
    for (int alev = 0; alev <= finest_amr_lev; ++alev)
    {
        if (alev < finest_amr_lev) {
            r = std::max(r, rhs[alev].norm0(fine_mask[alev],0,0,local));
        } else {
            r = std::max(r, rhs[alev].norm0(0,0,local));
        }
    }
    return r;
}

void
MLMG::buildFineMask ()
{
    BL_PROFILE("MLMG::buildFineMask()");

    fine_mask.clear();
    fine_mask.resize(namrlevs-1);
    
    const auto& amrrr = linop.AMRRefRatio();
    for (int alev = 0; alev < finest_amr_lev; ++alev)
    {
        fine_mask[alev].define(rhs[alev].boxArray(), rhs[alev].DistributionMap(), 1, 0);
        fine_mask[alev].setVal(1);

        BoxArray baf = rhs[alev+1].boxArray();
        baf.coarsen(amrrr[alev]);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(fine_mask[alev], MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            auto& fab = fine_mask[alev][mfi];

            const std::vector< std::pair<int,Box> >& isects = baf.intersections(fab.box());

            for (int ii = 0; ii < isects.size(); ++ii)
            {
                fab.setVal(0,isects[ii].second,0);
            }
        }
    }
}

void
MLMG::prepareForSolve (const Vector<MultiFab*>& a_sol, const Vector<MultiFab const*>& a_rhs)
{
    BL_PROFILE("MLMG::prepareForSolve()");

    AMREX_ASSERT(a_sol[0]->nGrow() > 0);
    AMREX_ASSERT(namrlevs <= a_sol.size());
    AMREX_ASSERT(namrlevs <= a_rhs.size());

    timer.assign(ntimers, 0.0);

    linop.prepareForSolve();
    
    sol = a_sol;
    
    rhs.resize(namrlevs);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        rhs[alev].define(a_rhs[alev]->boxArray(), a_rhs[alev]->DistributionMap(), 1, 0);
        MultiFab::Copy(rhs[alev], *a_rhs[alev], 0, 0, 1, 0);
        linop.applyMetricTerm(alev, 0, rhs[alev]);
    }

    if (linop.isCellCentered())
    {
        const auto& amrrr = linop.AMRRefRatio();
        for (int falev = finest_amr_lev; falev > 0; --falev)
        {
            amrex::average_down(*sol[falev], *sol[falev-1], 0, 1, amrrr[falev-1]);
            amrex::average_down( rhs[falev],  rhs[falev-1], 0, 1, amrrr[falev-1]);
        }

        // enforce solvability if appropriate
        if (linop.isSingular(0))
        {
            Real offset = rhs[0].sum() / linop.Geom(0,0).Domain().d_numPts();
            if (verbose >= 4) {
                amrex::Print() << "MLMG: Subtracting " << offset << " from rhs\n";
            }
            for (int alev = 0; alev < namrlevs; ++alev) {
                rhs[alev].plus(-offset, 0, 1);
            }
        }
    }
    else
    {
        const auto& amrrr = linop.AMRRefRatio();
        for (int falev = finest_amr_lev; falev > 0; --falev)
        {
            AMREX_ALWAYS_ASSERT(amrrr[0] == 2);
            amrex::average_down(rhs[falev],  rhs[falev-1], 0, 1, 2);
        }        
    }

    const int nc = 1;
    int ng = linop.isCellCentered() ? 0 : 1;
    linop.make(res, nc, ng);
    linop.make(rescor, nc, ng);

    ng = 1;
    cor.resize(namrlevs);
    for (int alev = 0; alev <= finest_amr_lev; ++alev)
    {
        const int nmglevs = linop.NMGLevels(alev);
        cor[alev].resize(nmglevs);
        for (int mglev = 0; mglev < nmglevs; ++mglev)
        {
            cor[alev][mglev].reset(new MultiFab(res[alev][mglev].boxArray(),
                                                res[alev][mglev].DistributionMap(),
                                                nc, ng));
            cor[alev][mglev]->setVal(0.0);
        }
    }

    cor_hold.resize(std::max(namrlevs-1,1));
    {
        const int alev = 0;
        const int nmglevs = linop.NMGLevels(alev);
        cor_hold[alev].resize(nmglevs);
        for (int mglev = 0; mglev < nmglevs-1; ++mglev)
        {
            cor_hold[alev][mglev].reset(new MultiFab(cor[alev][mglev]->boxArray(),
                                                     cor[alev][mglev]->DistributionMap(),
                                                     nc, ng));
            cor_hold[alev][mglev]->setVal(0.0);
        }
    }
    for (int alev = 1; alev < finest_amr_lev; ++alev)
    {
        cor_hold[alev].resize(1);
        cor_hold[alev][0].reset(new MultiFab(cor[alev][0]->boxArray(),
                                             cor[alev][0]->DistributionMap(),
                                             nc, ng));
        cor_hold[alev][0]->setVal(0.0);
    }

    buildFineMask();
}

void
MLMG::getGradSolution (const Vector<std::array<MultiFab*,AMREX_SPACEDIM> >& a_grad_sol)
{
    BL_PROFILE("MLMG::getGradSolution()");
    for (int alev = 0; alev <= finest_amr_lev; ++alev) {
        linop.compGrad(alev, a_grad_sol[alev], *sol[alev]);
    }
}

void
MLMG::getFluxes (const Vector<std::array<MultiFab*,AMREX_SPACEDIM> >& a_flux)
{
    BL_PROFILE("MLMG::getFluxes()");
    const Real betainv = 1.0 / linop.getBScalar();
    for (int alev = 0; alev <= finest_amr_lev; ++alev) {
        linop.compFlux(alev, a_flux[alev], *sol[alev]);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            linop.unapplyMetricTerm(alev, 0, *a_flux[alev][idim]);
            if (betainv != 1.0) {
                a_flux[alev][idim]->mult(betainv);
            }
        }
    }
}

void
MLMG::compResidual (const Vector<MultiFab*>& a_res, const Vector<MultiFab*>& a_sol,
                    const Vector<MultiFab const*>& a_rhs)
{
    BL_PROFILE("MLMG::compResidual()");

    linop.prepareForSolve();
    
    const auto& amrrr = linop.AMRRefRatio();

    for (int alev = finest_amr_lev; alev >= 0; --alev) {
        const MultiFab* crse_bcdata = (alev > 0) ? a_sol[alev-1] : nullptr;
        const MultiFab* prhs = a_rhs[alev];
#if (AMREX_SPACEDIM != 3)
        MultiFab rhstmp(prhs->boxArray(), prhs->DistributionMap(), 1, 0);
        MultiFab::Copy(rhstmp, *prhs, 0, 0, 1, 0);
        linop.applyMetricTerm(alev, 0, rhstmp);
        prhs = &rhstmp;
#endif
        linop.solutionResidual(alev, *a_res[alev], *a_sol[alev], *prhs, crse_bcdata);
        if (alev < finest_amr_lev) {
            linop.reflux(alev, *a_res[alev], *a_sol[alev], *prhs,
                         *a_res[alev+1], *a_sol[alev+1]);
            if (linop.isCellCentered()) {
                amrex::average_down(*a_res[alev+1], *a_res[alev], 0, 1, amrrr[alev]);
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
MLMG::fillSolutionBC ()
{
    BL_PROFILE("MLMG::fillSolutionBC()");

    for (int alev = 0; alev <= finest_amr_lev; ++alev)
    {
        const MultiFab* crse_bcdata = (alev > 0) ? sol[alev-1] : nullptr;
        linop.fillSolutionBC(alev, *sol[alev], crse_bcdata);
    }
}

}
