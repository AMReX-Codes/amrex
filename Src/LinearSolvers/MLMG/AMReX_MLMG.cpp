
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_MLMG_F.H>
#include <AMReX_MLABecLaplacian.H>

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

    bool is_nsolve = linop.m_parent;

    Real solve_start_time = amrex::second();

    Real composite_norminf;

    prepareForSolve(a_sol, a_rhs);

    computeMLResidual(finest_amr_lev);

    int ncomp = linop.getNComp();

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

    if (!is_nsolve && resnorm0 <= res_target)
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

        const int niters = do_fixed_number_of_iters ? do_fixed_number_of_iters : max_iters;
        for (int iter = 0; iter < niters; ++iter)
        {
            oneIter(iter);

            converged = false;

            // Test convergence on the fine amr level
            computeResidual(finest_amr_lev);

            if (is_nsolve) continue;

            Real fine_norminf = ResNormInf(finest_amr_lev);
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
        if (!converged && do_fixed_number_of_iters == 0) {
            amrex::Print() << "MLMG: Failed to converge after " << max_iters << " iterations."
                           << " resid, resid/" << norm_name << " = "
                           << composite_norminf << ", "
                           << composite_norminf/max_norm << "\n";
            amrex::Abort("MLMG failed");
        }
        timer[iter_time] = amrex::second() - iter_start_time;
    }

    int ng_back = final_fill_bc ? 1 : 0;
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (a_sol[alev] != sol[alev])
        {
            MultiFab::Copy(*a_sol[alev], *sol[alev], 0, 0, ncomp, ng_back);
        }
    }

    timer[solve_time] = amrex::second() - solve_start_time;
    if (verbose >= 1) {
        ParallelReduce::Max<Real>(timer.data(), timer.size(), 0,
                                  ParallelContext::CommunicatorSub());
        if (ParallelContext::MyProcSub() == 0)
        {
            amrex::AllPrint() << "MLMG: Timers: Solve = " << timer[solve_time]
                              << " Iter = " << timer[iter_time]
                              << " Bottom = " << timer[bottom_time] << "\n";
        }
    }

    return composite_norminf;
}

// in  : Residual (res) on the finest AMR level
// out : sol on all AMR levels
void
MLMG::oneIter (int iter)
{
    BL_PROFILE("MLMG::oneIter()");

    int ncomp = linop.getNComp();

    for (int alev = finest_amr_lev; alev > 0; --alev)
    {
        miniCycle(alev);

        MultiFab::Add(*sol[alev], *cor[alev][0], 0, 0, ncomp, 0);

        // compute residual for the coarse AMR level
        computeResWithCrseSolFineCor(alev-1,alev);

        if (alev != finest_amr_lev) {
            std::swap(cor_hold[alev][0], cor[alev][0]); // save it for the up cycle
        }
    }

    // coarest amr level
    {
        // enforce solvability if appropriate
        if (linop.isSingular(0))
        {
            if (linop.isCellCentered())
            {
                Real npinv = 1.0 / linop.Geom(0,0).Domain().d_numPts();
                Vector<Real> offset(ncomp);
                for (int c = 0; c < ncomp; ++c) {
                    offset[c] = res[0][0].sum(c, true) * npinv;
                }
                ParallelAllReduce::Sum(offset.data(), ncomp, ParallelContext::CommunicatorSub());
                for (int c = 0; c < ncomp; ++c) {
                    res[0][0].plus(-offset[c], c, 1);
                }
            }
            else
            {
                Real offset = getNodalSum(0, 0, res[0][0]);
                res[0][0].plus(-offset, 0, 1);
            }
        }

        if (iter < max_fmg_iters) {
            mgFcycle ();
        } else {
            mgVcycle (0, 0);
        }

        MultiFab::Add(*sol[0], *cor[0][0], 0, 0, ncomp, 0);
    }

    for (int alev = 1; alev <= finest_amr_lev; ++alev)
    {
        // (Fine AMR correction) = I(Coarse AMR correction)
        interpCorrection(alev);

        MultiFab::Add(*sol[alev], *cor[alev][0], 0, 0, ncomp, 0);

        if (alev != finest_amr_lev) {
	  MultiFab::Add(*cor_hold[alev][0], *cor[alev][0], 0, 0, ncomp, 0);
        }

        // Update fine AMR level correction
        computeResWithCrseCorFineCor(alev);

        miniCycle(alev);

        MultiFab::Add(*sol[alev], *cor[alev][0], 0, 0, ncomp, 0);

        if (alev != finest_amr_lev) {
	  MultiFab::Add(*cor[alev][0], *cor_hold[alev][0], 0, 0, ncomp, 0);
        }
    }

    averageDownAndSync();
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
                         res[alev+1][mglev], *sol[alev+1], rhs[alev+1]);
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

    int ncomp = linop.getNComp();

    MultiFab& crse_sol = *sol[calev];
    const MultiFab& crse_rhs = rhs[calev];
    MultiFab& crse_res = res[calev][0];

    MultiFab& fine_sol = *sol[falev];
    const MultiFab& fine_rhs = rhs[falev];
    MultiFab& fine_cor = *cor[falev][0];
    MultiFab& fine_res = res[falev][0];
    MultiFab& fine_rescor = rescor[falev][0];
    
    const MultiFab* crse_bcdata = nullptr;
    if (calev > 0) {
        crse_bcdata = sol[calev-1];
    }
    linop.solutionResidual(calev, crse_res, crse_sol, crse_rhs, crse_bcdata);

    linop.correctionResidual(falev, 0, fine_rescor, fine_cor, fine_res, BCMode::Homogeneous);
    MultiFab::Copy(fine_res, fine_rescor, 0, 0, ncomp, 0);

    linop.reflux(calev, crse_res, crse_sol, crse_rhs, fine_res, fine_sol, fine_rhs);

    if (linop.isCellCentered()) {
        const int amrrr = linop.AMRRefRatio(calev);
        amrex::average_down(fine_res, crse_res, 0, ncomp, amrrr);
    }
}

// Compute fine AMR level residual fine_res = fine_res - L(fine_cor) with coarse providing BC.
void
MLMG::computeResWithCrseCorFineCor (int falev)
{
    BL_PROFILE("MLMG::computeResWithCrseCorFineCor()");

    int ncomp = linop.getNComp();

    const MultiFab& crse_cor = *cor[falev-1][0];

    MultiFab& fine_cor = *cor[falev][0];
    MultiFab& fine_res = res[falev][0];
    MultiFab& fine_rescor = rescor[falev][0];

    // fine_rescor = fine_res - L(fine_cor)
    linop.correctionResidual(falev, 0, fine_rescor, fine_cor, fine_res,
                             BCMode::Inhomogeneous, &crse_cor);
    MultiFab::Copy(fine_res, fine_rescor, 0, 0, ncomp, 0);
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
        if (verbose >= 4)
        {
            Real norm = res[amrlev][mglev].norm0();
            amrex::Print() << "AT LEVEL "                << mglev << "\n"
                           << "   DN: Norm before smooth " << norm << "\n";
        }

        cor[amrlev][mglev]->setVal(0.0);
        bool skip_fillboundary = true;
        for (int i = 0; i < nu1; ++i) {
            linop.smooth(amrlev, mglev, *cor[amrlev][mglev], res[amrlev][mglev],
                         skip_fillboundary);
            skip_fillboundary = false;
        }

        // rescor = res - L(cor)
        computeResOfCorrection(amrlev, mglev);

        if (verbose >= 4)
        {
            Real norm = rescor[amrlev][mglev].norm0();
            amrex::Print() << "   DN: Norm after  smooth " << norm << "\n";
        }

        // res_crse = R(rescor_fine); this provides res/b to the level below
        linop.restriction(amrlev, mglev+1, res[amrlev][mglev+1], rescor[amrlev][mglev]);
    }
    BL_PROFILE_VAR_STOP(blp_down);

    BL_PROFILE_VAR_START(blp_bottom);
    if (amrlev == 0)
    {
        if (verbose >= 4)
        {
            Real norm = res[amrlev][mglev_bottom].norm0();
            amrex::Print() << "AT LEVEL "                << mglev_bottom << "\n"
                           << "   DN: Norm before bottom " << norm << "\n";
        }
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
        if (verbose >= 4)
        {
            computeResOfCorrection(amrlev, mglev);
            Real norm = rescor[amrlev][mglev].norm0();
            amrex::Print() << "AT LEVEL "                << mglev << "\n"
                           << "   UP: Norm before smooth " << norm << "\n";
        }
        for (int i = 0; i < nu2; ++i) {
            linop.smooth(amrlev, mglev, *cor[amrlev][mglev], res[amrlev][mglev]);
        }
        if (verbose >= 4)
        {
            computeResOfCorrection(amrlev, mglev);
            Real norm = rescor[amrlev][mglev].norm0();
            amrex::Print() << "AT LEVEL "                << mglev << "\n"
                           << "   UP: Norm after  smooth " << norm << "\n";
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
    const int ncomp = linop.getNComp();

    for (int mglev = 1; mglev <= mg_bottom_lev; ++mglev)
    {
        amrex::average_down(res[amrlev][mglev-1], res[amrlev][mglev], 0, ncomp, ratio);
    }

    bottomSolve();

    for (int mglev = mg_bottom_lev-1; mglev >= 0; --mglev)
    {
        // cor_fine = I(cor_crse)
        interpCorrection (amrlev, mglev);

        // rescor = res - L(cor)
        computeResOfCorrection(amrlev, mglev);
        // res = rescor; this provides b to the vcycle below
        MultiFab::Copy(res[amrlev][mglev], rescor[amrlev][mglev], 0,0,ncomp,0);

        // save cor; do v-cycle; add the saved to cor
        std::swap(cor[amrlev][mglev], cor_hold[amrlev][mglev]);
        mgVcycle(amrlev, mglev);
        MultiFab::Add(*cor[amrlev][mglev], *cor_hold[amrlev][mglev], 0, 0, ncomp, 0);
    }
}

// Interpolate correction from coarse to fine AMR level.
void
MLMG::interpCorrection (int alev)
{
    BL_PROFILE("MLMG::interpCorrection_1");

    const int ncomp = linop.getNComp();

    const MultiFab& crse_cor = *cor[alev-1][0];
    MultiFab& fine_cor = *cor[alev][0];

    BoxArray ba = fine_cor.boxArray();
    const int amrrr = linop.AMRRefRatio(alev-1);
    IntVect refratio{amrrr};
    ba.coarsen(refratio);

    const Geometry& crse_geom = linop.Geom(alev-1,0);

    const int ng = linop.isCellCentered() ? 1 : 0;
    MultiFab cfine(ba, fine_cor.DistributionMap(), ncomp, ng);
    cfine.setVal(0.0);
    cfine.ParallelCopy(crse_cor, 0, 0, ncomp, 0, ng, crse_geom.periodicity());

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
                                     &refratio[0],
				     &ncomp);
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
                                         BL_TO_FORTRAN_ANYD(cfine[mfi]),
					 &ncomp);
                fine_cor[mfi].copy(tmpfab, fbx, 0, fbx, 0, ncomp);
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

    const int ncomp = linop.getNComp();

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
        const int ng = linop.isCellCentered() ? crse_cor.nGrow() : 0;
        cfine.define(cba, fine_cor.DistributionMap(), ncomp, ng);
        cfine.setVal(0.0);
        cfine.ParallelCopy(crse_cor, 0, 0, ncomp, 0, ng, crse_geom.periodicity());
        cmf = & cfine;
    }

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
                                     BL_TO_FORTRAN_ANYD(  (*cmf)[mfi]),
                                     &refratio,&ncomp);
        }
    }
    else
    {
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
                                         BL_TO_FORTRAN_ANYD((*cmf)[mfi]),
					 &ncomp);
                fine_cor[mfi].copy(tmpfab, fbx, 0, fbx, 0, ncomp);
            }
        }
    }
}

// (Fine MG level correction) += I(Coarse MG level correction)
void
MLMG::addInterpCorrection (int alev, int mglev)
{
    BL_PROFILE("MLMG::addInterpCorrection()");

    const int ncomp = linop.getNComp();

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
        const int ng = 0;
        cfine.define(cba, fine_cor.DistributionMap(), ncomp, ng);
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

    a_rhs.setVal(0.0);
    a_rhs.ParallelCopy(res[0].back());

    a_solver.solve({&a_sol}, {&a_rhs}, -1.0, -1.0);

    cor[0].back()->ParallelCopy(a_sol);
}

void
MLMG::actualBottomSolve ()
{
    BL_PROFILE("MLMG::actualBottomSolve()");

    const int ncomp = linop.getNComp();

    if (!linop.isBottomActive()) return;

    Real bottom_start_time = amrex::second();

    ParallelContext::push(linop.BottomCommunicator());

    const int amrlev = 0;
    const int mglev = linop.NMGLevels(amrlev) - 1;
    MultiFab& x = *cor[amrlev][mglev];
    MultiFab& b = res[amrlev][mglev];

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
        MultiFab* bottom_b = &b;
        MultiFab raii_b;
        if (linop.isBottomSingular())
        {
            raii_b.define(b.boxArray(), b.DistributionMap(), ncomp, b.nGrow());
            MultiFab::Copy(raii_b,b,0,0,ncomp,b.nGrow());
            bottom_b = &raii_b;

            Vector<Real> offset(ncomp);
            if (linop.isCellCentered())
            {
                Real npinv = 1.0 / linop.Geom(amrlev,mglev).Domain().d_numPts();
                for (int c = 0; c < ncomp; ++c) {
                    offset[c] = bottom_b->sum(c,true) * npinv;
                }
                ParallelAllReduce::Sum(offset.data(), ncomp, linop.BottomCommunicator());
            }
            else
            {
                AMREX_ASSERT_WITH_MESSAGE(ncomp==1, "ncomp > 1 not supported for singular nodal problem");
                offset[0] = getNodalSum(amrlev, mglev, *bottom_b);
            }

            for (int c = 0; c < ncomp; ++c) {
                bottom_b->plus(-offset[c], c, 1);
            }
        }

        if (bottom_solver == BottomSolver::hypre)
        {
            bottomSolveWithHypre(x, *bottom_b);
        }
        else
        {
            MLCGSolver cg_solver(linop);
            cg_solver.setVerbose(bottom_verbose);
            cg_solver.setMaxIter(bottom_maxiter);
            
            const Real cg_rtol = 1.e-4;
            const Real cg_atol = -1.0;
            int ret = cg_solver.solve(x, *bottom_b, cg_rtol, cg_atol);
            if (ret != 0 && verbose >= 1) {
                amrex::Print() << "MLMG: Bottom solve failed.\n";
            }
            const int n = ret==0 ? nub : nuf;
            for (int i = 0; i < n; ++i) {
                linop.smooth(amrlev, mglev, x, b);
            }
        }
    }

    ParallelContext::pop();

    timer[bottom_time] += amrex::second() - bottom_start_time;
}

// Compute single-level masked inf-norm of Residual (res).
Real
MLMG::ResNormInf (int alev, bool local)
{
    BL_PROFILE("MLMG::ResNormInf()");
    const int mglev = 0;
    Real norm = 0.0;
    for (int n = 0; n < linop.getNComp(); n++)
      {
	Real newnorm = 0.0;
	if (fine_mask[alev]) {
	  newnorm = res[alev][mglev].norm0(*fine_mask[alev],n,0,local);
	} else {
	  newnorm = res[alev][mglev].norm0(n,0,local);
	}
	if (newnorm > norm) norm = newnorm;
      }
    return norm;
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
    Real r = 0.0;
    for (int alev = 0; alev <= finest_amr_lev; ++alev)
    {
        if (alev < finest_amr_lev) {
            r = std::max(r, rhs[alev].norm0(*fine_mask[alev],0,0,local));
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

    const int ncomp = linop.getNComp();

    fine_mask.clear();
    fine_mask.resize(namrlevs);
    
    const auto& amrrr = linop.AMRRefRatio();
    for (int alev = 0; alev < finest_amr_lev; ++alev)
    {
        fine_mask[alev].reset(new iMultiFab(rhs[alev].boxArray(), rhs[alev].DistributionMap(), ncomp, 0));
        fine_mask[alev]->setVal(1);

        BoxArray baf = rhs[alev+1].boxArray();
        baf.coarsen(amrrr[alev]);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*fine_mask[alev], MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            auto& fab = (*fine_mask[alev])[mfi];

            const std::vector< std::pair<int,Box> >& isects = baf.intersections(fab.box());

            for (int ii = 0; ii < isects.size(); ++ii)
            {
                fab.setVal(0,isects[ii].second,0);
            }
        }
    }

    if (!linop.isCellCentered()) {
        for (int alev = 0; alev < finest_amr_lev; ++alev) {
            linop.fixUpResidualMask(alev, *fine_mask[alev]);
        }
    }
}

void
MLMG::prepareForSolve (const Vector<MultiFab*>& a_sol, const Vector<MultiFab const*>& a_rhs)
{
    BL_PROFILE("MLMG::prepareForSolve()");

    AMREX_ASSERT(namrlevs <= a_sol.size());
    AMREX_ASSERT(namrlevs <= a_rhs.size());

    timer.assign(ntimers, 0.0);

    const int ncomp = linop.getNComp();

    if (!linop_prepared) {
        linop.prepareForSolve();
        linop_prepared = true;
    }

    sol.resize(namrlevs);
    sol_raii.resize(namrlevs);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (a_sol[alev]->nGrow() == 1)
        {
            sol[alev] = a_sol[alev];
        }
        else
        {
            sol_raii[alev].reset(new MultiFab(a_sol[alev]->boxArray(),
                                              a_sol[alev]->DistributionMap(), ncomp, 1));
            sol_raii[alev]->setVal(0.0);
            MultiFab::Copy(*sol_raii[alev], *a_sol[alev], 0, 0, ncomp, 0);
            sol[alev] = sol_raii[alev].get();
        }
    }
    
    rhs.resize(namrlevs);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
      rhs[alev].define(a_rhs[alev]->boxArray(), a_rhs[alev]->DistributionMap(), ncomp, 0);
      MultiFab::Copy(rhs[alev], *a_rhs[alev], 0, 0, ncomp, 0);
      linop.applyMetricTerm(alev, 0, rhs[alev]);
    }

    for (int falev = finest_amr_lev; falev > 0; --falev)
    {
        linop.averageDownSolutionRHS(falev-1, *sol[falev-1], rhs[falev-1], *sol[falev], rhs[falev]);
    }
    
    // enforce solvability if appropriate
    if (linop.isSingular(0))
    {
        if (linop.isCellCentered())
        {
            Real npinv = 1.0 / linop.Geom(0,0).Domain().d_numPts();
            Vector<Real> offset(ncomp);
            for (int c = 0; c < ncomp; ++c) {
                offset[c] = rhs[0].sum(c,true) * npinv;
            }
            ParallelAllReduce::Sum(offset.data(), ncomp, ParallelContext::CommunicatorSub());
            if (verbose >= 4) {
                for (int c = 0; c < ncomp; ++c) {
                    amrex::Print() << "MLMG: Subtracting " << offset[c] 
                                   << " from rhs component " << c << "\n";
                }
            }
            for (int alev = 0; alev < namrlevs; ++alev) {
                for (int c = 0; c < ncomp; ++c) {
                    rhs[alev].plus(-offset[c], c, 1);
                }
            }
        }
        else
        {
            Real offset = getNodalSum(0, 0, rhs[0]);
            for (int alev = 0; alev < namrlevs; ++alev) {
                rhs[alev].plus(-offset, 0, 1);
            }
        }
    }

    int ng = linop.isCellCentered() ? 0 : 1;
    linop.make(res, ncomp, ng);
    linop.make(rescor, ncomp, ng);
    for (int alev = 0; alev <= finest_amr_lev; ++alev)
    {
        const int nmglevs = linop.NMGLevels(alev);
        for (int mglev = 0; mglev < nmglevs; ++mglev)
        {
            rescor[alev][mglev].setVal(0.0);
        }
    }

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
                                                ncomp, ng));
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
                                                     ncomp, ng));
            cor_hold[alev][mglev]->setVal(0.0);
        }
    }
    for (int alev = 1; alev < finest_amr_lev; ++alev)
    {
        cor_hold[alev].resize(1);
        cor_hold[alev][0].reset(new MultiFab(cor[alev][0]->boxArray(),
                                             cor[alev][0]->DistributionMap(),
                                             ncomp, ng));
        cor_hold[alev][0]->setVal(0.0);
    }

    buildFineMask();

    if (linop.m_parent) do_nsolve = false;  // no embeded N-Solve
    if (linop.m_domain_covered[0]) do_nsolve = false;
    if (linop.doAgglomeration()) do_nsolve = false;
    if (AMREX_SPACEDIM != 3) do_nsolve = false;

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
    ns_linop = std::move(linop.makeNLinOp(nsolve_grid_size));

    const int ncomp = linop.getNComp();

    const BoxArray& ba = (*ns_linop).m_grids[0][0];
    const DistributionMapping& dm =(*ns_linop).m_dmap[0][0]; 

    ns_sol.reset(new MultiFab(ba, dm, ncomp, 1));
    ns_rhs.reset(new MultiFab(ba, dm, ncomp, 0));
    ns_sol->setVal(0.0);
    ns_rhs->setVal(0.0);

    ns_linop->setLevelBC(0, ns_sol.get());
    
    ns_mlmg.reset(new MLMG(*ns_linop));
    ns_mlmg->setVerbose(0);
    ns_mlmg->setFixedIter(1);
    ns_mlmg->setMaxFmgIter(20);
    ns_mlmg->setBottomSolver(BottomSolver::smoother);
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

    const int ncomp = linop.getNComp();
   
    sol.resize(namrlevs);
    sol_raii.resize(namrlevs);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (a_sol[alev]->nGrow() == 1)
        {
            sol[alev] = a_sol[alev];
        }
        else
        {
            if (sol_raii[alev] == nullptr)
            {
                sol_raii[alev].reset(new MultiFab(a_sol[alev]->boxArray(),
                                                  a_sol[alev]->DistributionMap(), ncomp, 1));
            }
            MultiFab::Copy(*sol_raii[alev], *a_sol[alev], 0, 0, ncomp, 0);
            sol[alev] = sol_raii[alev].get();
        }
    }

    if (!linop_prepared) {
        linop.prepareForSolve();
        linop_prepared = true;
    }
    
    const auto& amrrr = linop.AMRRefRatio();

    for (int alev = finest_amr_lev; alev >= 0; --alev) {
        const MultiFab* crse_bcdata = (alev > 0) ? sol[alev-1] : nullptr;
        const MultiFab* prhs = a_rhs[alev];
#if (AMREX_SPACEDIM != 3)
        MultiFab rhstmp(prhs->boxArray(), prhs->DistributionMap(), ncomp, 0);
        MultiFab::Copy(rhstmp, *prhs, 0, 0, ncomp, 0);
        linop.applyMetricTerm(alev, 0, rhstmp);
        prhs = &rhstmp;
#endif
        linop.solutionResidual(alev, *a_res[alev], *sol[alev], *prhs, crse_bcdata);
        if (alev < finest_amr_lev) {
            linop.reflux(alev, *a_res[alev], *sol[alev], *prhs,
                         *a_res[alev+1], *sol[alev+1], *a_rhs[alev+1]);
            if (linop.isCellCentered()) {
                amrex::average_down(*a_res[alev+1], *a_res[alev], 0, ncomp, amrrr[alev]);
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

    for (int alev = 0; alev < namrlevs; ++alev)
    {
        if (a_in[alev]->nGrow() == 1)
        {
            in[alev] = a_in[alev];
        }
        else
        {
            in_raii[alev].define(a_in[alev]->boxArray(),
                                 a_in[alev]->DistributionMap(),
                                 a_in[alev]->nComp(), 1);
            MultiFab::Copy(in_raii[alev], *a_in[alev], 0, 0, a_in[alev]->nComp(), 0);
            in[alev] = &(in_raii[alev]);
        }
        rh[alev].define(a_in[alev]->boxArray(),
                        a_in[alev]->DistributionMap(),
                        a_in[alev]->nComp(), 0);
        rh[alev].setVal(0.0);
    }

    if (!linop_prepared) {
        linop.prepareForSolve();
        linop_prepared = true;
    }

    const auto& amrrr = linop.AMRRefRatio();

    for (int alev = finest_amr_lev; alev >= 0; --alev) {
        const MultiFab* crse_bcdata = (alev > 0) ? in[alev-1] : nullptr;
        linop.solutionResidual(alev, *out[alev], *in[alev], rh[alev], crse_bcdata);
        if (alev < finest_amr_lev) {
            linop.reflux(alev, *out[alev], *in[alev], rh[alev],
                         *out[alev+1], *in[alev+1], rh[alev+1]);
            if (linop.isCellCentered()) {
                amrex::average_down(*out[alev+1], *out[alev], 0, out[alev]->nComp(), amrrr[alev]);
            }
        }
    }

#if (AMREX_SPACEDIM != 3)
    for (int alev = 0; alev <= finest_amr_lev; ++alev) {
        linop.unapplyMetricTerm(alev, 0, *out[alev]);
    }
#endif

    for (int alev = 0; alev <= finest_amr_lev; ++alev) {
        out[alev]->negate();
    }
}

void
MLMG::averageDownAndSync ()
{
    const auto& amrrr = linop.AMRRefRatio();

    int ncomp = linop.getNComp();

    if (linop.isCellCentered())
    {
        for (int falev = finest_amr_lev; falev > 0; --falev)
        {
            amrex::average_down(*sol[falev], *sol[falev-1], 0, ncomp, amrrr[falev-1]);
        }
    }
    else
    {
        linop.nodalSync(finest_amr_lev, 0, *sol[finest_amr_lev]);

        for (int falev = finest_amr_lev; falev > 0; --falev)
        {
            const auto& fmf = *sol[falev];
            auto&       cmf = *sol[falev-1];

            MultiFab tmpmf(amrex::coarsen(fmf.boxArray(), amrrr[falev-1]),
                           fmf.DistributionMap(), ncomp, 0);
            amrex::average_down(fmf, tmpmf, 0, ncomp, amrrr[falev-1]);
            cmf.ParallelCopy(tmpmf, 0, 0, ncomp);
            linop.nodalSync(falev-1, 0, cmf);
        }
    }
}

Real
MLMG::getNodalSum (int amrlev, int mglev, MultiFab& mf) const
{
    MultiFab one(mf.boxArray(), mf.DistributionMap(), 1, 0);
    one.setVal(1.0);
    const bool local = true;
    Real s1 = linop.xdoty(amrlev, mglev, mf, one, local);
    Real s2 = linop.xdoty(amrlev, mglev, one, one, local);
    ParallelAllReduce::Sum<Real>({s1,s2}, linop.Communicator(amrlev,mglev));
    return s1/s2;
}

void
MLMG::bottomSolveWithHypre (MultiFab& x, const MultiFab& b)
{
#if !defined(AMREX_USE_HYPRE)
    amrex::Abort("bottomSolveWithHypre is called without building with Hypre");
#else

    const int ncomp = linop.getNComp();

    if (hypre_solver == nullptr)  // We should reuse the setup
    {
        const BoxArray& ba = linop.m_grids[0].back();
        const DistributionMapping& dm = linop.m_dmap[0].back();
        const Geometry& geom = linop.m_geom[0].back();
        MPI_Comm comm = linop.BottomCommunicator();

        hypre_solver.reset(new HypreABecLap2(ba, dm, geom, comm));
        hypre_solver->setVerbose(bottom_verbose);

        hypre_solver->setScalars(linop.getAScalar(), linop.getBScalar());

        auto ac = linop.getACoeffs(0, linop.NMGLevels(0)-1);
        if (ac)
        {
            hypre_solver->setACoeffs(*ac);
        }
        else
        {
            MultiFab alpha(ba,dm,ncomp,0);
            alpha.setVal(0.0);
            hypre_solver->setACoeffs(alpha);
        }

        auto bc = linop.getBCoeffs(0, linop.NMGLevels(0)-1);
        if (bc[0])
        {
            hypre_solver->setBCoeffs(bc);
        }
        else
        {
            std::array<MultiFab,AMREX_SPACEDIM> beta;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                beta[idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                                  dm, ncomp, 0);
                beta[idim].setVal(1.0);
            }
            hypre_solver->setBCoeffs(amrex::GetArrOfConstPtrs(beta));
        }

        hypre_bndry.reset(new MLMGBndry(ba, dm, ncomp, geom));
        hypre_bndry->setHomogValues();
        const Real* dx = linop.m_geom[0][0].CellSize();
        int crse_ratio = linop.m_coarse_data_crse_ratio > 0 ? linop.m_coarse_data_crse_ratio : 1;
        RealVect bclocation(AMREX_D_DECL(0.5*dx[0]*crse_ratio,
                                         0.5*dx[1]*crse_ratio,
                                         0.5*dx[2]*crse_ratio));
        hypre_bndry->setLOBndryConds(linop.m_lobc, linop.m_hibc, -1, bclocation);
    }

    hypre_solver->solve(x, b, 1.e-4, -1., bottom_maxiter, *hypre_bndry);

#endif
}

}
