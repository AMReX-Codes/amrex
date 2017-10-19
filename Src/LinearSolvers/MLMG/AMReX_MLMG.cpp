
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MG_F.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_BC_TYPES.H>

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

void
MLMG::solve (const Vector<MultiFab*>& a_sol, const Vector<MultiFab const*>& a_rhs,
             Real a_tol_real, Real a_tol_abs)
{
    Real solve_start_time = ParallelDescriptor::second();

    prepareForSolve(a_sol, a_rhs);

    bool bndryregister_updated = true;  // We did that when setBC.
    computeMLResidual(finest_amr_lev, bndryregister_updated);

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
    const Real res_target = std::max(a_tol_abs, a_tol_real*max_norm);

    if (resnorm0 <= res_target)
    {
        if (verbose >= 1) {
            amrex::Print() << "MLMG: No iterations needed\n";
        }
    }
    else
    {
        Real iter_start_time = ParallelDescriptor::second();
        for (int iter = 0; iter < max_iters; ++iter)
        {
            oneIter(iter);

            bool converged = false;
            Real composite_norminf;

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
                                   << " composite resid/" << norm_name << " = "
                                   << composite_norminf/max_norm << "\n";
                }

                break;
            }
        }
        timer[iter_time] = ParallelDescriptor::second() - iter_start_time;
    }

    timer[solve_time] = ParallelDescriptor::second() - solve_start_time;
    if (verbose >= 1) {
        ParallelDescriptor::ReduceRealMax(timer.data(), timer.size());
        amrex::Print() << "MLMG: Timers: Solve = " << timer[solve_time]
                       << " Iter = " << timer[iter_time]
                       << " Bottom = " << timer[bottom_time] << "\n";
    }
}

// in  : Residual (res) on the finest AMR level
// out : sol on all AMR levels
void
MLMG::oneIter (int iter)
{
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

    const auto& amrrr = linop.AMRRefRatio();
    for (int falev = finest_amr_lev; falev > 0; --falev)
    {
        amrex::average_down(*sol[falev], *sol[falev-1], 0, 1, amrrr[falev-1]);
    }
}

// Compute multi-level Residual (res) up to amrlevmax.
void
MLMG::computeMLResidual (int amrlevmax, bool bndryregister_updated)
{
    const int mglev = 0;
    for (int alev = amrlevmax; alev >= 0; --alev) {
        const MultiFab* crse_bcdata = nullptr;
        if (alev > 0 && !bndryregister_updated) {
            crse_bcdata = sol[alev-1];
        }
        linop.solutionResidual(alev, res[alev][mglev], *sol[alev], rhs[alev], crse_bcdata);
        if (alev < finest_amr_lev) {
            linop.reflux(alev, res[alev][mglev], *sol[alev], *sol[alev+1]);
        }
    }
}

// Compute single AMR level residual without masking.
void
MLMG::computeResidual (int alev)
{
    MultiFab& x = *sol[alev];
    const MultiFab& b = rhs[alev];
    MultiFab& r = res[alev][0];

    const MultiFab* crse_bcdata = nullptr;
    if (alev > 0) {
        crse_bcdata = sol[alev-1]; // TODO: don't have to do this everytime. - wqz
    }
    linop.solutionResidual(alev, r, x, b, crse_bcdata);
}

// Compute coarse AMR level composite residual with coarse solution and fine correction
void
MLMG::computeResWithCrseSolFineCor (int calev, int falev)
{
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

    linop.reflux(calev, crse_res, crse_sol, fine_sol);

    const int amrrr = linop.AMRRefRatio(calev);
    amrex::average_down(fine_res, crse_res, 0, 1, amrrr);
}

// Compute fine AMR level residual fine_res = fine_res - L(fine_cor) with coarse providing BC.
void
MLMG::computeResWithCrseCorFineCor (int falev)
{
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
    const int mglev = 0;
    mgVcycle(amrlev, mglev);
}

// in   : Residual (res) 
// out  : Correction (cor) from bottom to this function's local top
void
MLMG::mgVcycle (int amrlev, int mglev_top)
{
    const int mglev_bottom = linop.NMGLevels(amrlev) - 1;

    for (int mglev = mglev_top; mglev < mglev_bottom; ++mglev)
    {
        cor[amrlev][mglev]->setVal(0.0);
        for (int i = 0; i < nu1; ++i) {
            linop.smooth(amrlev, mglev, *cor[amrlev][mglev], res[amrlev][mglev]);
        }

        // rescor = res - L(cor)
        computeResOfCorrection(amrlev, mglev);

        // res_crse = R(rescor_fine); this provides res/b to the level below
        const int ratio = 2;
        amrex::average_down(rescor[amrlev][mglev], res[amrlev][mglev+1], 0, 1, ratio);
    }

    if (amrlev == 0)
    {
        bottomSolve();
    }
    else
    {
        cor[amrlev][mglev_bottom]->setVal(0.0);
        for (int i = 0; i < nuf; ++i) {
            linop.smooth(amrlev, mglev_bottom, *cor[amrlev][mglev_bottom], res[amrlev][mglev_bottom]);
        }
    }
    
    for (int mglev = mglev_bottom-1; mglev >= 0; --mglev)
    {
        // cor_fine += I(cor_crse)
        addInterpCorrection(amrlev, mglev);
        for (int i = 0; i < nu2; ++i) {
            linop.smooth(amrlev, mglev, *cor[amrlev][mglev], res[amrlev][mglev]);
        }
    }
}

// FMG cycle on the coarest AMR level.
// in:  Residual on the top MG level (i.e., 0)
// out: Correction (cor) on all MG levels
void
MLMG::mgFcycle ()
{
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
    const MultiFab& crse_cor = *cor[alev-1][0];
    MultiFab& fine_cor = *cor[alev][0];

    BoxArray ba = fine_cor.boxArray();
    const int amrrr = linop.AMRRefRatio(alev-1);
    IntVect refratio{amrrr};
    ba.coarsen(refratio);
    
    MultiFab cfine(ba, fine_cor.DistributionMap(), 1, 0);
    cfine.copy(crse_cor);

    Geometry g1, g2;
    Vector<BCRec> bcr;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(fine_cor, MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
    {
        amrex::pc_interp.interp(cfine[mfi], 0, fine_cor[mfi], 0, 1,
                                mfi.tilebox(), refratio, g1, g2, bcr, 0, 0);
    }
}

// Interpolate correction between MG levels
// inout: Correction (cor) on coarse MG lev.  (out due to FillBoundary)
// out  : Correction (cor) on fine MG lev.
void
MLMG::interpCorrection (int alev, int mglev)
{
    MultiFab& crse_cor = *cor[alev][mglev+1];
    MultiFab& fine_cor = *cor[alev][mglev  ];
    const Geometry& crse_geom = linop.Geom(alev,mglev+1);
    const Geometry& fine_geom = linop.Geom(alev,mglev  );
    crse_cor.FillBoundary(crse_geom.periodicity());
    Vector<BCRec> bcr(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (crse_geom.isPeriodic(idim)) {
            bcr[0].setLo(idim,BCType::int_dir);
            bcr[0].setHi(idim,BCType::int_dir);
        } else {
            bcr[0].setLo(idim,BCType::ext_dir);
            bcr[0].setHi(idim,BCType::ext_dir);
        }
    }
    IntVect refratio{2};
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(fine_cor, MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
    {
        amrex::cell_cons_interp.interp(crse_cor[mfi], 0, fine_cor[mfi], 0, 1,
                                       mfi.tilebox(), refratio, crse_geom, fine_geom,
                                       bcr, 0, 0);
    }    
}

// (Fine MG level correction) += I(Coarse MG level correction)
void
MLMG::addInterpCorrection (int alev, int mglev)
{
    const MultiFab& crse_cor = *cor[alev][mglev+1];
    MultiFab&       fine_cor = *cor[alev][mglev  ];
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(crse_cor,true); mfi.isValid(); ++mfi)
    {
        const Box&         bx = mfi.tilebox();
        const int          nc = fine_cor.nComp();
        const FArrayBox& cfab = crse_cor[mfi];
        FArrayBox&       ffab = fine_cor[mfi];

        FORT_INTERP(ffab.dataPtr(),
                    ARLIM(ffab.loVect()), ARLIM(ffab.hiVect()),
                    cfab.dataPtr(),
                    ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
                    bx.loVect(), bx.hiVect(), &nc);
    }    
}

// Compute rescor = res - L(cor)
// in   : res
// inout: cor (out due to FillBoundary in linop.correctionResidual)
// out  : rescor
void
MLMG::computeResOfCorrection (int amrlev, int mglev)
{
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
    Real bottom_start_time = ParallelDescriptor::second();

    const int amrlev = 0;
    const int mglev = linop.NMGLevels(amrlev) - 1;
    MultiFab& x = *cor[amrlev][mglev];
    MultiFab& b = res[amrlev][mglev];
    
    x.setVal(0.0);

    if (bottom_solver == BottomSolver::smoother)
    {
        for (int i = 0; i < nuf; ++i) {
            linop.smooth(amrlev, mglev, x, b);
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
        cg_solver.solve(x, b, cg_rtol, cg_atol);

        for (int i = 0; i < nub; ++i) {
            linop.smooth(amrlev, mglev, x, b);
        }
    }
    timer[bottom_time] += ParallelDescriptor::second() - bottom_start_time;
}

// Compute single-level masked inf-norm of Residual (res).
Real
MLMG::ResNormInf (int alev, bool local)
{
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
    const int mglev = 0;
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
    }

    const auto& amrrr = linop.AMRRefRatio();
    for (int falev = finest_amr_lev; falev > 0; --falev)
    {
        amrex::average_down(*sol[falev], *sol[falev-1], 0, 1, amrrr[falev-1]);
        amrex::average_down( rhs[falev],  rhs[falev-1], 0, 1, amrrr[falev-1]);
    }

    const int nc = 1;
    int ng = 0;
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

}
