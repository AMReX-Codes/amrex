
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MG_F.H>
#include <AMReX_MLCGSolver.H>

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

void
MLMG::oneIter (int iter)
{
    for (int alev = finest_amr_lev; alev > 0; --alev)
    {
        miniCycle(alev);

        MultiFab::Add(*sol[alev], *cor[alev][0], 0, 0, 1, 0);

        computeResWithCrseSolFineCor(alev-1,alev);

        if (alev != finest_amr_lev) {
            std::swap(cor_hold[alev][0], cor[alev][0]); // save it for the up cycle
        }
    }

    // coarest amr level
    {    
        // enforce solvability if appropriate

        if (iter < nfcycles) {
            mgFcycle ();
        } else {
            mgVcycle (0, 0);
        }

        MultiFab::Add(*sol[0], *cor[0][0], 0, 0, 1, 0);
    }

    for (int alev = 1; alev <= finest_amr_lev; ++alev)
    {
        interpCorrection(alev);

        MultiFab::Add(*sol[alev], *cor[alev][0], 0, 0, 1, 0);

        if (alev != finest_amr_lev) {
            MultiFab::Add(*cor_hold[alev][0], *cor[alev][0], 0, 0, 1, 0);
        }

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

void
MLMG::computeMLResidual (int amrlevmax, bool bndryregister_updated)
{
    const int mglev = 0;
    for (int alev = amrlevmax; alev >= 0; --alev) {
        if (alev > 0 && !bndryregister_updated) {
            linop.updateSolBC(alev, *sol[alev-1]);
        }
        linop.residual(alev, 0, res[alev][mglev], *sol[alev], rhs[alev], BCMode::Inhomogeneous);
        if (alev < finest_amr_lev) {
            linop.reflux(alev, res[alev][mglev], *sol[alev], *sol[alev+1]);
        }
    }
}

void
MLMG::computeResidual (int alev)
{
    MultiFab& x = *sol[alev];
    const MultiFab& b = rhs[alev];
    MultiFab& r = res[alev][0];

    if (alev > 0) {
        linop.updateSolBC(alev, *sol[alev-1]); // TODO: don't have to do this everytime. - wqz
    }
    linop.residual(alev, 0, r, x, b, BCMode::Inhomogeneous);
}

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
    
    if (calev > 0) {
        linop.updateSolBC(calev, *sol[calev-1]);
    }
    linop.residual(calev, 0, crse_res, crse_sol, crse_rhs, BCMode::Inhomogeneous);

    linop.residual(falev, 0, fine_rescor, fine_cor, fine_res, BCMode::Homogeneous);
    MultiFab::Copy(fine_res, fine_rescor, 0, 0, 1, 0);

    linop.reflux(calev, crse_res, crse_sol, fine_sol);

    const int amrrr = linop.AMRRefRatio(calev);
    amrex::average_down(fine_res, crse_res, 0, 1, amrrr);
}

void
MLMG::computeResWithCrseCorFineCor (int falev)
{
    const MultiFab& crse_cor = *cor[falev-1][0];

    MultiFab& fine_cor = *cor[falev][0];
    MultiFab& fine_res = res[falev][0];
    MultiFab& fine_rescor = rescor[falev][0];

    linop.updateCorBC(falev, crse_cor);
    // fine_rescor = fine_res - L(fine_cor)
    linop.correctionResidual(falev, fine_rescor, fine_cor, fine_res);
    MultiFab::Copy(fine_res, fine_rescor, 0, 0, 1, 0);
}

void
MLMG::miniCycle (int amrlev)
{
    const int mglev = 0;
    mgVcycle(amrlev, mglev);
}

void
MLMG::mgVcycle (int amrlev, int mglev)
{
    MultiFab& x = *cor[amrlev][mglev];
    MultiFab& b = res[amrlev][mglev];
    MultiFab& r = rescor[amrlev][mglev];
    const int mg_bottom_lev = linop.NMGLevels(amrlev) - 1;

    x.setVal(0.0);

    if (mglev == mg_bottom_lev)
    {
        if (amrlev == 0)
        {
            bottomSolve();
        }
        else
        {
            for (int i = 0; i < nuf; ++i) {
                linop.smooth(amrlev, mglev, x, b, BCMode::Homogeneous);
            }
        }
    }
    else
    {
        for (int i = 0; i < nu1; ++i) {
            linop.smooth(amrlev, mglev, x, b, BCMode::Homogeneous);
        }
        
        computeResOfCorrection(amrlev, mglev);

        const int ratio = 2;
        amrex::average_down(r, res[amrlev][mglev+1], 0, 1, ratio);

        for (int i = 0; i < nu0; ++i) {
            mgVcycle(amrlev, mglev+1);
        }

        addInterpCorrection(amrlev, mglev);

        for (int i = 0; i < nu2; ++i) {
            linop.smooth(amrlev, mglev, x, b, BCMode::Homogeneous);
        }
    }    
}

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

    cor[amrlev][mg_bottom_lev]->setVal(0.0);
    bottomSolve();

    for (int mglev = mg_bottom_lev-1; mglev >= 0; --mglev)
    {
        cor[amrlev][mglev]->setVal(0.0);  // interp next line performs add not assignment
        addInterpCorrection(amrlev, mglev); // interp from mglev+1 to mglev

        computeResOfCorrection(amrlev, mglev);
        MultiFab::Copy(res[amrlev][mglev], rescor[amrlev][mglev], 0,0,1,0);

        std::swap(cor[amrlev][mglev], cor_hold[amrlev][mglev]);
        mgVcycle(amrlev, mglev);
        MultiFab::Add(*cor[amrlev][mglev], *cor_hold[amrlev][mglev], 0, 0, 1, 0);
    }
}

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

void
MLMG::computeResOfCorrection (int amrlev, int mglev)
{
    MultiFab& x = *cor[amrlev][mglev];
    const MultiFab& b = res[amrlev][mglev];
    MultiFab& r = rescor[amrlev][mglev];
    linop.residual(amrlev, mglev, r, x, b, BCMode::Homogeneous);
}

void
MLMG::bottomSolve ()
{
    Real bottom_start_time = ParallelDescriptor::second();

    const int amrlev = 0;
    const int mglev = linop.NMGLevels(amrlev) - 1;
    MultiFab& x = *cor[amrlev][mglev];
    MultiFab& b = res[amrlev][mglev];
    
    if (bottom_solver == BottomSolver::smoother)
    {
        for (int i = 0; i < nuf; ++i) {
            linop.smooth(amrlev, mglev, x, b, BCMode::Homogeneous);
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
            linop.smooth(amrlev, mglev, x, b, BCMode::Homogeneous);
        }
    }
    timer[bottom_time] += ParallelDescriptor::second() - bottom_start_time;
}

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
        }
    }

    cor_hold.resize(namrlevs-1);
    {
        const int alev = 0;
        const int nmglevs = linop.NMGLevels(alev);
        cor_hold[alev].resize(nmglevs);
        for (int mglev = 0; mglev < nmglevs-1; ++mglev)
        {
            cor_hold[alev][mglev].reset(new MultiFab(cor[alev][mglev]->boxArray(),
                                                     cor[alev][mglev]->DistributionMap(),
                                                     nc, ng));
        }
    }
    for (int alev = 1; alev < finest_amr_lev; ++alev)
    {
        cor_hold[alev].resize(1);
        cor_hold[alev][0].reset(new MultiFab(cor[alev][0]->boxArray(),
                                             cor[alev][0]->DistributionMap(),
                                             nc, ng));
    }

    buildFineMask();
}

}
