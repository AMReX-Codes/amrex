
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MG_F.H>

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
    AMREX_ASSERT(a_sol[0]->nGrow() > 0);
    AMREX_ASSERT(namrlevs <= a_sol.size());
    AMREX_ASSERT(namrlevs <= a_rhs.size());

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
    linop.make(cor, nc, ng);

    cor_hold.resize(namrlevs-1);
    for (int alev = 1; alev < finest_amr_lev; ++alev)
    {
        cor_hold[alev].define(cor[alev][0].boxArray(),
                              cor[alev][0].DistributionMap(),
                              cor[alev][0].nComp(),
                              cor[alev][0].nGrow());
    }

    //
    // TODO: We need to fill the fine amr level ghost cells by interploating from the coarse
    // Example: ml_prolongation.f90
    // We can use the fillpatch meta data.
    //

    //
    // TODO: enforce solvability if appropriate
    //
    
    for (int alev = finest_amr_lev; alev >= 0; --alev)
    {
        //
        // TODO: compute residues res given sol and rhs
        // Example: compute_defect.f90
        // 
        if (alev < finest_amr_lev)
        {
            //
            // TODO: (1) Fix crse level residue at crse/fine boundary
            //       (2) retrict the fine res down to this crse level
            // Example: cc_ml_resid.f90
            //
        }
    }

    //
    // TODO: compute the intial inf-norm of res and rhs
    // Example: ml_norm.f90
    //

    //
    // TODO: need a multi-levle covergence test function
    // Example: ml_cc.f90
    // 
    if (false) // replace with the covergence test
    {
        if (verbose >= 1) {
            amrex::Print() << "MLMG: No iterations needed\n";
        }
    }
    else
    {
        for (int iter = 0; iter < max_iters; ++iter)
        {
            oneIter(iter);
        }
    }
}

void
MLMG::oneIter (int iter)
{
    // if converged?
    //    return

    computeResidual(finest_amr_lev);

    for (int alev = finest_amr_lev; alev > 0; --alev)
    {
        miniCycle(alev);

        MultiFab::Add(*sol[alev], cor[alev][0], 0, 0, 1, 0);

        computeResWithCrseSolFineCor(alev-1,alev);

        if (alev != finest_amr_lev) {
            MultiFab::Copy(cor_hold[alev], cor[alev][0], 0, 0, 1, 1); // save it for the up cycle
        }
    }

    // coarest amr level
    {    
        // enforce solvability if appropriate

        mgCycle ();

        MultiFab::Add(*sol[0], cor[0][0], 0, 0, 1, 0);
    }

    for (int alev = 1; alev <= finest_amr_lev; ++alev)
    {
        interpCorrection(alev);

        MultiFab::Add(*sol[alev], cor[alev][0], 0, 0, 1, 0);

        if (alev != finest_amr_lev) {
            MultiFab::Add(cor_hold[alev], cor[alev][0], 0, 0, 1, 1);
        }

        computeResWithCrseCorFineCor(alev);

        miniCycle (alev);

        MultiFab::Add(*sol[alev], cor[alev][0], 0, 0, 1, 0);

        if (alev != finest_amr_lev) {
            MultiFab::Add(cor[alev][0], cor_hold[alev], 0, 0, 1, 1);            
        }
    }

    const auto& amrrr = linop.AMRRefRatio();
    for (int falev = finest_amr_lev; falev > 0; --falev)
    {
        amrex::average_down(*sol[falev], *sol[falev-1], 0, 1, amrrr[falev-1]);
    }

    // ...
    if (verbose > 1) {
        for (int alev = 0; alev <= finest_amr_lev; ++alev) {
            Real resmax = res[alev][0].norm0();
            amrex::Print() << "MLMG: Iter " << iter << " Level " << alev 
                           << " max resid " << resmax << "\n";
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
        linop.updateSolBC(alev, *sol[alev-1]);
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
    MultiFab& fine_cor = cor[falev][0];
    MultiFab& fine_res = res[falev][0];
    MultiFab& fine_rescor = rescor[falev][0];
    
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
    const MultiFab& crse_cor = cor[falev-1][0];

    MultiFab& fine_cor = cor[falev][0];
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
MLMG::mgCycle ()
{
    const int amrlev = 0;
    const int mglev = 0;
    if (cycle_type == Cycle::Vcycle)
    {
        mgVcycle (amrlev, mglev);
    }
    else if (cycle_type == Cycle::Wcycle)
    {
        amrex::Abort("not implemented");
    }
    else if (cycle_type == Cycle::Fcycle)
    {
        amrex::Abort("not implemented");
    }
}

void
MLMG::mgVcycle (int amrlev, int mglev)
{
    MultiFab& x = cor[amrlev][mglev];
    MultiFab& b = res[amrlev][mglev];
    MultiFab& r = rescor[amrlev][mglev];
    const int mg_bottom_lev = linop.NMGLevels(amrlev) - 1;

    x.setVal(0.0);

    if (mglev == mg_bottom_lev)
    {
        if (amrlev == 0)
        {
            bottomSolve ();
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

        const int refratio = 2;
        amrex::average_down(r, res[amrlev][mglev+1], 0, 1, refratio);

        for (int i = 0; i < nu0; ++i) {
            mgVcycle(amrlev, mglev+1);
        }

        interpCorrection(amrlev, mglev);

        for (int i = 0; i < nu2; ++i) {
            linop.smooth(amrlev, mglev, x, b, BCMode::Homogeneous);
        }
    }    
}

void
MLMG::interpCorrection (int alev)
{
    const MultiFab& crse_cor = cor[alev-1][0];
    MultiFab& fine_cor = cor[alev][0];

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
MLMG::interpCorrection (int alev, int mglev)
{
    const MultiFab& crse_cor = cor[alev][mglev+1];
    MultiFab&       fine_cor = cor[alev][mglev  ];
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
    MultiFab& x = cor[amrlev][mglev];
    const MultiFab& b = res[amrlev][mglev];
    MultiFab& r = rescor[amrlev][mglev];
    linop.residual(amrlev, mglev, r, x, b, BCMode::Homogeneous);
}

void
MLMG::bottomSolve ()
{
    const int amrlev = 0;
    const int mglev = linop.NMGLevels(amrlev) - 1;
    MultiFab& x = cor[amrlev][mglev];
    MultiFab& b = res[amrlev][mglev];
    
    if (bottom_solver == BottomSolver::smoother)
    {
        for (int i = 0; i < nuf; ++i) {
            linop.smooth(amrlev, mglev, x, b, BCMode::Homogeneous);
        }
    }
    else
    {
        if (bottom_solver == BottomSolver::bicgstab)
        {
            amrex::Abort("MLMG:: bicgstab not implemented");
        }
        else if (bottom_solver == BottomSolver::cg)
        {
            amrex::Abort("MLMG:: cg not implemented");
        }
        else if (bottom_solver == BottomSolver::cabicgstab)
        {
            amrex::Abort("MLMG:: cabicgstab not implemented");
        }

        for (int i = 0; i < nub; ++i) {
            linop.smooth(amrlev, mglev, x, b, BCMode::Homogeneous);
        }
    }
}

}
