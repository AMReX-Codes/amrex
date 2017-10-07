
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>

namespace amrex {

MLMG::MLMG (MLLinOp& a_lp)
    : linop(a_lp)
{}

MLMG::~MLMG ()
{}

void
MLMG::solve (const Vector<MultiFab*>& sol, const Vector<MultiFab*>& rhs,
             Real a_tol_real, Real a_tol_abs)
{
    AMREX_ASSERT(sol[0]->nGrow() > 0);

    linop.prepareForSolve();

    const int namrlevs = linop.NAMRLevels();
    const int finest_amr_lev = namrlevs-1;
    AMREX_ASSERT(namrlevs <= sol.size());
    AMREX_ASSERT(namrlevs <= rhs.size());

    Array<MultiFab> res(namrlevs);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        const auto& ba = rhs[alev]->boxArray();
        const auto& dm = rhs[alev]->DistributionMap();
        res[alev].define(ba, dm, 1, 0);
    }
    
    const auto& amrrr = linop.AMRRefRatio();

    for (int falev = finest_amr_lev; falev > 0; --falev)
    {
        amrex::average_down(*sol[falev], *sol[falev-1], 0, 1, amrrr[falev-1]);
        amrex::average_down(*rhs[falev], *rhs[falev-1], 0, 1, amrrr[falev-1]);
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
    if (true) // replace with the covergence test
    {
        if (verbose >= 1) {
            amrex::Print() << "MLMG: No iterations needed\n";
        }
    }
    else
    {
        for (int iter = 0; iter < max_iters; ++iter)
        {
            // OneSolve
        }
    }
}

}

#if 0

void OneSolve ()
{
    const int namrlevs = linop.NAMRLevels();
    const int finest_amr_lev = namrlevs-1;

    uu = 0;
    uu_hold = 0;

    // Down the V-cycle
    for (int alev = finest_amr_lev; alev >= 0; --alev)
    {
        
    }
}

#endif
