
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>

namespace amrex {

MLMG::MLMG (MLLinOp& a_lp)
    : linop(a_lp),
      namrlevs(a_lp.NAMRLevels()),
      finest_amr_lev(a_lp.NAMRLevels()-1)
{}

MLMG::~MLMG ()
{}

void
MLMG::solve (const Vector<MultiFab*>& sol, const Vector<MultiFab const*>& a_rhs,
             Real a_tol_real, Real a_tol_abs)
{
    AMREX_ASSERT(sol[0]->nGrow() > 0);
    AMREX_ASSERT(namrlevs <= sol.size());
    AMREX_ASSERT(namrlevs <= a_rhs.size());

    linop.prepareForSolve();

    rhs.resize(namrlevs);
    res.resize(namrlevs);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        const auto& ba = a_rhs[alev]->boxArray();
        const auto& dm = a_rhs[alev]->DistributionMap();
        rhs[alev].define(ba, dm, 1, 0);
        res[alev].define(ba, dm, 1, 0);
        MultiFab::Copy(rhs[alev], *a_rhs[alev], 0, 0, 1, 0);
    }
    
    const auto& amrrr = linop.AMRRefRatio();

    for (int falev = finest_amr_lev; falev > 0; --falev)
    {
        amrex::average_down(*sol[falev], *sol[falev-1], 0, 1, amrrr[falev-1]);
        amrex::average_down( rhs[falev],  rhs[falev-1], 0, 1, amrrr[falev-1]);
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
            oneIter();

            // test convergence
        }
    }
}

void
MLMG::oneIter ()
{
    // if converged?
    //    return

    // uu: correction


    for (int alev = finest_amr_lev; alev > 0; --alev)
    {
        // Update uu[alev] with res[alev] as rhs
        // Example: mg.f90
        // uu = 0
        miniCycle ();

        // Add sol[n] += cor[n]

        // This can be put in a function.
        // compute residual on the coarse amrlevel
        // compute current level residual (using correction)
        // update crse residual with crse/fine residual and restriction of fine residual

        // uu_hold = uu // save it for the up cycle
        // zero uu[alev] ???  don't think this is needed because of interpolation going up.
    }

    // coarest amr level
    {    
        // enforce solvability if appropriate
        // Update uu[0] with res[0] as rhs
        mgCycle (); // Example: mg.f90
        
        // Add sol[0] += cor[0]
    }

    for (int alev = 1; alev <= finest_amr_lev; ++alev)
    {
        // interpolate uu[alev] from uu[alev-1] && uu_hold[alev-1]

        // soln[alev] += uu

        // if (alev < finest_amr_lev) call plus_plus(uu_hold(n), uu(n), nghost(uu(n)))

        // This can be put in a function
        // Interpolate uu[alev-1] to supply boundary conditions for new residual calculation
        // compute defect: res = res - L(uu)

        // uu = 0
        minicycle ();
     
        // soln += uu
    }

    for (int alev = finest_amr_lev; alev > 0; --alev)
    {
        // restrict sol[alev] to sol[alev-1]
    }

    for (int alev = 1; alev <= finest_amr_lev; ++alev)
    {
        // Interpolate soln to supply boundary conditions 
    }
}

void
MLMG::miniCycle ()
{
}

void
MLMG::mgCycle ()
{
}

}

