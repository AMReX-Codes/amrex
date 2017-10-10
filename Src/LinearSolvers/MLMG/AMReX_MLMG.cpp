
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

    rhs_o.resize(namrlevs);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        rhs_o[alev].define(a_rhs[alev]->boxArray(), a_rhs[alev]->DistributionMap(), 1, 0);
        MultiFab::Copy(rhs_o[alev], *a_rhs[alev], 0, 0, 1, 0);        
    }

    const auto& amrrr = linop.AMRRefRatio();

    for (int falev = finest_amr_lev; falev > 0; --falev)
    {
        amrex::average_down(*sol[falev], *sol[falev-1], 0, 1, amrrr[falev-1]);
        amrex::average_down(rhs_o[falev], rhs_o[falev-1], 0, 1, amrrr[falev-1]);
    }

    const int nc = 1;
    int ng = 0;
    linop.make(rhs_c, nc, ng);
    linop.make(res, nc, ng);
    ng = 1;
    linop.make(cor, nc, ng);

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
            oneIter(sol);

            // test convergence
        }
    }
}

void
MLMG::oneIter (const Vector<MultiFab*>& sol)
{
    // if converged?
    //    return

    for (int alev = finest_amr_lev; alev > 0; --alev)
    {
        computeResidual(alev);

        miniCycle(alev);

        MultiFab::Add(*sol[alev], cor[alev][0], 0, 0, 1, 0);

        // This can be put in a function.
        // compute residual on the coarse amrlevel
        // compute current level residual (using correction)
        // update crse residual with crse/fine residual and restriction of fine residual

        if (alev != finest_amr_lev) {
//            MultiFab::Copy(cor_hold[alev], cor[alev][0], 0, 0, 1, 0); // save it for the up cycle
        }
    }

    // coarest amr level
    {    
        computeResidual(0);

        // enforce solvability if appropriate

        mgCycle ();

        MultiFab::Add(*sol[0], cor[0][0], 0, 0, 1, 0);
        
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
        miniCycle (alev);
     
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
MLMG::computeResidual (int alev)
{

}

void
MLMG::miniCycle (int alev)
{
    auto& xs = cor[alev];
    const auto& bs = res[alev];

    for (auto& x : xs) x.setVal(0.0);

    for (int i = 0; i < nu1; ++i) {
        int mglev = 0;
        // linop.smooth(alev, mglev, xs[mglev], bs[mglev]);
    }

    // for ref ratio of 4 ...
    
}

void
MLMG::mgCycle ()
{
    auto& xs = cor[0];
    const auto& bs = res[0];

    for (auto& x : xs) x.setVal(0.0);

    for (int i = 0; i < nu1; ++i) {
        int mglev = 0;
        // linop.smooth(alev, mglev, xs[mglev], bs[mglev]);
    }

    // compute defect 
    
}

}

