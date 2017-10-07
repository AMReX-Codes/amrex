
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>

namespace amrex {

MLMG::MLMG (MLLinOp& a_lp)
    : m_lp(a_lp)
{}

MLMG::~MLMG ()
{}

void
MLMG::solve (const Vector<MultiFab*>& sol, const Vector<MultiFab*>& rhs,
             Real a_tol_real, Real a_tol_abs)
{
    AMREX_ASSERT(sol[0]->nGrow() > 0);

    m_lp.prepareForSolve();

    const int namrlevs = m_lp.NAMRLevels();
    AMREX_ASSERT(namrlevs <= sol.size());
    AMREX_ASSERT(namrlevs <= rhs.size());
    
    const auto& amrrr = m_lp.AMRRefRatio();

    for (int falev = namrlevs-1; falev > 0; --falev)
    {
        amrex::average_down(*sol[falev], *sol[falev-1], 0, 1, amrrr[falev-1]);
        amrex::average_down(*rhs[falev], *rhs[falev-1], 0, 1, amrrr[falev-1]);
    }
}

}
