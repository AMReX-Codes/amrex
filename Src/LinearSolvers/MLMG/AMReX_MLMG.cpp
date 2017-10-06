
#include <AMReX_MLMG.H>

namespace amrex {

MLMG::MLMG (MLLinOp& a_lp)
    : m_lp(a_lp)
{}

MLMG::~MLMG ()
{}

void
MLMG::solve (const Vector<MultiFab*>& sol, const Vector<MultiFab const*>& rhs,
             Real a_tol_real, Real a_tol_abs)
{
    AMREX_ASSERT(sol[0]->nGrow() > 0);
}

}
