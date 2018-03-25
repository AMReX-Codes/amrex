#include "AMReX_EBSimpleSolver.H"

namespace amrex
{

  EBSimpleSolver::EBSimpleSolver()
  {
    m_isDefined = false;

    m_operator = NULL;
    m_homogeneous = true;

    m_numSmooths = 4;
  }

  EBSimpleSolver::~EBSimpleSolver()
  {
  }

  void EBSimpleSolver::setHomogeneous(bool a_homogeneous)
  {
    m_homogeneous = a_homogeneous;
  }

  void EBSimpleSolver::define(LinearOp<FabArray<EBCellFAB> >* a_operator,
                              bool                             a_homogeneous)
  {
    m_isDefined = true;

    m_operator = dynamic_cast <MGLevelOp<FabArray<EBCellFAB> >* > (a_operator);
    if (m_operator == NULL)
    {
      amrex::Error("EBSimpleSolver::define - operator not an MGLevelOp");
    }

    m_homogeneous = a_homogeneous;
  }

  void EBSimpleSolver::setNumSmooths(const int& a_numSmooths)
  {
    m_numSmooths = a_numSmooths;
  }

  void EBSimpleSolver::solve(FabArray<EBCellFAB>&       a_phi,
                             const FabArray<EBCellFAB>& a_rhs)
  {
    BL_ASSERT(m_isDefined);

    m_operator->relax(a_phi,a_rhs,m_numSmooths);
  }
}
