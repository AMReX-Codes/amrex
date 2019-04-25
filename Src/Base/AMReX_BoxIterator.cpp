#include "AMReX_BoxIterator.H"

namespace amrex
{
  void BoxIterator::define (const Box& a_bx) noexcept
  {
    if (a_bx.ok() && a_bx.smallEnd() <= a_bx.bigEnd())
      {
        m_current = a_bx.smallEnd();
        m_boxLo   = a_bx.smallEnd();
        m_boxHi   = a_bx.bigEnd();
      }
    else
      {
        m_current = IntVect::TheUnitVector();
        m_boxLo   = IntVect::TheUnitVector();
        m_boxHi   = IntVect::TheZeroVector();
      }
  }

  void BoxIterator::setBox(const Box& a_bx) noexcept
  {
    define(a_bx);
  }

}
