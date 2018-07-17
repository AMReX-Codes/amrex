#include "AMReX_BLassert.H"
#include "AMReX_LoHiSide.H"

namespace amrex
{
  Side::LoHiSide flip(const Side::LoHiSide& a_side)
  {
    AMREX_ASSERT((a_side == Side::Lo) || (a_side == Side::Hi));

    return (a_side == Side::Lo) ? Side::Hi : Side::Lo;
  }

  Side::LoHiSide Side::flip(const Side::LoHiSide& a_side)
  {
    AMREX_ASSERT((a_side == Side::Lo) || (a_side == Side::Hi));

    return (a_side == Side::Lo) ? Side::Hi : Side::Lo;
  }

  int sign(const Side::LoHiSide& a_side)
  {
    AMREX_ASSERT((a_side == Side::Lo) || (a_side == Side::Hi));

    return (a_side == Side::Lo) ? -1 : 1;
  }

  SideIterator::SideIterator()
    :
    m_current(-1)
  {
    reset();
  }

  void SideIterator::begin()
  {
    m_current = 0;
  }

  void SideIterator::next()
  {
    ++m_current;
  }

  void SideIterator::operator ++ ()
  {
    ++m_current;
  }

  Side::LoHiSide SideIterator::operator () () const
  {
    switch (m_current)
      {
      case 0:
        return Side::Lo;
        //break;

      case 1:
        return Side::Hi;
        //break;

      default:
        return Side::Invalid;
        //break;
      }
  }

  bool SideIterator::ok() const
  {
    return ((m_current > -1) && (m_current < Side::NUMSIDES));
  }
}
