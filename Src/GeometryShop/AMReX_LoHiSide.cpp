

/*
 *      .o.       ooo        ooooo ooooooooo.             ooooooo  ooooo 
 *     .888.      `88.       .888' `888   `Y88.            `8888    d8'  
 *    .8"888.      888b     d'888   888   .d88'  .ooooo.     Y888..8P    
 *   .8' `888.     8 Y88. .P  888   888ooo88P'  d88' `88b     `8888'     
 *  .88ooo8888.    8  `888'   888   888`88b.    888ooo888    .8PY888.    
 * .8'     `888.   8    Y     888   888  `88b.  888    .o   d8'  `888b   
 *o88o     o8888o o8o        o888o o888o  o888o `Y8bod8P' o888o  o88888o 
 *
 */

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
