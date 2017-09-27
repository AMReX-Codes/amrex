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

#include "AMReX_BoxIterator.H"

namespace amrex
{
  void BoxIterator::define (const Box& a_bx)
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

  void BoxIterator::setBox(const Box& a_bx)
  {
    define(a_bx);
  }

}
