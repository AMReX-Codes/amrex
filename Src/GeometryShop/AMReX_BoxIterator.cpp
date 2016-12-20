#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"
#include "NamespaceHeader.H"

void BoxIterator::define (const Box& a_bx)
{
  // this algorithm breaks for screwy empty boxes
  if (a_bx.smallEnd() <= a_bx.bigEnd())
  {
    m_current = a_bx.smallEnd();
    m_boxLo   = a_bx.smallEnd();
    m_boxHi   = a_bx.bigEnd();
  }
  else
  {
    m_current = IntVect::Unit;
    m_boxLo   = IntVect::Unit;
    m_boxHi   = IntVect::Zero;
  }
}

void BoxIterator::setBox(const Box& a_bx)
{
  define(a_bx);
}

#include "NamespaceFooter.H"
