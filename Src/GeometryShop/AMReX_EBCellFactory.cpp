#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL
#include "EBCellFactory.H"
#include "NamespaceHeader.H"
/***************/
/***************/
EBCellFactory::~EBCellFactory()
{
}
/***************/
/***************/
EBCellFactory::EBCellFactory(const EBISLayout& a_ebisl)
{
  m_ebisl = a_ebisl;
}
/***************/
/***************/
EBCellFAB*
EBCellFactory::create(const Box& a_box, int a_ncomps,
                      const DataIndex& a_dit) const
{
  return new EBCellFAB(m_ebisl[a_dit], a_box, a_ncomps);
}
/***************/
/***************/
#include "NamespaceFooter.H"
