#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CELLEDGE_H_
#define _CELLEDGE_H_


#include "MayDay.H"
#include "RealVect.H"
#include "IntVect.H"

#include "NamespaceHeader.H"

class CellEdge
{
public:

  /*
   * Don't use this kind of construction anymore
  CellEdge(const IntVect& a_cell,
           const int&     a_dir,
           const bool*    a_lohi);
           */

  CellEdge(const IntVect& a_node0,
           const IntVect& a_node1);

  CellEdge(const IntVect& a_node0,
           const int      a_dir);

  ~CellEdge();

  void shift(const int a_dir);

  void shift(const int a_dir,
             const int a_offset);

  // test CellEdge equality (for unordered_map) with _both_ nodes _and_ the direction
  bool operator==(const CellEdge& a_edge) const;

  //IntVect m_cell; // don't use this anymore
  int     m_dir;
  //bool    m_lohi[SpaceDim-1]; // don't use this anymore
  IntVect m_node0;
  IntVect m_node1;

private:

  CellEdge()
  {
    MayDay::Abort("CellEdge uses strong construction");
  }

  void operator=(const CellEdge& a_inputCellEdge)
  {
    MayDay::Abort("CellEdge doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
