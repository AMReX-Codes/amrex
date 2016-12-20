#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CellEdge.H"

#include "NamespaceHeader.H"

// basically just for constructors

/*
CellEdge::CellEdge(const IntVect& a_cell,
                   const int&     a_dir,
                   const bool*    a_lohi) // note should be bool[SpaceDim-1]
{

  m_cell = a_cell;

  if (a_dir < 0 || a_dir > SpaceDim)
    MayDay::Abort("CellEdge constructor: direction must be between 0 and SpaceDim");
  m_dir = a_dir;

  for (int i = 0; i < SpaceDim-1; i++)
    m_lohi[i] = a_lohi[i];

  // calculate IntVects of nodes (for node-centered data)
  m_node0 = a_cell;
  int ilohi = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    if (idir != a_dir)
    {
      m_node0[idir] += (int) a_lohi[ilohi]; // cast to int and add to the IntVect
      ilohi++;
    }
  }

  m_node1 = m_node0;
  m_node1[a_dir] += 1; // increment along direction of edge

}
*/

CellEdge::CellEdge(const IntVect& a_node0,
                   const IntVect& a_node1)
{

  // initialize
  m_dir = -1;
  int oneNormDiff = 0;

  // look for the incremented direction
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    oneNormDiff += abs( a_node1[idir] - a_node0[idir] );
    if ( (a_node1[idir]-a_node0[idir]) == 1 )
    {
      m_node0 = a_node0; // positive increment in idir direction
      m_node1 = a_node1;
      m_dir = idir;
    }
    else if ( (a_node1[idir]-a_node0[idir]) == -1 )
    {
      m_node0 = a_node1; // negative increment in idir direction
      m_node1 = a_node0; // flip node0 and node1
      m_dir = idir;
    }
  }

  // check for invalid nodes
  if ( m_dir==-1 || (oneNormDiff != 1) )
    MayDay::Abort("CellEdge: Constructor: Cannot make CellEdge with non-adjacent nodes");

  // check for invalid cell
  /*
  IntVect tmp0 = a_node0 - a_cell;
  IntVect tmp1 = a_node1 - a_cell;
  for (int idir = 0; idir < SpaceDim; idir++)
    if ( tmp0[idir]>1 || tmp0[idir]<0 || tmp1[idir]>1 || tmp1[idir]<0 )
      MayDay::Abort("Cannot make CellEdge with a cell and two nodes that are not part of the cell");
  */

  // set lohi
  /*
  for (int idir = 0; idir < SpaceDim-1; idir++)
  {
    if (idir != m_dir)
    {
      m_lohi[idir] = abs(tmp0[idir])==1; // for all directions other than m_dir, basically cast tmp0 to boolean
    }
  }
  */

}

CellEdge::CellEdge(const IntVect& a_node0,
                   const int      a_dir)
{
  if (Abs(a_dir) > (SpaceDim-1))
  {
    MayDay::Abort("CellEdge: Contructor: invalid direction");
  }
  else if (a_dir<0)
  {
    // negative direction - flip
    m_dir = abs(a_dir);
    m_node0 = a_node0; m_node0[a_dir]--; // move back one
    m_node1 = a_node0;
  }
  else
  {
    m_dir = a_dir;
    m_node0 = a_node0;
    m_node1 = a_node0; m_node1[a_dir]++; // move forward one
    // don't specify a cell or a lohi
  }
}

void CellEdge::shift(const int a_dir)
{
  shift(a_dir,1);
}

void CellEdge::shift(const int a_dir,
                     const int a_offset)
{
  m_node0[a_dir]+=a_offset;
  m_node1[a_dir]+=a_offset;
}

CellEdge::~CellEdge()
{
}

// test CellEdge equality (for unordered_map) with _both_ nodes _and_ the direction
bool CellEdge::operator==(const CellEdge& a_edge) const
{
  return ((this->m_node0==a_edge.m_node0) && (this->m_node1==a_edge.m_node1) && (this->m_dir==a_edge.m_dir));
}

#include "NamespaceFooter.H"

