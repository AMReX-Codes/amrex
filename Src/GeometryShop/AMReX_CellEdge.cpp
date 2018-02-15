
/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#include "AMReX_CellEdge.H"


namespace amrex
{
// basically just for constructors


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
      amrex::Abort("CellEdge: Constructor: Cannot make CellEdge with non-adjacent nodes");

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
}

