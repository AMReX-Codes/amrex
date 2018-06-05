#ifndef AMREX_CELLEDGE_H_
#define AMREX_CELLEDGE_H_

#include "AMReX_IntVect.H"

namespace amrex
{
///
/**
   class used in the STL reader.   Basically this defines a face using only IntVects.  
   In the platonically ideal universe, where everything is resolved and nothing is ever multivalued,
   this is semantically the same thing as a FaceIndex.
*/
  class CellEdge
  {
  public:


    ///
    CellEdge(const IntVect& a_node0,
             const IntVect& a_node1);

    ///
    CellEdge(const IntVect& a_node0,
             const int      a_dir);
    ///
    ~CellEdge();

    ///
    void shift(const int a_dir);

    ///
    void shift(const int a_dir,
               const int a_offset);

    /// test CellEdge equality (for unordered_map) with _both_ nodes _and_ the direction
    bool operator==(const CellEdge& a_edge) const;


    int     m_dir;
    IntVect m_node0;
    IntVect m_node1;

  private:

    CellEdge()
      {
      }

    void operator=(const CellEdge& a_inputCellEdge)
      {
      }
  };
}

#endif
