#ifndef AMREX_STLBOX_H_
#define AMREX_STLBOX_H_

#include "AMReX_Box.H"
#include "AMReX_RealVect.H"
#include "AMReX_IntVect.H"
#include "AMReX_Box.H"
using std::shared_ptr;

#include "AMReX_STLMesh.H"
#include "AMReX_STLUtil.H"


namespace amrex
{
using namespace STLUtil;
/*
 * This is just a class to hold data that associates
 * a Chombo Box with an STL mesh.  It has, e.g., maps
 * between Chombo cells and the STL vertices that lie
 * within them.  The class also includes some get
 * functions and print functions for its members
 */

  class STLBox
  {
  public:

    STLBox(shared_ptr<STLMesh> a_stlmesh,
           const Box&             a_region,
           const Box          &   a_domain,
           const RealVect&        a_origin,
           const RealVect&        a_dx);

    void SetMeshBox(shared_ptr<STLMesh> a_stlmesh,
                    const Box&             a_region,
                    const Box          &   a_domain,
                    const RealVect&        a_origin,
                    const RealVect&        a_dx);

    shared_ptr<STLMesh> m_msh; // pointer to mesh

    // box that we're exploring
    Box           m_region;
    Box           m_domain;
    RealVect      m_origin;
    RealVect      m_dx;

    // for each cell, a Vector of the vertices and a Vector of the triangles
    // that wholly or partially are contained by the cell
    // vertices and triangles are specified by index (see STLStruct.H for TriInCell)
    stlCellMap m_cellmap;

    // for each Chombo node, whether it is 'inside' or 'outside'
    // according to the normal vectors specified in the STL file and mesh.triangles.normal
    stlNodeMap m_nodemap;

    // for each Chombo edge that intersects the boundary, store the
    // intersection point (since we already compute it, this just makes
    // it faster to retrieve for queries)
    stlEdgeMap m_edgemap;

  };

}
#endif

