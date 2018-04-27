#ifndef AMREX_GeomIntersectUtils_H_
#define AMREX_GeomIntersectUtils_H_

#include <AMReX_IntVect.H>
#include <AMReX_Print.H>
#include <list>
#include <map>

namespace amrex
{
  /*
    A data structure with which to associate data with the
    edges of cells.  An edge connects two indices in an index
    space, typically separated by a unit vector in one of the 
    available directions.
   */
  struct Edge
  {
    Edge(const IntVect& lhs, const IntVect& rhs);

    bool operator== (const Edge& rhs) const;
    bool operator< (const Edge& rhs) const;

    IntVect IV_l,IV_r;
    size_t ID;
  };

  std::ostream& operator<<(std::ostream& os, const Edge& e);

  /*
    A data set on edges
   */
  typedef std::pair<Edge, amrex::Vector<amrex::Real>> Node;
  typedef std::map<Edge, amrex::Vector<amrex::Real>> NodeMap;
  typedef NodeMap::iterator NodeMapIt;



  struct ContourFragment
  {
    ContourFragment()  {}
    ContourFragment(size_t len) : m_pts(len) {}
    const NodeMapIt& operator[] (int n) const { return m_pts[n]; }
    NodeMapIt& operator[] (int n) { return m_pts[n]; }
    size_t size() const {return m_pts.size();}

  protected:
    Vector<NodeMapIt> m_pts;
  };

  /*
    Connects three nodes from the above map in order to form a 
    triangle.  Typically used as one element in a triangulated 
    surface, representing an isocontour or embedded boundary 
    surface.
   */
  struct Triangle
    : public ContourFragment
  {
    Triangle() : ContourFragment(3) {}
    Triangle(const NodeMapIt& p0, const NodeMapIt& p1, const NodeMapIt& p2) : ContourFragment(3)
    {
      (*this)[0] = p0;
      (*this)[1] = p1;
      (*this)[2] = p2;
    }
  };

  /*
    Connects two nodes from the above map in order to form a 
    line segment.  Typically used as one element in a polyline
    representing an isocontour or embedded boundary surface.
   */
  struct Segment
    : public ContourFragment
  {
    Segment() : ContourFragment(2) {}
    Segment(const NodeMapIt& p0, const NodeMapIt& p1) : ContourFragment(2)
    {
      (*this)[0] = p0;
      (*this)[1] = p1;
    }
  };

  struct GeomIntersectUtils
  {
    static int MarchingCubesEdgeTable[256];
    static int MarchingCubesTriTable[256][16];
  };
}
#endif
