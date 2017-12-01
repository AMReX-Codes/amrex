#ifndef _GeomIntersectUtils_H_
#define _GeomIntersectUtils_H_

#include <AMReX_RealVect.H>
#include <AMReX_IntVect.H>
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
    Data (a node) on a set of edges
   */
  typedef std::map<Edge, RealVect> NodeMap;
  typedef NodeMap::iterator NodeMapIt;

  /*
    Connects three nodes from the above map in order to form a 
    triangle.  Typically used as one element in a triangulated 
    surface, representing an isocontour or embedded boundary 
    surface.
   */
  struct Triangle
  {
    Triangle() : p(3) {}
    Real Area();
    const NodeMapIt& operator[] (int n) const { return p[n]; }
    NodeMapIt& operator[] (int n) { return p[n]; }
    Array<NodeMapIt> p;
  };

  /*
    Connects two nodes from the above map in order to form a 
    line segment.  Typically used as one element in a polyline
    representing an isocontour or embedded boundary surface.
   */
  struct Segment
  {
    Segment() {}
    Segment(const NodeMapIt& p1, const NodeMapIt& p2) : p({p1, p2}) {}
    const NodeMapIt& ID_l() const {return p[0];}
    const NodeMapIt& ID_r() const {return p[1];}
    const NodeMapIt& operator[] (int n) const { return p[n]; }
    NodeMapIt& operator[] (int n) { return p[n]; }
    void flip ();
  protected:
    std::array<NodeMapIt,2> p;
  };
}
#endif
