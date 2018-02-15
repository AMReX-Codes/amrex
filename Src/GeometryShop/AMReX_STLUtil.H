#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _STLUTIL_H_
#define _STLUTIL_H_

/**********************************************************************/
// to switch between ordered/sorted maps and unordered hash tables,
// define this to be anything, and make sure to compile with -std=c++0x
/**********************************************************************/

// #define STL_UNORDERED_MAP 1

#ifndef STL_UNORDERED_MAP
#include <map>
#else
#include <unordered_map>
#endif

#include "RealVect.H"
#include "IntVect.H"
#include "MayDay.H"
#include "CellEdge.H"

using namespace std;

#include "NamespaceHeader.H"

namespace STLUtil
{
  // returns true if a[0]<b[0];
  // else if a[0]==b[0], returns true if a[1]<b[1], etc.
  // allows sorting by each component of an IntVect sequentially
  // "Strict Weak Ordering" of IntVects
  //
  // Yes, this is just like IntVect::lexLT except that
  // you can specify a direction to compare first
  // i.e. if you want to sort by y, then z and x,
  // you can create a comparator via
  // IVCompareSWO myComparator(1);
  struct IVCompareSWO
  {
    IVCompareSWO(int a_dir) {
      m_dir = a_dir;
      if (m_dir<0 || m_dir>(SpaceDim-1))
        MayDay::Abort("IVCompareSWO: cannot compare IntVects along a direction that doesn't exist");
    }

    IVCompareSWO() {
      m_dir = 0; // default to x,y,z order
    }

    inline bool operator()(const IntVect& a, const IntVect& b) const
    {
      if (a[m_dir%SpaceDim] < b[m_dir%SpaceDim])
        return true;
      else if (SpaceDim>1 && \
          a[ m_dir   %SpaceDim] == b[ m_dir   %SpaceDim] && \
          a[(m_dir+1)%SpaceDim] <  b[(m_dir+1)%SpaceDim])
        return true;
      else if (SpaceDim>2 && \
          a[ m_dir   %SpaceDim] == b[ m_dir   %SpaceDim] && \
          a[(m_dir+1)%SpaceDim] == b[(m_dir+1)%SpaceDim] && \
          a[(m_dir+2)%SpaceDim] <  b[(m_dir+2)%SpaceDim])
        return true;
      else
        return false;
    }

  private:
    int m_dir; // compare along this direction first
  };

  // a functor to compare edges for map/sort purposes
  // first sort based on node0, then compare the direction of the edge
  struct EdgeCompareSWO
  {
    inline bool operator()(const CellEdge& a,const CellEdge& b) const
    {
      if (a.m_node0==b.m_node0)
        return (a.m_dir<b.m_dir); // compare directions, first x, then y,z

      IVCompareSWO comparator;
      return comparator(a.m_node0,b.m_node0);
    }
  };

  // a functor to compare pairs of real & int based on the real
  struct RealIntCompare
  {
    inline bool operator()(const pair<Real,int>& a,const pair<Real,int>& b)
    {
      return a.first<b.first;
    }
  };

#ifdef STL_UNORDERED_MAP
  // for combining hashes, stolen from Boost
  template <class T>
  inline void hash_combine(std::size_t& seed, const T& v)
  {
      std::hash<T> hasher;
      seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
  }

  // now make a hash function for intvects
  struct IVHash
    : std::unary_function<IntVect, std::size_t>
  {
    inline std::size_t operator()(const IntVect& a) const
    {
      std::size_t seed = 0;
      for (int i=0; i<SpaceDim; i++)
        hash_combine(seed,a[i]);
      return seed;
    }
  };

  // use this to make a hash f'n for a CellEdge, using just m_node0 and m_dir
  // yeah, this could be defined as CellEdge::hash_value(), but whatever
  struct CellEdgeHash
    : std::unary_function<CellEdge, std::size_t>
  {
  inline std::size_t operator()(const CellEdge& a) const
    {
      IVHash node0hasher;
      std::size_t seed = node0hasher(a.m_node0);
      hash_combine(seed,a.m_dir);
      return seed;
    }
  };
#endif

  typedef struct {
    // single cell has vertices and triangles in it
    Vector<int> vertices;
    Vector<int> triangles;
  } TriInCell;

  // map definitions, depending on whether or not we can use unordered maps
#ifndef STL_UNORDERED_MAP
  typedef map<IntVect, TriInCell, IVCompareSWO> CellMap;
  typedef map<IntVect, bool, IVCompareSWO>      NodeMap;
  typedef map<CellEdge, RealVect, EdgeCompareSWO> EdgeMap;
#else
  typedef unordered_map<IntVect, TriInCell, IVHash> CellMap;
  typedef unordered_map<IntVect, bool, IVHash>      NodeMap;
  typedef unordered_map<CellEdge, RealVect, CellEdgeHash> EdgeMap;
#endif

  // printing functions
  void PMap(const CellMap& m); // print a cell map
  void PMap(const NodeMap& m); // print a node map
  void PMap(const pair<IntVect, TriInCell>& p); // print a cell map entry
  void PIV(const IntVect& iv); // print an intvect
  void PRV(const RealVect& iv); // print a realvect
  void PVec(const Vector<int>& v); // print a vector of things
  void PVec(const Vector<IntVect>& v);
  void PVec(const Vector<RealVect>& v);
  void PVec(const Vector< Vector<IntVect> >& v);
  void PVec(const Vector< Vector<int> >& v);

  // convert IntVect to it's physical location
  RealVect IVToRV(const IntVect& iv,
                  const RealVect& a_origin,
                  const RealVect& a_dx);

  // get signum of the values in an IntVect
  IntVect RVSign(const RealVect& rv);

  typedef CellMap::iterator CellMapIt;
  typedef NodeMap::iterator NodeMapIt;
  typedef EdgeMap::iterator EdgeMapIt;

}

#include "NamespaceFooter.H"
#endif

