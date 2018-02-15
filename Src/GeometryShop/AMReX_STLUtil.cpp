#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "STLUtil.H"

#include "NamespaceHeader.H"

/*
 * Printing functions
 */

namespace STLUtil
{

  void PMap(const CellMap& m)
  {
    // print out map
    CellMap::const_iterator it;
    pout() << "Cell Map has " << m.size() << " cells\n";
    for ( it = m.begin(); it != m.end(); it++)
    {
      pout() << "cell "; PIV(it->first); pout() << "\n";
      pout() << " verts: "; PVec(it->second.vertices); pout() << "\n";
      pout() << " tris : "; PVec(it->second.triangles); pout() << "\n";
    }
  }

  void PMap(const pair<IntVect, TriInCell>& p)
  {
    // print out map
    pout() << "cell "; PIV(p.first); pout() << "\n";
    pout() << " verts: "; PVec(p.second.vertices); pout() << "\n";
    pout() << " tris : "; PVec(p.second.triangles); pout() << "\n";
  }

  void PMap(const NodeMap& m)
  {
    // print out map
    NodeMap::const_iterator it;
    pout() << "Node map has " << m.size() << " nodes\n";
    for (it=m.begin(); it!=m.end(); it++)
    {
      pout() << "node "; PIV(it->first); pout() << ": " << it->second << "\n";
    }
  }

  void PIV(const IntVect& iv)
  {
    pout() << "(";
    for (int i=0; i<(SpaceDim-1); i++)
      pout() << iv[i] << ",";
    pout() << iv[SpaceDim-1] << ")"; // last element without comma
  }

  void PRV(const RealVect& iv)
  {
    pout() << "(";
    for (int i=0; i<(SpaceDim-1); i++)
      pout() << iv[i] << ",";
    pout() << iv[SpaceDim-1] << ")"; // last element without comma
  }

  void PVec(const Vector<int>& v)
  {
    if (v.size()<1)
      return;
    for (int i=0; i<(v.size()-1); i++)
      pout() << v[i] << ",";
    pout() << v[v.size()-1]; // last element without comma
  }

  void PVec(const Vector<IntVect>& v)
  {
    if (v.size()<1)
      return;
    for (int i=0; i<v.size(); i++)
    {
      pout() << "\n   " << i << ": "; PIV(v[i]);
    }
  }

  void PVec(const Vector<RealVect>& v)
  {
    if (v.size()<1)
      return;
    for (int i=0; i<v.size(); i++)
    {
      pout() << "\n   " << i << ": "; PRV(v[i]);
    }
  }

  void PVec(const Vector< Vector<IntVect> >& v)
  {
    if (v.size()<1)
      return;
    for (int i=0; i<v.size(); i++)
    {
      pout() << "\n  " << i << ": ";
      PVec(v[i]);
    }
  }

  void PVec(const Vector< Vector<int> >& v)
  {
    if (v.size()<1)
      return;
    for (int i=0; i<v.size(); i++)
    {
      pout() << "\n  " << i << ": ";
      PVec(v[i]);
    }
  }

  // convert IntVect to it's physical location
  RealVect IVToRV(const IntVect& iv,
                  const RealVect& a_origin,
                  const RealVect& a_dx)
  {
    RealVect out = a_origin;
    for (int idir=0; idir<SpaceDim; idir++)
    {
      out[idir] += ((Real) iv[idir]) * a_dx[idir];
    }
    return out;
  }

  // get signum of the values in an IntVect
  IntVect RVSign(const RealVect& rv)
  {
    IntVect rval;
    for (int i=0; i<SpaceDim; i++)
      rval[i] = (0 < rv[i]) - (rv[i] < 0);
    return rval;
  }

}


#include "NamespaceFooter.H"
