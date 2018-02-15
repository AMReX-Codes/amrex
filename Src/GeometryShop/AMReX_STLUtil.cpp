
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



#include "AMReX_STLUtil.H"
#include "AMReX_Print.H"


/*
 * Printing functions
 */

namespace amrex
{
namespace STLUtil
{

  void PMap(const stlCellMap& m)
  {
    // print out map
    stlCellMap::const_iterator it;
    amrex::Print() << "Cell Map has " << m.size() << " cells\n";
    for ( it = m.begin(); it != m.end(); it++)
    {
      amrex::Print() << "cell "; PIV(it->first); amrex::Print() << "\n";
      amrex::Print() << " verts: "; PVec(it->second.vertices); amrex::Print() << "\n";
      amrex::Print() << " tris : "; PVec(it->second.triangles); amrex::Print() << "\n";
    }
  }

  void PMap(const pair<IntVect, TriInCell>& p)
  {
    // print out map
    amrex::Print() << "cell "; PIV(p.first); amrex::Print() << "\n";
    amrex::Print() << " verts: "; PVec(p.second.vertices); amrex::Print() << "\n";
    amrex::Print() << " tris : "; PVec(p.second.triangles); amrex::Print() << "\n";
  }

  void PMap(const stlNodeMap& m)
  {
    // print out map
    stlNodeMap::const_iterator it;
    amrex::Print() << "Node map has " << m.size() << " nodes\n";
    for (it=m.begin(); it!=m.end(); it++)
    {
      amrex::Print() << "node "; PIV(it->first); amrex::Print() << ": " << it->second << "\n";
    }
  }

  void PIV(const IntVect& iv)
  {
    amrex::Print() << "(";
    for (int i=0; i<(SpaceDim-1); i++)
      amrex::Print() << iv[i] << ",";
    amrex::Print() << iv[SpaceDim-1] << ")"; // last element without comma
  }

  void PRV(const RealVect& iv)
  {
    amrex::Print() << "(";
    for (int i=0; i<(SpaceDim-1); i++)
      amrex::Print() << iv[i] << ",";
    amrex::Print() << iv[SpaceDim-1] << ")"; // last element without comma
  }

  void PVec(const Vector<int>& v)
  {
    if (v.size()<1)
      return;
    for (int i=0; i<(v.size()-1); i++)
      amrex::Print() << v[i] << ",";
    amrex::Print() << v[v.size()-1]; // last element without comma
  }

  void PVec(const Vector<IntVect>& v)
  {
    if (v.size()<1)
      return;
    for (int i=0; i<v.size(); i++)
    {
      amrex::Print() << "\n   " << i << ": "; PIV(v[i]);
    }
  }

  void PVec(const Vector<RealVect>& v)
  {
    if (v.size()<1)
      return;
    for (int i=0; i<v.size(); i++)
    {
      amrex::Print() << "\n   " << i << ": "; PRV(v[i]);
    }
  }

  void PVec(const Vector< Vector<IntVect> >& v)
  {
    if (v.size()<1)
      return;
    for (int i=0; i<v.size(); i++)
    {
      amrex::Print() << "\n  " << i << ": ";
      PVec(v[i]);
    }
  }

  void PVec(const Vector< Vector<int> >& v)
  {
    if (v.size()<1)
      return;
    for (int i=0; i<v.size(); i++)
    {
      amrex::Print()     << "\n  " << i << ": ";
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
}


