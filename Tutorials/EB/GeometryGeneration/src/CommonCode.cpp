#include <WriteEBPlotFile.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxIterator.H>
#include "AMReX_ParmParse.H"
namespace amrex
{
/****/
  string convertInt(int number)
  {
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
  }
/****/
  string convertReal(Real number)
  {
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
  }


/****/
  void
  getFinestDomain(Box&       a_domain,
                  RealVect&      a_dx)
  {
    ParmParse pp;
    std::vector<int> n_cell(SpaceDim);
    pp.getarr("n_cell",n_cell,0,SpaceDim);

    BL_ASSERT(n_cell.size() == SpaceDim);
    IntVect lo = IntVect::TheZeroVector();
    IntVect hi;
    for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
      {
        amrex::Print() << " bogus number of cells input = " << n_cell[ivec];
        amrex::Error();
      }
      hi[ivec] = n_cell[ivec] - 1;
    }

    a_domain = Box(lo, hi);

    Real prob_hi;
    pp.get("domain_length",prob_hi);
    for (int idir=0;idir<SpaceDim;idir++)
    {
      a_dx[idir] = prob_hi/n_cell[idir];
    }
  }

/****/
  void 
  setDataToSomething(MultiFab         & a_data,
                     const EBLevelGrid& a_eblg)
  {
    a_data.setVal(-1.0);
    Real pi = 4.*atan(1.0);
    for(MFIter mfi(a_data); mfi.isValid(); ++mfi)
    {
      EBISBox ebisBox = a_eblg.getEBISL()[mfi];
      EBGraph ebgraph = ebisBox.getEBGraph();
      Box valid              = a_eblg.getDBL()  [mfi];
      IntVectSet ivs(valid);
      Real dx = ebisBox.getDomain().size()[0];
      dx = 1.0/dx;
      for(VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit)
      {
        Real x = vofit().gridIndex()[0]*dx;
        Real y = vofit().gridIndex()[1]*dx;
        Real val  = 10*sin(pi*x)*sin(pi*y);
        a_data[mfi](vofit().gridIndex(), 0) = val;
        ++val;
      }
    }
  }
}
/****/
