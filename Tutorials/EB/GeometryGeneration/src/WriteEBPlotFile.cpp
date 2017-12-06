#include <WriteEBPlotFile.H>
#include <WritePlotFile.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_EBLevelGrid.H>
#include "CommonCode.H"
#include "AMReX_MFIter.H"


namespace amrex
{
  string convertIntA(int number)
  {
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
  }


///write a single level plot file with an extra variable called vfrac which holds volume fractions
  void
  WriteSingleLevelEBPlotFile (const std::string               & a_filename,
                              const amrex::MultiFab           & a_mf,
                              const EBLevelGrid               & a_eblg,
                              amrex::Vector<std::string>        a_names)
  {
    Vector<std::string> allnames;
    if(a_names.size() == a_mf.nComp())
    {
      allnames = a_names;
    }
    else
    {
      allnames.resize(a_mf.nComp());
      for(int icomp = 0; icomp < a_mf.nComp(); icomp++)
      {
        allnames[icomp] = string("var") + convertIntA(icomp);
      }
    }
    allnames.push_back(string("vfrac"));

    DistributionMapping dm = a_mf.DistributionMap();
    BoxArray            ba = a_mf.boxArray();
    MultiFab allmf(ba, dm, a_mf.nComp()+1, 0);
    //copy data over
    allmf.copy(a_mf, 0, 0, a_mf.nComp());

    //fill in volume fraction data
    for(MFIter mfi(allmf); mfi.isValid(); ++mfi)
    {
      const EBISBox& ebisBox = a_eblg.getEBISL()[mfi];
      const Box& box = ba[mfi];
      for(BoxIterator bit(box); bit.ok(); ++bit)
      {
        VolIndex vof(bit(), 0); //assuming no multivalued here.
        Real vfrac = ebisBox.volFrac(vof);
        allmf[mfi](bit(), a_mf.nComp()) = vfrac;
      }
    }

    //output to plot file.   If I could figure out how to manufacture an AmrData, I would do multilevel as well.
    Geometry geom(a_eblg.getDomain());
    Real bgVal = 1.0e30;
    IntVect refrat = 2*IntVect::Unit;
    writePlotFile(a_filename.c_str(), allmf, geom, refrat, bgVal, allnames);
  }


}
