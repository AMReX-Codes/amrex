#ifndef _WriteEBPlotFile_H_
#define _WriteEBPlotFile_H_

#include "WritePlotFile.H"
#include "AMReX_MultiFab.H"
#include "AMReX_Geometry.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_REAL.H"
#include <iostream>

using std::string;

namespace amrex
{
///for fancy filenames
  extern string convertInt(int number);

///for fancy filenames
  extern string convertReal(Real number);


///read stuff from parmparse
  extern  void
  getFinestDomain(Box&       a_domain,
                  RealVect&      a_dx);


//puts positive data at all valid points, sets covered values to -1
  extern void 
  setDataToSomething(MultiFab         & a_mf,
                     const EBLevelGrid& a_eblg);
}
#endif

