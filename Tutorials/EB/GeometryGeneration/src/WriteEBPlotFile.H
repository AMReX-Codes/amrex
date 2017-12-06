#ifndef _WriteEBPlotFile_H_
#define _WriteEBPlotFile_H_

//#include <WritePlotFile.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_EBLevelGrid.H>
#include <CommonCode.H>
#include <iostream>

namespace amrex
{
  ///write a single level plot file with an extra variable called vfrac which holds volume fractions
  void
  WriteEBPlotFile(const std::string               & a_filename,
                  const amrex::MultiFab           & a_mf,
                  const EBLevelGrid               & a_eblg,
                  amrex::Vector<std::string> a_names);

}

#endif

