#ifndef _WriteEBPlotFile_H____
#define _WriteEBPlotFile_H____

#include <WritePlotFile.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_EBLevelGrid.H>
#include <iostream>

namespace amrex
{
  ///write a single level plot file with an extra variable called vfrac which holds volume fractions
  extern void
  WriteSingleLevelEBPlotFile(const std::string               & a_filename,
                             const amrex::MultiFab           & a_mf,
                             const EBLevelGrid               & a_eblg,
                             amrex::Vector<std::string> a_names = amrex::Vector<std::string>());

}

#endif

