// --------------------------------------------------------------------
// AppendToPlotFile.H
// --------------------------------------------------------------------
#ifndef BL_APPENDTOPLOTFILE_H
#define BL_APPENDTOPLOTFILE_H

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_AmrData.H>
#include <AMReX_ParallelDescriptor.H>

void
AppendToPlotFile(amrex::AmrData&                        amrData,
                 const amrex::Vector<amrex::MultiFab*>& mfout,
                 const std::string&                     oFile,
                 const amrex::Vector<std::string>&      nnames,
                 const std::string&                     mfBaseName,
                 const std::string&                     NewHeaderName,
                 bool                                   verbose);

#endif
// --------------------------------------------------------------------
