// --------------------------------------------------------------------
// AMReX_BLWritePlotFile.H
// --------------------------------------------------------------------
#ifndef BL_WRITEPLOTFILE_H
#define BL_WRITEPLOTFILE_H

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>

using amrex::Real;

namespace amrex{

void WritePlotfile(const std::string &pfversion,
                   const amrex::Vector<amrex::MultiFab> &data,
                   const Real time,
                   const amrex::Vector<Real> &probLo,
                   const amrex::Vector<Real> &probHi,
                   const amrex::Vector<int>  &refRatio,
                   const amrex::Vector<amrex::Box> &probDomain,
                   const amrex::Vector<amrex::Vector<Real> > &dxLevel,
                   const int coordSys,
                   const std::string &oFile,
                   const amrex::Vector<std::string> &names,
                   const bool verbose,
		   const bool isCartGrid = false,
		   const Real *vfeps = 0,
		   const int  *levelSteps = 0);


// this function makes a 2d plotfile from a 3d one
// the 3d file must have width = 1 in one direction
// compile in 3d

void WritePlotfile2DFrom3D(const std::string &pfversion,
                           const amrex::Vector<amrex::MultiFab> &data,
                           const Real time,
                           const amrex::Vector<Real> &probLo,
                           const amrex::Vector<Real> &probHi,
                           const amrex::Vector<int>  &refRatio,
                           const amrex::Vector<amrex::Box> &probDomain,
                           const amrex::Vector<amrex::Vector<Real> > &dxLevel,
                           const int coordSys,
                           const std::string &oFile,
                           const amrex::Vector<std::string> &names,
                           const bool verbose,
		           const bool isCartGrid = false,
		           const Real *vfeps = 0,
		           const int  *levelSteps = 0);

} //namespace amrex
#endif
// --------------------------------------------------------------------
// --------------------------------------------------------------------
