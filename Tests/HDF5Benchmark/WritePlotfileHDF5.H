#ifndef WRITEPLOTFILEHDF5_H_
#define WRITEPLOTFILEHDF5_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

#ifdef BL_HDF5
void WriteMultiLevelPlotfileHDF5 (const std::string &plotfilename,
                                  int nlevels,
                                  const amrex::Vector<const amrex::MultiFab*> &mf,
                                  const amrex::Vector<std::string> &varnames,
                                  const amrex::Vector<amrex::Geometry> &geom,
                                  amrex::Real time,
                                  amrex::Real dt,
				  const amrex::Vector<amrex::IntVect> &ref_ratio);
#endif

#endif
