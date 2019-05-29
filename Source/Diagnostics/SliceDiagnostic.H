#ifndef WARPX_SliceDiagnostic_H_
#define WARPX_SliceDiagnostic_H_

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Geometry.H>

#include <WarpX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Vector.H>
#include <AMReX_BLassert.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFabUtil_C.H>



std::unique_ptr<amrex::MultiFab> CreateSlice( const amrex::MultiFab& mf, 
               const amrex::Vector<amrex::Geometry> &dom_geom, 
               amrex::RealBox &slice_realbox, 
               amrex::IntVect &slice_cr_ratio );

void CheckSliceInput( const amrex::RealBox real_box, 
     amrex::RealBox &slice_cc_nd_box, amrex::RealBox &slice_realbox, 
     amrex::IntVect &slice_cr_ratio, amrex::Vector<amrex::Geometry> dom_geom,
     amrex::IntVect const SliceType, amrex::IntVect &slice_lo, 
     amrex::IntVect &slice_hi, amrex::IntVect &interp_lo);

void InterpolateSliceValues( amrex::MultiFab& smf, 
     amrex::IntVect interp_lo, amrex::RealBox slice_realbox, 
     amrex::Vector<amrex::Geometry> geom, int ncomp, int nghost, 
     amrex::IntVect slice_lo, amrex::IntVect slice_hi, 
     amrex::IntVect SliceType, const amrex::RealBox real_box);

void InterpolateLo( const amrex::Box& bx, amrex::FArrayBox &fabox, 
     amrex::IntVect slice_lo, amrex::Vector<amrex::Geometry> geom, 
     int idir, amrex::IntVect IndType, amrex::RealBox slice_realbox, 
     int srccomp, int ncomp, int nghost, const amrex::RealBox real_box);

#endif

