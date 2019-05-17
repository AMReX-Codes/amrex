#ifndef WARPX_SliceDiagnostic_H_
#define WARPX_SliceDiagnostic_H_

#include <vector>

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

using namespace amrex;



std::unique_ptr<MultiFab> CreateSlice( const amrex::MultiFab& mf, 
                          const amrex::Vector<Geometry> &dom_geom, 
                          amrex::RealBox &slice_realbox, 
                          amrex::IntVect &slice_cr_ratio );


void CheckSliceInput( const RealBox real_box, RealBox &slice_cc_nd_box, 
     RealBox &slice_realbox, IntVect &cr_ratio, IntVect slice_cr_ratio, 
     Vector<Geometry> dom_geom, IntVect const SliceType, IntVect &slice_lo, 
     IntVect &slice_hi, IntVect &interp_lo);

void InterpolateSliceValues( MultiFab& smf, IntVect interp_lo, 
     RealBox slice_realbox, Vector<Geometry> geom, int ncomp, 
     int nghost, IntVect slice_lo, IntVect slice_hi, 
     IntVect SliceType, const RealBox real_box);

void InterpolateLo( const Box& bx, FArrayBox &fabox, IntVect slice_lo,
     Vector<Geometry> geom, int idir, IntVect IndType, RealBox slice_realbox, 
     int srccomp, int ncomp, int nghost, const RealBox real_box);

#endif

