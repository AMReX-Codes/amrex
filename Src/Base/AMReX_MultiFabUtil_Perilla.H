#ifndef AMREX_MultiFabUtil_PERILLAH_
#define AMREX_MultiFabUtil_PERILLAH_

#ifdef USE_PERILLA
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Amr.H>
#include <RegionGraph.H>

namespace amrex
{

  class RGIter;

    void average_down_push(RGIter& rgi, Amr& amr, MultiFab& S_fine, MultiFab& S_crse, MultiFab& crse_S_fine, RegionGraph* RG_fine, RegionGraph* RG_crse,
			 const Geometry& fine_geom, const Geometry& crse_geom,
			 int scomp, int ncomp, const IntVect& ratio, int f);
    void average_down_pull(RGIter& rgi, MultiFab& S_fine, MultiFab& S_crse, RegionGraph* RG_fine, RegionGraph* RG_crse, const Geometry& fine_geom, const Geometry& crse_geom,
			   int scomp, int ncomp, const IntVect& ratio, int f);
  
    void average_down_push(RGIter& rgi, Amr& amr, MultiFab& S_fine, MultiFab& S_crse, MultiFab& crse_S_fine, RegionGraph* RG_fine, RegionGraph* RG_crse,
			   const Geometry& fine_geom, const Geometry& crse_geom,
			   int scomp, int ncomp, const int ratio, int f);
    void average_down_pull(RGIter& rgi, MultiFab& S_fine, MultiFab& S_crse, RegionGraph* RG_fine, RegionGraph* RG_crse, const Geometry& fine_geom, const Geometry& crse_geom,
			   int scomp, int ncomp, const int ratio, int f);
  
    void average_down_push(RGIter& rgi, Amr& amr, MultiFab& S_fine, MultiFab& S_crse, MultiFab& crse_S_fine, RegionGraph* RG_fine, RegionGraph* RG_crse,
			   int scomp, int ncomp, const IntVect& ratio, int f);
    void average_down_pull(RGIter& rgi, MultiFab& S_fine, MultiFab& S_crse, RegionGraph* RG_fine, RegionGraph* RG_crse,
			   int scomp, int ncomp, const IntVect& ratio, int f);
  
    void average_down_push(RGIter& rgi, Amr& amr, MultiFab& S_fine, MultiFab& S_crse, MultiFab& crse_S_fine, RegionGraph* RG_fine, RegionGraph* RG_crse,
			   int scomp, int ncomp,       int      ratio, int f);  
    void average_down_pull(RGIter& rgi, MultiFab& S_fine, MultiFab& S_crse, RegionGraph* RG_fine, RegionGraph* RG_crse,
			   int scomp, int ncomp,       int      ratio, int f);

}
#endif

#endif
