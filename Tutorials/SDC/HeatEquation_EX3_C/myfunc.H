#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include "AMReX_SDCstruct.H"
// #include <AMReX_LO_BCTYPES.H>

using namespace amrex;

void main_main ();

void sweep (MultiFab& phi_old,
	    MultiFab& phi_new,
	    std::array<MultiFab, AMREX_SPACEDIM>& flux,
	    Vector<MultiFab>& phi_sdc,
	    Vector<MultiFab>& res_sdc,
	    //	    Vector<Vector<MultiFab> >& f_sdc,
	    Vector<MultiFab>& f_sdc,	    
	    Real dt,
	    const Geometry& geom, 
	    const BoxArray& grids,
	    const DistributionMapping& dmap, 
	    const Vector<BCRec>& bc,
            SDCstruct &sdcmats);

#endif
