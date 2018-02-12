#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
// #include <AMReX_LO_BCTYPES.H>

using namespace amrex;

void main_main ();

void advance (MultiFab& phi_old,
              MultiFab& phi_new,
	      Real dt,
              const Geometry& geom, 
	      const BoxArray& grids,
	      const DistributionMapping& dmap, 
              const Vector<BCRec>& bc);

#endif
