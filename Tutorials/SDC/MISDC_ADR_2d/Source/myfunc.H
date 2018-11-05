#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include "AMReX_SDCstruct.H"
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

// #include <AMReX_LO_BCTYPES.H>

using namespace amrex;

void main_main ();

void SDC_advance (MultiFab& phi_old,
		  MultiFab& phi_new,
		  std::array<MultiFab, AMREX_SPACEDIM>& flux,
		  Real dt,
		  const Geometry& geom, 
		  const Vector<BCRec>& bc,
        	  MLMG&  mlmg,
		  MLABecLaplacian& mlabec,
		  SDCstruct &sdcmats,Real a,Real d,Real r);
void SDC_feval (std::array<MultiFab, AMREX_SPACEDIM>& flux,
		const Geometry& geom, 
		const Vector<BCRec>& bc,
		SDCstruct &sdcmats,
		Real a,Real d,Real r,		
		int sdc_m,int npiece);
void SDC_fcomp(MultiFab& rhs,
	       const Geometry& geom, 
	       const Vector<BCRec>& bc,
	       SDCstruct &sdcmats,
	       MLMG &mlmg,
	       MLABecLaplacian& mlabec,
	       Real dt,Real d,Real r,
	       int sdc_m,int npiece);

#endif
