#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

const int SDC_NNODES=3;
const int SDC_NPIECES=2;

void main_main ();


void sweep (MultiFab& phi_old,
	    MultiFab& phi_new,
	    std::array<MultiFab, AMREX_SPACEDIM>& flux,
	    Vector<MultiFab>& phi_sdc,
	    Vector<Vector<MultiFab> >& f_sdc,
	    Real dt,
	    const Geometry& geom,
	    Real qnodes [],
	    Real qmat [][SDC_NNODES],
	    Real qmatE [][SDC_NNODES],
	    Real qmatI [][SDC_NNODES]);


#endif
