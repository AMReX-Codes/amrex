#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <array>

void main_main ();
void advance (amrex::MultiFab& old_phi, amrex::MultiFab& new_phi,
	      std::array<amrex::MultiFab, AMREX_SPACEDIM>& flux,
	      amrex::Real dt, const amrex::Geometry& geom, int nsteps);

void fill_physbc (amrex::MultiFab& phi, const amrex::Geometry& geom);

#endif
