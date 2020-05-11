#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

void main_main ();

void advance (amrex::MultiFab& phi_old,
              amrex::MultiFab& phi_new,
	      amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& flux,
	      amrex::Real dt,
              const amrex::Geometry& geom);

#endif
