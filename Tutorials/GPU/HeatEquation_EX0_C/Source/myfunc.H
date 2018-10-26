#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void main_main ();

void advance (MultiFab& phi_old,
              MultiFab& phi_new,
	      std::array<MultiFab, AMREX_SPACEDIM>& flux,
	      Real dt,
              Geometry& geom);

/*
template<class L>
__global__ void magic (L f0) { f0(); }

template<class L, class... Lambdas>
__global__ void magic (L f0, Lambdas... fs) { f0(); magic_device(fs...); }

template<class L>
__device__ void magic_device (L f0) { f0(); }

template<class L, class... Lambdas>
__device__ void magic_device (L f0, Lambdas... fs) { f0(); magic_device(fs...); }
*/

#endif
