#ifndef CNS_K_H_
#define CNS_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

AMREX_GPU_DEVICE
void cns_initdata (amrex::Box const& bx, amrex::FArrayBox& state, amrex::GeometryData const& geomdata);

AMREX_GPU_DEVICE
amrex::Real cns_estdt (amrex::Box const& bx, amrex::FArrayBox const& state,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx);

AMREX_GPU_DEVICE
void cns_ctoprim (amrex::Box const& bx, amrex::FArrayBox const& u, amrex::FArrayBox& q);

#endif
