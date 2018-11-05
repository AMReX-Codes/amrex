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
void cns_compute_temperature (amrex::Box const& bx, amrex::FArrayBox& state);

AMREX_GPU_DEVICE
void cns_ctoprim (amrex::Box const& bx, amrex::FArrayBox const& u, amrex::FArrayBox& q);

AMREX_GPU_DEVICE
void cns_flux_to_dudt (amrex::Box const& bx, amrex::FArrayBox& dudt,
                       AMREX_D_DECL(amrex::FArrayBox const& fx,
                                    amrex::FArrayBox const& fy,
                                    amrex::FArrayBox const& fz),
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv);

AMREX_GPU_DEVICE
void cns_slope_x (amrex::Box const& bx, amrex::FArrayBox& slope, amrex::FArrayBox const& q);

AMREX_GPU_DEVICE
void cns_slope_y (amrex::Box const& bx, amrex::FArrayBox& slope, amrex::FArrayBox const& q);

AMREX_GPU_DEVICE
void cns_slope_z (amrex::Box const& bx, amrex::FArrayBox& slope, amrex::FArrayBox const& q);

AMREX_GPU_DEVICE
void cns_riemann_x (amrex::Box const& bx, amrex::FArrayBox& flux,
                    amrex::FArrayBox const& slope, amrex::FArrayBox const& q);

AMREX_GPU_DEVICE
void cns_riemann_y (amrex::Box const& bx, amrex::FArrayBox& flux,
                    amrex::FArrayBox const& slope, amrex::FArrayBox const& q);

AMREX_GPU_DEVICE
void cns_riemann_z (amrex::Box const& bx, amrex::FArrayBox& flux,
                    amrex::FArrayBox const& slope, amrex::FArrayBox const& q);

#endif
