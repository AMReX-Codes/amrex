#ifndef AMREX_MULITIFAB_UTIL_C_H_
#define AMREX_MULITIFAB_UTIL_C_H_

#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE
void amrex_int_to_real (Box const& bx, FArrayBox& rfab, IArrayBox const& ifab,
                        int rcomp, int icomp, int ncomp);


AMREX_GPU_HOST_DEVICE
void amrex_fill_slice_interp (Box const& bx, FArrayBox & slicefab, FArrayBox const& fullfab,
                              int scomp, int fcomp, int ncomp,
                              int dir, Real coord, GeometryData const& gd);

AMREX_GPU_HOST_DEVICE
void amrex_avg_nd_to_cc (Box const& bx, FArrayBox& ccfab, FArrayBox const& ndfab,
                         int cccomp, int ndcomp, int ncomp);

AMREX_GPU_HOST_DEVICE
void amrex_avg_eg_to_cc (Box const& bx, FArrayBox& ccfab,
                         AMREX_D_DECL(FArrayBox const& exfab,
                                      FArrayBox const& eyfab,
                                      FArrayBox const& ezfab));

AMREX_GPU_HOST_DEVICE
void amrex_avg_fc_to_cc (Box const& bx, FArrayBox& ccfab,
                         AMREX_D_DECL(FArrayBox const& exfab,
                                      FArrayBox const& eyfab,
                                      FArrayBox const& ezfab),
                         GeometryData const& gd);

}

#endif

