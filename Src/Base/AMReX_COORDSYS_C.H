#ifndef AMREX_COORDSYS_F_H_
#define AMREX_COORDSYS_F_H_

#include <AMReX_Gpu.H>
#include <AMReX_Array.H>
#include <AMReX_FArrayBox.H>

namespace amrex {

    AMREX_GPU_HOST_DEVICE
    void amrex_setvol (Box const& bx, FArrayBox& vol,
                       GpuArray<Real,AMREX_SPACEDIM> const& offset,
                       GpuArray<Real,AMREX_SPACEDIM> const& dx, const int coord);

    AMREX_GPU_HOST_DEVICE
    void amrex_setarea (Box const& bx, FArrayBox& area,
                        GpuArray<Real,AMREX_SPACEDIM> const& offset,
                        GpuArray<Real,AMREX_SPACEDIM> const& dx, const int dir, const int coord);

    AMREX_GPU_HOST_DEVICE
    void amrex_setdloga (Box const& bx, FArrayBox& dloga,
                         GpuArray<Real,AMREX_SPACEDIM> const& offset,
                         GpuArray<Real,AMREX_SPACEDIM> const& dx, const int dir, const int coord);
}

#endif
