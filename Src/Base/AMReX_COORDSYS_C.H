#ifndef AMREX_COORDSYS_F_H_
#define AMREX_COORDSYS_F_H_

#include <AMReX_CudaQualifiers.H>
#include <AMReX_Array.H>
#include <AMReX_FArrayBox.H>

namespace amrex {

    AMREX_CUDA_HOST_DEVICE
    void amrex_setvol (Box const& bx, FArrayBox& vol,
                       CudaArray<Real,AMREX_SPACEDIM> const& offset,
                       CudaArray<Real,AMREX_SPACEDIM> const& dx, const int coord);

    AMREX_CUDA_HOST_DEVICE
    void amrex_setarea (Box const& bx, FArrayBox& area,
                        CudaArray<Real,AMREX_SPACEDIM> const& offset,
                        CudaArray<Real,AMREX_SPACEDIM> const& dx, const int dir, const int coord);

    AMREX_CUDA_HOST_DEVICE
    void amrex_setdloga (Box const& bx, FArrayBox& dloga,
                         CudaArray<Real,AMREX_SPACEDIM> const& offset,
                         CudaArray<Real,AMREX_SPACEDIM> const& dx, const int dir, const int coord);
}

#endif
