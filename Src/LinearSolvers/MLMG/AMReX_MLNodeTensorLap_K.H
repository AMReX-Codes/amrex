#ifndef AMREX_MLNODETENSORLAP_K_H_
#define AMREX_MLNODETENSORLAP_K_H_

#include <AMReX_MLNodeTensorLaplacian.H>
#include <AMReX_LO_BCTYPES.H>

namespace amrex {

namespace {

    constexpr int ixx = 0;
    constexpr int ixy = 1;
#if (AMREX_SPACEDIM == 2)
    constexpr int iyy = 2;
#else
    constexpr int ixz = 2;
    constexpr int iyy = 3;
    constexpr int iyz = 4;
    constexpr int izz = 5;
#endif
    constexpr int nelems = MLNodeTensorLaplacian::nelems;
}

}

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLNodeTensorLap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLNodeTensorLap_2D_K.H>
#else
#include <AMReX_MLNodeTensorLap_3D_K.H>
#endif

#endif
