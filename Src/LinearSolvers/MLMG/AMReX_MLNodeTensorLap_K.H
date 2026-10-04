#ifndef AMREX_MLNODETENSORLAP_K_H_
#define AMREX_MLNODETENSORLAP_K_H_
#include <AMReX_Config.H>

#include <AMReX_MLNodeTensorLaplacian.H>
#include <AMReX_LO_BCTYPES.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLNodeTensorLap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLNodeTensorLap_2D_K.H>
#else
#include <AMReX_MLNodeTensorLap_3D_K.H>
#endif

#endif
