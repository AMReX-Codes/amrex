#ifndef AMREX_MLNODELAP_K_H_
#define AMREX_MLNODELAP_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_LO_BCTYPES.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBCellFlag.H>
#endif

namespace amrex {
namespace {

    constexpr int crse_cell = 0;
    constexpr int fine_cell = 1;
    constexpr int crse_node = 0;
    constexpr int crse_fine_node = 1;
    constexpr int fine_node = 2;
    constexpr double eps = 1.e-100;

}}

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLNodeLap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLNodeLap_2D_K.H>
#else
#include <AMReX_MLNodeLap_3D_K.H>
#endif

#endif
