#ifndef AMREX_MLNODELAP_K_H_
#define AMREX_MLNODELAP_K_H_

#include <AMReX_FArrayBox.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBCellFlag.H>
#endif

#include <AMReX_MLNodeLap_ND_K.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLNodeLap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLNodeLap_2D_K.H>
#else
#include <AMReX_MLNodeLap_3D_K.H>
#endif

#endif
