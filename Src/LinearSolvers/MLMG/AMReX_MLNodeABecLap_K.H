#ifndef AMREX_MLNODEABECLAP_K_H_
#define AMREX_MLNODEABECLAP_K_H_
#include <AMReX_Config.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLNodeABecLap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLNodeABecLap_2D_K.H>
#else
#include <AMReX_MLNodeABecLap_3D_K.H>
#endif

#endif
