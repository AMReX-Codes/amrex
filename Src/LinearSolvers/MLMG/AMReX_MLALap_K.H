#ifndef AMREX_MLALAP_K_H_
#define AMREX_MLALAP_K_H_

#include <AMReX_FArrayBox.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLALap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLALap_2D_K.H>
#else
#include <AMReX_MLALap_3D_K.H>
#endif

#endif
