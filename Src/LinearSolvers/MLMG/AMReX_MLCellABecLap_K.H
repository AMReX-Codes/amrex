#ifndef AMREX_MLCELLABECLAP_K_H_
#define AMREX_MLCELLABECLAP_K_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLCellABecLap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLCellABecLap_2D_K.H>
#else
#include <AMReX_MLCellABecLap_3D_K.H>
#endif

#endif
