#ifndef AMREX_FILL_FC_H_
#define AMREX_FILL_FC_H_
#include <AMReX_Config.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_FilFC_1D_C.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_FilFC_2D_C.H>
#else
#include <AMReX_FilFC_3D_C.H>
#endif
#endif
