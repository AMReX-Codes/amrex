#ifndef AMReX_extrapolater_K_H_
#define AMReX_extrapolater_K_H_

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_extrapolater_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_extrapolater_2D_K.H>
#else
#include <AMReX_extrapolater_3D_K.H>
#endif

#endif
