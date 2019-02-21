#ifndef AMREX_MLPOISSON_K_H_
#define AMREX_MLPOISSON_K_H_

#include <AMReX_FArrayBox.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLPoisson_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLPoisson_2D_K.H>
#else
#include <AMReX_MLPoisson_3D_K.H>
#endif

#endif
