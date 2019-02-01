#ifndef AMREX_MG_K_H_
#define AMREX_MG_K_H_

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MG_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MG_2D_K.H>
#else
#include <AMReX_MG_3D_K.H>
#endif

#endif
