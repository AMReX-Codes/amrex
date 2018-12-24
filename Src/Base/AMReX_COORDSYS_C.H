#ifndef AMREX_COORDSYS_C_H_
#define AMREX_COORDSYS_C_H_

#if (AMREX_SPACEDIM == 1)
#include <AMReX_COORDSYS_1D_C.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_COORDSYS_2D_C.H>
#else
#include <AMReX_COORDSYS_3D_C.H>
#endif

#endif
