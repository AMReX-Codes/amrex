#ifndef AMREX_INTERP_BNDRYDATA_K_H_
#define AMREX_INTERP_BNDRYDATA_K_H_

#if (AMREX_SPACEDIM == 1)
#include <AMReX_InterpBndryData_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_InterpBndryData_2D_K.H>
#else
#include <AMReX_InterpBndryData_3D_K.H>
#endif

#endif
