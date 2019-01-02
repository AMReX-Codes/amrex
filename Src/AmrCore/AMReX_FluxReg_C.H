#ifndef AMREX_FLUXREG_C_H_
#define AMREX_FLUXREG_C_H_

#if (AMREX_SPACEDIM == 1)
#include <AMReX_FluxReg_1D_C.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_FluxReg_2D_C.H>
#else
#include <AMReX_FluxReg_3D_C.H>
#endif

#endif
