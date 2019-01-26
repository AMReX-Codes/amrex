#ifndef AMREX_YAFLUXREGISTER_K_H_
#define AMREX_YAFLUXREGISTER_K_H_

#if (AMREX_SPACEDIM == 1)
#include <AMReX_YAFluxRegister_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_YAFluxRegister_2D_K.H>
#else
#include <AMReX_YAFluxRegister_3D_K.H>
#endif

#endif
