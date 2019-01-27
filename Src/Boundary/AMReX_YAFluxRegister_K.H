#ifndef AMREX_YAFLUXREGISTER_K_H_
#define AMREX_YAFLUXREGISTER_K_H_

namespace {
    // must be same as in AMReX_YAFluxRegister.H
    constexpr int amrex_yafluxreg_crse_cell = 0;
    constexpr int amrex_yafluxreg_crse_fine_boundary_cell = 1;
    constexpr int amrex_yafluxreg_fine_cell = 2;
}

#if (AMREX_SPACEDIM == 1)
#include <AMReX_YAFluxRegister_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_YAFluxRegister_2D_K.H>
#else
#include <AMReX_YAFluxRegister_3D_K.H>
#endif

#endif
