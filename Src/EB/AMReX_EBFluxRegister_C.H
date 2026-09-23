#ifndef AMREX_EB_FLUXREGISTER_C_H_
#define AMREX_EB_FLUXREGISTER_C_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_YAFluxRegister_K.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {
    //! Cells with smaller volume fraction are treated as covered (mirrors the eb2.small_volfrac default).
#ifdef AMREX_USE_FLOAT
    inline constexpr Real eb_flux_reg_small_volfrac = Real(1.e-5);
#else
    inline constexpr Real eb_flux_reg_small_volfrac = Real(1.e-14);
#endif
}

#if (AMREX_SPACEDIM == 2)
#include <AMReX_EBFluxRegister_2D_C.H>
#elif (AMREX_SPACEDIM == 3)
#include <AMReX_EBFluxRegister_3D_C.H>
#endif

#endif
