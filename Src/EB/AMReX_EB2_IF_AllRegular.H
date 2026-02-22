#ifndef AMREX_EB2_IF_ALL_REGULAR_H_
#define AMREX_EB2_IF_ALL_REGULAR_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

/**
 * \file AMReX_EB2_IF_AllRegular.H
 *
 * Implicit-function helper that marks the entire domain as regular (no EB).
 */

namespace amrex::EB2 {

// For all implicit functions, >0: body; =0: boundary; <0: fluid

//! Trivial implicit function that always returns the "fluid" sign.
class AllRegularIF
    : public GPUable
{
public:
    //! Always returns a negative value so EB2 treats every cell as regular.
    constexpr Real operator() (const RealArray&) const noexcept { return -1.0; }

    //! Overload for AMREX_D_DECL signature to avoid temporary RealArray creation.
    constexpr Real operator() (AMREX_D_DECL(Real, Real, Real)) const noexcept { return -1.0; }
};

}

#endif
