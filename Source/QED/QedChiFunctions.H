#ifndef WARPX_amrex_qed_chi_functions_h_
#define WARPX_amrex_qed_chi_functions_h_

#include "QedWrapperCommons.H"

#define PXRMP_CORE_ONLY
#include <chi_functions.hpp>

namespace QedUtils{
    //Thin wrappers around functions used to calculate chi parameters

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real chi_photon(
        amrex::Real px, amrex::Real py, amrex::Real pz,
        amrex::Real ex, amrex::Real ey, amrex::Real ez,
        amrex::Real bx, amrex::Real by, amrex::Real bz)
    {
        const amrex::Real dummy_lambda = 1.0;
        return picsar::multi_physics::chi_photon(
            px, py, pz, ex, ey, ez, bx, by, bz, dummy_lambda);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real chi_lepton(
        amrex::Real px, amrex::Real py, amrex::Real pz,
        amrex::Real ex, amrex::Real ey, amrex::Real ez,
        amrex::Real bx, amrex::Real by, amrex::Real bz)
    {
        const amrex::Real dummy_lambda = 1.0;
        return picsar::multi_physics::chi_lepton(
            px, py, pz, ex, ey, ez, bx, by, bz, dummy_lambda);
    }
    //_________
};

#endif //WARPX_amrex_qed_chi_functions_h_