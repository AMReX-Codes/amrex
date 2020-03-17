/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_amrex_qed_chi_functions_h_
#define WARPX_amrex_qed_chi_functions_h_

/**
 * This header contains wrappers around functions provided by
 * the PICSAR QED library to calculate the 'chi' parameter
 * for photons or electrons and positrons.
 */

#include "QedWrapperCommons.H"

//#define PXRMP_CORE_ONLY allows importing only the 'core functions' of the
//QED library.
#define PXRMP_CORE_ONLY
#include <chi_functions.hpp>

namespace QedUtils{
    /**
    * Function to calculate the 'chi' parameter for photons.
    * Suitable for GPU kernels.
    * @param[in] px,py,pz components of momentum (SI units)
    * @param[in] ex,ey,ez components of electric field (SI units)
    * @param[in] bx,by,bz components of magnetic field (SI units)
    * @return chi parameter
    */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real chi_photon(
        amrex::Real px, amrex::Real py, amrex::Real pz,
        amrex::Real ex, amrex::Real ey, amrex::Real ez,
        amrex::Real bx, amrex::Real by, amrex::Real bz)
    {
        //laser wavelength is unused if SI units are set
        const amrex::Real dummy_lambda = 1.0;
        return picsar::multi_physics::chi_photon(
            px, py, pz, ex, ey, ez, bx, by, bz, dummy_lambda);
    }

    /**
    * Function to calculate the 'chi' parameter for electrons or positrons.
    * Suitable for GPU kernels.
    * @param[in] px,py,pz components of momentum (SI units)
    * @param[in] ex,ey,ez components of electric field (SI units)
    * @param[in] bx,by,bz components of magnetic field (SI units)
    * @return chi parameter
    */
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real chi_lepton(
        amrex::Real px, amrex::Real py, amrex::Real pz,
        amrex::Real ex, amrex::Real ey, amrex::Real ez,
        amrex::Real bx, amrex::Real by, amrex::Real bz)
    {
        //laser wavelength is unused if SI units are set
        const amrex::Real dummy_lambda = 1.0;
        return picsar::multi_physics::chi_lepton(
            px, py, pz, ex, ey, ez, bx, by, bz, dummy_lambda);
    }
    //_________
};

#endif //WARPX_amrex_qed_chi_functions_h_
