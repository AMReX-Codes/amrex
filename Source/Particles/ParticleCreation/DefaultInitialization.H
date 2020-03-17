/* Copyright 2019-2020 Andrew Myers, Axel Huebl,
 * Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef DEFAULTINITIALIZATION_H_
#define DEFAULTINITIALIZATION_H_

#include "AMReX_REAL.H"
#include "AMReX_GpuContainers.H"

#include <map>
#include <string>
#include <cmath>

/**
 * \brief This set of initialization policies describes what happens
 * when we need to create a new particle due to an elementary process.
 * For example, when an ionization event creates an electron, these
 * policies control the initial values of the electron's components.
 * These can always be over-written later.
 *
 * The specific meanings are as follows:
 *     Zero         - set the component to zero
 *     One          - set the component to one
 *     RandomExp    - a special flag for the optical depth component used by
 *                    certain QED processes, which gets a random initial value
 *                    extracted from an exponential distribution
 *
 */
enum struct InitializationPolicy {Zero=0, One, RandomExp};

/**
 * \brief This map sets the initialization policy for each particle component
 * used in WarpX.
 */
static std::map<std::string, InitializationPolicy> initialization_policies = {
    {"w",     InitializationPolicy::Zero },
    {"ux",    InitializationPolicy::Zero },
    {"uy",    InitializationPolicy::Zero },
    {"uz",    InitializationPolicy::Zero },
    {"Ex",    InitializationPolicy::Zero },
    {"Ey",    InitializationPolicy::Zero },
    {"Ez",    InitializationPolicy::Zero },
    {"Bx",    InitializationPolicy::Zero },
    {"By",    InitializationPolicy::Zero },
    {"Bz",    InitializationPolicy::Zero },
#ifdef WARPX_DIM_RZ
    {"theta", InitializationPolicy::Zero},
#endif

#ifdef WARPX_QED
    {"optical_depth_BW",   InitializationPolicy::RandomExp},
    {"optical_depth_QSR",   InitializationPolicy::RandomExp}
#endif

};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::ParticleReal initializeRealValue (const InitializationPolicy policy) noexcept
{
    switch (policy) {
        case InitializationPolicy::Zero : return 0.0;
        case InitializationPolicy::One  : return 1.0;
        case InitializationPolicy::RandomExp : {
            return -log(amrex::Random());
        }
        default : {
            amrex::Abort("Initialization Policy not recognized");
            return 1.0;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int initializeIntValue (const InitializationPolicy policy) noexcept
{
    switch (policy) {
        case InitializationPolicy::Zero : return 0;
        case InitializationPolicy::One  : return 1;
        default : {
            amrex::Abort("Initialization Policy not recognized");
            return 1;
        }
    }
}

#endif
