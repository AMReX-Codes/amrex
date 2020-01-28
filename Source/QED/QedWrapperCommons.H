/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_amrex_qed_wrapper_commons_h_
#define WARPX_amrex_qed_wrapper_commons_h_

/**
 * This header contains some common #define directives and a
 * 'dummy' class used by the QED library wrappers and related
 * components.
 */

#include <AMReX_AmrCore.H>
#include <AMReX_GpuQualifiers.H>

/**
 * PICSAR uses PXRMP_GPU to decorate methods which should be
 * compiled for GPU. The user has to set it to the right value
 * (AMREX_GPU_DEVICE in this case).
 * PXRMP_WITH_SI_UNITS sets the library to use International
 * System units.
 */
#define PXRMP_GPU AMREX_GPU_HOST_DEVICE
#define PXRMP_WITH_SI_UNITS
//_________________________

/**
 * A namespace called 'QedUtils' is used to encapsulate
 * free functions (defined elsewhere) and an
 * empty datastructure (DummyStruct), which is re-used by several
 * components.
 */
namespace QedUtils{
    struct DummyStruct{};
};
//_________________________


#endif //WARPX_amrex_qed_wrapper_commons_h_
