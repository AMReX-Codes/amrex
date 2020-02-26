/* Copyright 2019 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_COMPLEXFORFFT_H_
#define WARPX_COMPLEXFORFFT_H_

#include "Utils/WarpX_Complex.H"


// Check the complex type on GPU/CPU
#ifdef AMREX_USE_GPU

#include <cufft.h>
static_assert( sizeof(Complex) == sizeof(cuDoubleComplex),
    "The complex types in WarpX and cuFFT do not match.");

#else
#   ifdef WARPX_USE_PSATD
#       include <fftw3.h>
static_assert( sizeof(Complex) == sizeof(fftw_complex),
    "The complex types in WarpX and FFTW do not match.");
#   endif // WARPX_USE_PSATD

#endif

#endif //WARPX_COMPLEXFORFFT_H_
