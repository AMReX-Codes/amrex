#ifndef WARPX_COMPLEXFORFFT_H_
#define WARPX_COMPLEXFORFFT_H_

#include <WarpX_Complex.H>

// Check the complex type on GPU/CPU
#ifdef AMREX_USE_GPU

#include <cufft.h>
static_assert( sizeof(Complex) == sizeof(cuDoubleComplex),
    "The complex types in WarpX and cuFFT do not match.");

#else

#include <fftw3.h>
static_assert( sizeof(Complex) == sizeof(fftw_complex),
    "The complex types in WarpX and FFTW do not match.");

#endif

#endif //WARPX_COMPLEXFORFFT_H_
