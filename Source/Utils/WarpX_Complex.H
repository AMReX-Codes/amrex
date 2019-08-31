#ifndef WARPX_COMPLEX_H_
#define WARPX_COMPLEX_H_

#include <AMReX_REAL.H>

// Define complex type on GPU/CPU
#ifdef AMREX_USE_GPU

#include <thrust/complex.h>

using Complex = thrust::complex<amrex::Real>;

#ifdef WARPX_USE_PSATD
#include <cufft.h>
static_assert( sizeof(Complex) == sizeof(cuDoubleComplex),
    "The complex types in WarpX and cuFFT do not match.");
#endif // WARPX_USE_PSATD

#else

#include <complex>
using Complex = std::complex<amrex::Real>;

#ifdef WARPX_USE_PSATD
#include <fftw3.h>
static_assert( sizeof(Complex) == sizeof(fftw_complex),
    "The complex types in WarpX and FFTW do not match.");
#endif // WARPX_USE_PSATD

#endif // AMREX_USE_GPU
static_assert(sizeof(Complex) == sizeof(amrex::Real[2]),
              "Unexpected complex type.");

#endif //WARPX_COMPLEX_H_
