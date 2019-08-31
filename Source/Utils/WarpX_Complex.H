#ifndef WARPX_COMPLEX_H_
#define WARPX_COMPLEX_H_

#include <AMReX_REAL.H>

// Define complex type on GPU/CPU
#ifdef AMREX_USE_GPU

#include <thrust/complex.h>
#include <cufft.h>
using Complex = thrust::complex<amrex::Real>;
static_assert( sizeof(Complex) == sizeof(cuDoubleComplex),
    "The complex types in WarpX and cuFFT do not match.");

#else

#include <complex>
#include <fftw3.h>
using Complex = std::complex<amrex::Real>;
static_assert( sizeof(Complex) == sizeof(fftw_complex),
    "The complex types in WarpX and FFTW do not match.");

#endif
static_assert(sizeof(Complex) == sizeof(amrex::Real[2]),
              "Unexpected complex type.");

#endif //WARPX_COMPLEX_H_
