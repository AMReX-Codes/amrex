#ifndef WARPX_COMPLEX_H_
#define WARPX_COMPLEX_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>

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

// wrapper around math functions, to run on CPU or accelerator.
namespace MathFunc
{
    // exp function
    template<typename T>
    AMREX_GPU_HOST_DEVICE T exp (const T& val){
#ifdef AMREX_USE_GPU
        return thrust::exp(val);
#else
        return std::exp(val);
#endif
    }
  
    // sqrt function
    template<typename T>
    AMREX_GPU_HOST_DEVICE T sqrt (const T& val){
#ifdef AMREX_USE_GPU
        return thrust::sqrt(val);
#else
        return std::sqrt(val);
#endif
    }

    // power function
    template<typename T1, typename T2>
    AMREX_GPU_HOST_DEVICE T1 pow (const T1& val, const T2& power){
#ifdef AMREX_USE_GPU
        return thrust::pow(val, power);
#else
        return std::pow(val, power);
#endif
    }
}

#endif //WARPX_COMPLEX_H_
