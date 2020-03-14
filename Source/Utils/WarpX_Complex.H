/* Copyright 2019-2020 Andrew Myers, David Grote, Maxence Thevenet
 * Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_COMPLEX_H_
#define WARPX_COMPLEX_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuComplex.H>

#include <complex>

// Define complex type on GPU/CPU
#ifdef AMREX_USE_GPU

using Complex = amrex::GpuComplex<amrex::Real>;

#    ifdef WARPX_USE_PSATD
#    include <cufft.h>
#        ifdef AMREX_USE_FLOAT
static_assert( sizeof(Complex) == sizeof(cuComplex),
    "The complex types in WarpX and cuFFT do not match.");
#        else
static_assert( sizeof(Complex) == sizeof(cuDoubleComplex),
    "The complex types in WarpX and cuFFT do not match.");
#        endif
#    endif // WARPX_USE_PSATD

#else

using Complex = std::complex<amrex::Real>;

#    ifdef WARPX_USE_PSATD
#    include <fftw3.h>
#        ifdef AMREX_USE_FLOAT
static_assert( sizeof(Complex) == sizeof(fftwf_complex),
    "The complex types in WarpX and FFTW do not match.");
#        else
static_assert( sizeof(Complex) == sizeof(fftw_complex),
    "The complex types in WarpX and FFTW do not match.");
#        endif
#    endif // WARPX_USE_PSATD

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
        return amrex::exp(val);
#else
        return std::exp(val);
#endif
    }

    // sqrt function
    template<typename T>
    AMREX_GPU_HOST_DEVICE T sqrt (const T& val){
#ifdef AMREX_USE_GPU
        return amrex::sqrt(val);
#else
        return std::sqrt(val);
#endif
    }

    // power function
    template<typename T1, typename T2>
    AMREX_GPU_HOST_DEVICE T1 pow (const T1& val, const T2& power){
#ifdef AMREX_USE_GPU
        return amrex::pow(val, power);
#else
        return std::pow(val, power);
#endif
    }
}

#endif //WARPX_COMPLEX_H_
