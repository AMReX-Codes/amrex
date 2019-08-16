#ifndef WARPX_COMPLEX_H_
#define WARPX_COMPLEX_H_

#include <AMReX_REAL.H>

// Define complex type on GPU/CPU
#ifdef AMREX_USE_GPU

#include <thrust/complex.h>
using Complex = thrust::complex<amrex::Real>;

#else

#include <complex>
using Complex = std::complex<amrex::Real>;

#endif

static_assert(sizeof(Complex) == sizeof(amrex::Real[2]),
              "Unexpected complex type.");

#endif //WARPX_COMPLEX_H_
