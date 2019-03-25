#ifndef AMREX_LO_UTIL_K_H_
#define AMREX_LO_UTIL_K_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace amrex {

//  polyInterpCoeff:
//
//  This routine returns the Lagrange interpolating coefficients for a
//  polynomial through N points, evaluated at xInt (see Numerical Recipes,
//  v2, p102, e.g.):
//
//          (x-x2)(x-x3)...(x-xN)              (x-x1)(x-x2)...(x-x(N-1))
//  P(x) = ----------------------- y1  + ... + ------------------------  yN
//         (x1-x2)(x1-x3)...(x1-xN)            (x1-x2)(x1-x3)...(x1-xN)
//
//  P(xInt) = sum_(i=1)^(N) y[i]*c[i]

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void poly_interp_coeff (Real xInt, Real const* AMREX_RESTRICT x, int N, Real* AMREX_RESTRICT c) noexcept
{
    for (int j = 0; j < N; ++j) {
        Real num = 1.0, den = 1.0;
        for (int i = 0; i < N; ++i) {
            if (i != j) {
                num *= xInt-x[i];
                den *= x[j]-x[i];
            }
        }
        c[j] = num / den;
    }
}

}

#endif
