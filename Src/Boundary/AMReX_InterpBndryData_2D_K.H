#ifndef AMREX_INTERP_BNDRYDATA_2D_K_H_
#define AMREX_INTERP_BNDRYDATA_2D_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_LOUtil_K.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_x_o1 (int islo, int ic, int jc, int /*kc*/, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r) noexcept
{
    Real c = crse(ic,jc,0,n+nc);
    const int i = ic*r.x + islo*(r.x-1);
    for (int joff = 0; joff < r.y; ++joff) {
        int  j    =        jc*r.y  + joff;
        bdry(i,j,0,n+nb) = c;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_y_o1 (int islo, int ic, int jc, int /*kc*/, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r) noexcept
{
    Real c = crse(ic,jc,0,n+nc);
    const int j = jc*r.y + islo*(r.y-1);
    for (int ioff = 0; ioff < r.x; ++ioff) {
        int  i    =        ic*r.x  + ioff;
        bdry(i,j,0,n+nb) = c;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_x_o3 (int islo, int ic, int jc, int /*kc*/, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    const int i = ic*r.x + islo*(r.x-1);
    int j = jc*r.y;

    Real x[3], y[3], c[3];
    x[0] = 0.0;
    y[0] = crse(ic,jc,0,n+nc);
    int NN = 1;

    if (mask(i,j-1,0) == not_covered) {
        x[NN] = -1.0;
        y[NN] = crse(ic,jc-1,0,n+nc);
        ++NN;
    } else if (mask(i,j+2*r.y,0) == not_covered) {
        x[NN] = 2.0;
        y[NN] = crse(ic,jc+2,0,n+nc);
        ++NN;
    }

    if (mask(i,j+r.y,0) == not_covered) {
        x[NN] = 1.0;
        y[NN] = crse(ic,jc+1,0,n+nc);
        ++NN;
    } else if (mask(i,j-r.y-1,0) == not_covered) {
        x[NN] = -2.0;
        y[NN] = crse(ic,jc-2,0,n+nc);
        ++NN;
    }

    if ( (mask(i,j-1,0) != not_covered) and (mask(i,j+r.y,0) != not_covered) ) {
        NN = 1;
    }

    for (int off = 0; off < r.y; ++off) {
        Real xInt = ((off-r.y/2) + 0.5) / r.y;
        poly_interp_coeff(xInt, x, NN, c);
        Real b = 0.0;
        for (int m = 0; m < NN; ++m) {
            b += c[m]*y[m];
        }
        bdry(i,j+off,0,n+nb) = b;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_y_o3 (int islo, int ic, int jc, int /*kc*/, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    int i = ic*r.x;
    const int j = jc*r.y + islo*(r.y-1);

    Real x[3], y[3], c[3];
    x[0] = 0.0;
    y[0] = crse(ic,jc,0,n+nc);
    int NN = 1;

    if (mask(i-1,j,0) == not_covered) {
        x[NN] = -1.0;
        y[NN] = crse(ic-1,jc,0,n+nc);
        ++NN;
    } else if (mask(i+2*r.x,j,0) == not_covered) {
        x[NN] = 2.0;
        y[NN] = crse(ic+2,jc,0,n+nc);
        ++NN;
    }

    if (mask(i+r.x,j,0) == not_covered) {
        x[NN] = 1.0;
        y[NN] = crse(ic+1,jc,0,n+nc);
        ++NN;
    } else if (mask(i-r.x-1,j,0) == not_covered) {
        x[NN] = -2.0;
        y[NN] = crse(ic-2,jc,0,n+nc);
        ++NN;
    }

    if ( (mask(i-1,j,0) != not_covered) and (mask(i+r.x,j,0) != not_covered) ) {
        NN = 1;
    }

    for (int off = 0; off < r.x; ++off) {
        Real xInt = ((off-r.x/2) + 0.5) / r.x;
        poly_interp_coeff(xInt, x, NN, c);
        Real b = 0.0;
        for (int m = 0; m < NN; ++m) {
            b += c[m]*y[m];
        }
        bdry(i+off,j,0,n+nb) = b;
    }
}

}
#endif
