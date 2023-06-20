#ifndef AMREX_INTERP_BNDRYDATA_2D_K_H_
#define AMREX_INTERP_BNDRYDATA_2D_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_LOUtil_K.H>

namespace amrex {

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpbndrydata_o1 (int i, int j, int /*k*/, int n,
                         Array4<T> const& bdry, int nb,
                         Array4<T const> const& crse, int nc, Dim3 const& r) noexcept
{
    int ic = amrex::coarsen(i,r.x);
    int jc = amrex::coarsen(j,r.y);
    bdry(i,j,0,n+nb) = crse(ic,jc,0,n+nc);
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpbndrydata_x_o3 (int i, int j, int /*k*/, int n,
                           Array4<T> const& bdry, int nb,
                           Array4<T const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    int ic = amrex::coarsen(i,r.x);
    int jc = amrex::coarsen(j,r.y);

    T x[3], y[3], c[3];
    x[0] = T(0.0);
    y[0] = crse(ic,jc,0,n+nc);
    int NN = 1;

    if (mask(i,j-r.y,0) == not_covered) {
        x[NN] = T(-1.0);
        y[NN] = crse(ic,jc-1,0,n+nc);
        ++NN;
    } else if (mask.contains(i,(jc+2)*r.y,0) &&
               mask         (i,(jc+2)*r.y,0) == not_covered && crse.contains(ic,jc+2,0)) {
        x[NN] = T(2.0);
        y[NN] = crse(ic,jc+2,0,n+nc);
        ++NN;
    }

    if (mask(i,j+r.y,0) == not_covered) {
        x[NN] = T(1.0);
        y[NN] = crse(ic,jc+1,0,n+nc);
        ++NN;
    } else if (mask.contains(i,jc*r.y-r.y-1,0) &&
               mask         (i,jc*r.y-r.y-1,0) == not_covered && crse.contains(ic,jc-2,0)) {
        x[NN] = T(-2.0);
        y[NN] = crse(ic,jc-2,0,n+nc);
        ++NN;
    }

    if ( (mask(i,j-r.y,0) != not_covered) && (mask(i,j+r.y,0) != not_covered) ) {
        NN = 1;
    }

    T xInt = -T(0.5) + (j-jc*r.y+T(0.5))/r.y;
    poly_interp_coeff(xInt, x, NN, c);
    T b = T(0.0);
    for (int m = 0; m < NN; ++m) {
        b += c[m]*y[m];
    }
    bdry(i,j,0,n+nb) = b;
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpbndrydata_y_o3 (int i, int j, int /*k*/, int n,
                           Array4<T> const& bdry, int nb,
                           Array4<T const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    int ic = amrex::coarsen(i,r.x);
    int jc = amrex::coarsen(j,r.y);

    T x[3], y[3], c[3];
    x[0] = T(0.0);
    y[0] = crse(ic,jc,0,n+nc);
    int NN = 1;

    if (mask(i-r.x,j,0) == not_covered) {
        x[NN] = T(-1.0);
        y[NN] = crse(ic-1,jc,0,n+nc);
        ++NN;
    } else if (mask.contains((ic+2)*r.x,j,0) &&
               mask         ((ic+2)*r.x,j,0) == not_covered && crse.contains(ic+2,jc,0)) {
        x[NN] = T(2.0);
        y[NN] = crse(ic+2,jc,0,n+nc);
        ++NN;
    }

    if (mask(i+r.x,j,0) == not_covered) {
        x[NN] = T(1.0);
        y[NN] = crse(ic+1,jc,0,n+nc);
        ++NN;
    } else if (mask.contains(ic*r.x-r.x-1,j,0) &&
               mask         (ic*r.x-r.x-1,j,0) == not_covered && crse.contains(ic-2,jc,0)) {
        x[NN] = T(-2.0);
        y[NN] = crse(ic-2,jc,0,n+nc);
        ++NN;
    }

    if ( (mask(i-r.x,j,0) != not_covered) && (mask(i+r.x,j,0) != not_covered) ) {
        NN = 1;
    }

    T xInt = -T(0.5) + (i-ic*r.x+T(0.5))/r.x;
    poly_interp_coeff(xInt, x, NN, c);
    T b = T(0.0);
    for (int m = 0; m < NN; ++m) {
        b += c[m]*y[m];
    }
    bdry(i,j,0,n+nb) = b;
}

}
#endif
