#ifndef AMREX_INTERP_C_H_
#define AMREX_INTERP_C_H_
#include <AMReX_Config.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_Interp_1D_C.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_Interp_2D_C.H>
#else
#include <AMReX_Interp_3D_C.H>
#endif


namespace amrex {

//
// Fill fine values with piecewise-constant interpolation of coarse data.
// Operate only on faces that overlap--ie, only fill the fine faces that make up each
// coarse face, leave the in-between faces alone.
//
template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_face_interp_x (int fi, int fj, int fk, int n, Array4<T> const& fine,
                           Array4<T const> const& crse, Array4<int const> const& mask,
                           IntVect const& ratio) noexcept
{
    int ci = amrex::coarsen(fi, ratio[0]);
    if (ci*ratio[0] == fi) {
#if (AMREX_SPACEDIM >= 2)
        int cj = amrex::coarsen(fj, ratio[1]);
#else
        int cj = 0;
#endif
#if (AMREX_SPACEDIM == 3)
        int ck = amrex::coarsen(fk, ratio[2]);
#else
        int ck = 0;
#endif

        // Check solve mask to ensure we don't overwrite valid fine data.
        if (!mask || mask(ci, cj, ck, n)) {
            fine(fi, fj, fk, n) = crse(ci, cj, ck, n);
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_face_interp_y (int fi, int fj, int fk, int n, Array4<T> const& fine,
                           Array4<T const> const& crse, Array4<int const> const& mask,
                           IntVect const& ratio) noexcept
{
    int cj = amrex::coarsen(fj, ratio[1]);
    if (cj*ratio[1] == fj) {
        int ci = amrex::coarsen(fi, ratio[0]);
#if (AMREX_SPACEDIM == 3)
        int ck = amrex::coarsen(fk, ratio[2]);
#else
        int ck = 0;
#endif

        // Check solve mask to ensure we don't overwrite valid fine data.
        if (!mask || mask(ci, cj, ck, n)) {
            fine(fi, fj, fk, n) = crse(ci, cj, ck, n);
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_face_interp_z (int fi, int fj, int fk, int n, Array4<T> const& fine,
                           Array4<T const> const& crse, Array4<int const> const& mask,
                           IntVect const& ratio) noexcept
{
    int ck = amrex::coarsen(fk, ratio[2]);
    if (ck*ratio[2] == fk) {
        int ci = amrex::coarsen(fi, ratio[0]);
        int cj = amrex::coarsen(fj, ratio[1]);

        // Check solve mask to ensure we don't overwrite valid fine data.
        if (!mask || mask(ci, cj, ck, n)) {
            fine(fi, fj, fk, n) = crse(ci, cj, ck, n);
        }
    }
}

//
// Do linear in dir, pc transverse to dir, leave alone the fine values
// lining up with coarse edges--assume these have been set to hold the
// values you want to interpolate to the rest.
//
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_x (int i, int j, int k, int n, amrex::Array4<amrex::Real> const& fine,
                      IntVect const& ratio) noexcept
{
    const int ci = amrex::coarsen(i,ratio[0]);

    if (i-ci*ratio[0] != 0)
    {
        Real const w = static_cast<Real>(i-ci*ratio[0]) * (Real(1.)/Real(ratio[0]));
        int i1 = ci*ratio[0];
        int i2 = (ci+1)*ratio[0];
        fine(i,j,k,n) = (Real(1.)-w) * fine(i1,j,k,n) + w * fine(i2,j,k,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_y (int i, int j, int k, int n, amrex::Array4<amrex::Real> const& fine,
                      IntVect const& ratio) noexcept
{
    const int cj = amrex::coarsen(j,ratio[1]);

    if (j-cj*ratio[1] != 0)
    {
        Real const w = static_cast<Real>(j-cj*ratio[1]) * (Real(1.)/Real(ratio[1]));
        int j1 = cj*ratio[1];
        int j2 = (cj+1)*ratio[1];
        fine(i,j,k,n) = (Real(1.)-w) * fine(i,j1,k,n) + w * fine(i,j2,k,n);
    }

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_z (int i, int j, int k, int n, amrex::Array4<amrex::Real> const& fine,
                      IntVect const& ratio) noexcept
{
    const int ck = amrex::coarsen(k,ratio[2]);

    if (k-ck*ratio[2] != 0)
    {
        Real const w = static_cast<Real>(k-ck*ratio[2]) * (Real(1.)/Real(ratio[2]));
        int k1 = ck*ratio[2];
        int k2 = (ck+1)*ratio[2];
        fine(i,j,k,n) = (Real(1.)-w) * fine(i,j,k1,n) + w * fine(i,j,k2,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cell_quartic_interp_x (int i, int j, int k, int n, Array4<Real> const& fine,
                            Array4<Real const> const& crse) noexcept
{
    constexpr Array1D<Real,-2,2> c = {Real(0.01708984), Real(-0.12304688),
                                      Real(0.92285156), Real(0.20507812),
                                      Real(-0.02197266)};
    int ii = amrex::coarsen(i,2);
    int s = 2*(i-ii*2) - 1;  // if i == ii*2, s = -1; if i == ii*2+1, s = 1;
    fine(i,j,k,n) = c(-2*s)*crse(ii-2,j,k,n)
        +           c(  -s)*crse(ii-1,j,k,n)
        +           c(   0)*crse(ii  ,j,k,n)
        +           c(   s)*crse(ii+1,j,k,n)
        +           c( 2*s)*crse(ii+2,j,k,n);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cell_quartic_interp_y (int i, int j, int k, int n, Array4<Real> const& fine,
                            Array4<Real const> const& crse) noexcept
{
    constexpr Array1D<Real,-2,2> c = {Real(0.01708984), Real(-0.12304688),
                                      Real(0.92285156), Real(0.20507812),
                                      Real(-0.02197266)};
    int jj = amrex::coarsen(j,2);
    int s = 2*(j-jj*2) - 1;  // if j == jj*2, s = -1; if j == jj*2+1, s = 1;
    fine(i,j,k,n) = c(-2*s)*crse(i,jj-2,k,n)
        +           c(  -s)*crse(i,jj-1,k,n)
        +           c(   0)*crse(i,jj  ,k,n)
        +           c(   s)*crse(i,jj+1,k,n)
        +           c( 2*s)*crse(i,jj+2,k,n);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cell_quartic_interp_z (int i, int j, int k, int n, Array4<Real> const& fine,
                            Array4<Real const> const& crse) noexcept
{
    constexpr Array1D<Real,-2,2> c = {Real(0.01708984), Real(-0.12304688),
                                      Real(0.92285156), Real(0.20507812),
                                      Real(-0.02197266)};
    int kk = amrex::coarsen(k,2);
    int s = 2*(k-kk*2) - 1;  // if k == kk*2, s = -1; if k == kk*2+1, s = 1;
    fine(i,j,k,n) = c(-2*s)*crse(i,j,kk-2,n)
        +           c(  -s)*crse(i,j,kk-1,n)
        +           c(   0)*crse(i,j,kk  ,n)
        +           c(   s)*crse(i,j,kk+1,n)
        +           c( 2*s)*crse(i,j,kk+2,n);
}

}
#endif
