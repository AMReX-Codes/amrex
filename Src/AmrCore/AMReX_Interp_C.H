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
face_linear_face_interp_x (int ci, int cj, int ck, int n, Array4<T> const& fine,
                           Array4<T const> const& crse, Array4<int const> const& mask,
                           IntVect const& ratio) noexcept
{
    // Fill only overlapping faces.
    int fi = ci*ratio[0];
    int fj = 0;
    int fk = 0;

#if AMREX_SPACEDIM >= 2
    fj = cj*ratio[1];
#if AMREX_SPACEDIM == 3
    fk = ck*ratio[2];
#endif
#endif

    // Check mask to ensure we don't overwrite valid fine data.
    if (!mask(ci, cj, ck, n)) {
        return;
    }

    fine(fi,fj  ,fk  ,n) = crse(ci,cj,ck,n);
#if AMREX_SPACEDIM >= 2
    fine(fi,fj+1,fk  ,n) = crse(ci,cj,ck,n);
#if AMREX_SPACEDIM == 3
    fine(fi,fj  ,fk+1,n) = crse(ci,cj,ck,n);
    fine(fi,fj+1,fk+1,n) = crse(ci,cj,ck,n);
#endif
#endif

}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_face_interp_y (int ci, int cj, int ck, int n, Array4<T> const& fine,
                           Array4<T const> const& crse, Array4<int const> const& mask,
                           IntVect const& ratio) noexcept
{
    // Fill only overlapping faces.
    int fi = ci*ratio[0];
    int fj = cj*ratio[1];
    int fk = 0;

#if AMREX_SPACEDIM == 3
    fk = ck*ratio[2];
#endif

    // Check mask to ensure we don't overwrite valid fine data.
    if (!mask(ci, cj, ck, n)) {
        return;
    }

    fine(fi  ,fj,fk  ,n) = crse(ci,cj,ck,n);
    fine(fi+1,fj,fk  ,n) = crse(ci,cj,ck,n);
#if AMREX_SPACEDIM == 3
    fine(fi  ,fj,fk+1,n) = crse(ci,cj,ck,n);
    fine(fi+1,fj,fk+1,n) = crse(ci,cj,ck,n);
#endif

}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_face_interp_z (int ci, int cj, int ck, int n, Array4<T> const& fine,
                           Array4<T const> const& crse, Array4<int const> const& mask,
                           IntVect const& ratio) noexcept
{
    // Fill only overlapping faces.
    int fi = ci*ratio[0];
    int fj = cj*ratio[1];
    int fk = ck*ratio[2];

    // Check mask to ensure we don't overwrite valid fine data.
    if (!mask(ci, cj, ck, n)) {
        return;
    }

    fine(fi  ,fj  ,fk,n) = crse(ci,cj,ck,n);
    fine(fi+1,fj  ,fk,n) = crse(ci,cj,ck,n);
    fine(fi  ,fj+1,fk,n) = crse(ci,cj,ck,n);
    fine(fi+1,fj+1,fk,n) = crse(ci,cj,ck,n);
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
        Real const w = static_cast<Real>(i-ci*ratio[0]) * (Real(1.)/ratio[0]);
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
        Real const w = static_cast<Real>(j-cj*ratio[1]) * (Real(1.)/ratio[1]);
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
        Real const w = static_cast<Real>(k-ck*ratio[2]) * (Real(1.)/ratio[2]);
        int k1 = ck*ratio[2];
        int k2 = (ck+1)*ratio[2];
        fine(i,j,k,n) = (Real(1.)-w) * fine(i,j,k1,n) + w * fine(i,j,k2,n);
    }
}

}
#endif
