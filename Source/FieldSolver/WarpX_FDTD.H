#ifndef WARPX_FDTD_H_
#define WARPX_FDTD_H_

#include <AMReX_FArrayBox.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_computedivb(int i, int j, int k, int dcomp, Array4<Real> const& divB,
                       Array4<Real const> const& Bx, Array4<Real const> const& By, Array4<Real const> const& Bz,
                       Real dxinv, Real dyinv, Real dzinv
#ifdef WARPX_RZ
                       ,const Real rmin
#endif
                       )
{
#if (AMREX_SPACEDIM == 3)
    divB(i,j,k,dcomp) = (Bx(i+1,j  ,k  ) - Bx(i,j,k))*dxinv
               +        (By(i  ,j+1,k  ) - By(i,j,k))*dyinv
               +        (Bz(i  ,j  ,k+1) - Bz(i,j,k))*dzinv;
#else
#ifndef WARPX_RZ
    divB(i,j,0,dcomp) = (Bx(i+1,j  ,0) - Bx(i,j,0))*dxinv
               +        (Bz(i  ,j+1,0) - Bz(i,j,0))*dzinv;
#else
    Real ru = 1. + 0.5/(rmin*dxinv + i + 0.5);
    Real rd = 1. - 0.5/(rmin*dxinv + i + 0.5);
    divB(i,j,0,dcomp) = (ru*Bx(i+1,j,0) - rd*Bx(i,j,0))*dxinv
                       + (Bz(i,j+1,0) - Bz(i,j,0))*dzinv;
#endif
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_computedive(int i, int j, int k, int dcomp, Array4<Real> const& divE,
                       Array4<Real const> const& Ex, Array4<Real const> const& Ey, Array4<Real const> const& Ez,
                       Real dxinv, Real dyinv, Real dzinv
#ifdef WARPX_RZ
                       ,const Real rmin
#endif
                       )
{
#if (AMREX_SPACEDIM == 3)
    divE(i,j,k,dcomp) = (Ex(i,j,k) - Ex(i-1,j,k))*dxinv
               +        (Ey(i,j,k) - Ey(i,j-1,k))*dyinv
               +        (Ez(i,j,k) - Ez(i,j,k-1))*dzinv;
#else
#ifndef WARPX_RZ
    divE(i,j,0,dcomp) = (Ex(i,j,0) - Ex(i-1,j,0))*dxinv
                 +      (Ez(i,j,0) - Ez(i,j-1,0))*dzinv;
#else
    if (i == 0 && rmin == 0.) {
        // the bulk equation diverges on axis
        // (due to the 1/r terms). The following expressions regularize
        // these divergences.
        divE(i,j,0,dcomp) = 4.*Ex(i,j,0)*dxinv
                            + (Ez(i,j,0) - Ez(i,j-1,0))*dzinv;
    } else {
        Real ru = 1. + 0.5/(rmin*dxinv + i);
        Real rd = 1. - 0.5/(rmin*dxinv + i);
        divE(i,j,0,dcomp) = (ru*Ex(i,j,0) - rd*Ex(i-1,j,0))*dxinv
                           + (Ez(i,j,0) - Ez(i,j-1,0))*dzinv;
    }
#endif
#endif
}

#endif
