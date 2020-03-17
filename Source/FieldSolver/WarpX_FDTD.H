/* Copyright 2019 Axel Huebl, David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_FDTD_H_
#define WARPX_FDTD_H_

#include <AMReX_FArrayBox.H>

static void warpx_calculate_ckc_coefficients(amrex::Real dtsdx, amrex::Real dtsdy, amrex::Real dtsdz,
                                             amrex::Real &betaxy, amrex::Real &betaxz, amrex::Real &betayx,
                                             amrex::Real &betayz, amrex::Real &betazx, amrex::Real &betazy,
                                             amrex::Real &gammax, amrex::Real &gammay, amrex::Real &gammaz,
                                             amrex::Real &alphax, amrex::Real &alphay, amrex::Real &alphaz)
{
    using namespace amrex;

    // Cole-Karkkainen-Cowan push
    // computes coefficients according to Cowan - PRST-AB 16, 041303 (2013)
#if defined WARPX_DIM_3D
    const Real delta = std::max( { dtsdx,dtsdy,dtsdz } );
    const Real rx = (dtsdx/delta)*(dtsdx/delta);
    const Real ry = (dtsdy/delta)*(dtsdy/delta);
    const Real rz = (dtsdz/delta)*(dtsdz/delta);
    const Real beta = 0.125*(1. - rx*ry*rz/(ry*rz + rz*rx + rx*ry));
    betaxy = ry*beta;
    betaxz = rz*beta;
    betayx = rx*beta;
    betayz = rz*beta;
    betazx = rx*beta;
    betazy = ry*beta;
    gammax = ry*rz*(0.0625 - 0.125*ry*rz/(ry*rz + rz*rx + rx*ry));
    gammay = rx*rz*(0.0625 - 0.125*rx*rz/(ry*rz + rz*rx + rx*ry));
    gammaz = rx*ry*(0.0625 - 0.125*rx*ry/(ry*rz + rz*rx + rx*ry));
    alphax = 1. - 2.*betaxy - 2.*betaxz - 4.*gammax;
    alphay = 1. - 2.*betayx - 2.*betayz - 4.*gammay;
    alphaz = 1. - 2.*betazx - 2.*betazy - 4.*gammaz;

    betaxy *= dtsdx;
    betaxz *= dtsdx;
    betayx *= dtsdy;
    betayz *= dtsdy;
    betazx *= dtsdz;
    betazy *= dtsdz;
    alphax *= dtsdx;
    alphay *= dtsdy;
    alphaz *= dtsdz;
    gammax *= dtsdx;
    gammay *= dtsdy;
    gammaz *= dtsdz;
#elif defined WARPX_DIM_XZ
    const Real delta = std::max(dtsdx,dtsdz);
    const Real rx = (dtsdx/delta)*(dtsdx/delta);
    const Real rz = (dtsdz/delta)*(dtsdz/delta);
    betaxz = 0.125*rz;
    betazx = 0.125*rx;
    alphax = 1. - 2.*betaxz;
    alphaz = 1. - 2.*betazx;

    betaxz *= dtsdx;
    betazx *= dtsdz;
    alphax *= dtsdx;
    alphaz *= dtsdz;
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_computedivb(int i, int j, int k, int dcomp,
                       amrex::Array4<amrex::Real> const& divB,
                       amrex::Array4<amrex::Real const> const& Bx,
                       amrex::Array4<amrex::Real const> const& By,
                       amrex::Array4<amrex::Real const> const& Bz,
                       amrex::Real dxinv,
                       amrex::Real dyinv,
                       amrex::Real dzinv
#ifdef WARPX_DIM_RZ
                       , amrex::Real const rmin
#endif
                       )
{
#if defined WARPX_DIM_3D
    divB(i,j,k,dcomp) = (Bx(i+1,j  ,k  ) - Bx(i,j,k))*dxinv
               +        (By(i  ,j+1,k  ) - By(i,j,k))*dyinv
               +        (Bz(i  ,j  ,k+1) - Bz(i,j,k))*dzinv;
#elif defined WARPX_DIM_XZ
    divB(i,j,0,dcomp) = (Bx(i+1,j  ,0) - Bx(i,j,0))*dxinv
               +        (Bz(i  ,j+1,0) - Bz(i,j,0))*dzinv;
#elif defined WARPX_DIM_RZ
    const amrex::Real ru = 1. + 0.5/(rmin*dxinv + i + 0.5);
    const amrex::Real rd = 1. - 0.5/(rmin*dxinv + i + 0.5);
    divB(i,j,0,dcomp) = (ru*Bx(i+1,j,0) - rd*Bx(i,j,0))*dxinv
                       + (Bz(i,j+1,0) - Bz(i,j,0))*dzinv;
#endif
}

#endif
