/* Copyright 2019 Axel Huebl, David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_FDTD_H_
#define WARPX_FDTD_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ex_yee(int i, int j, int k,
                       amrex::Array4<amrex::Real> const& Ex,
                       amrex::Array4<amrex::Real const> const& By,
                       amrex::Array4<amrex::Real const> const& Bz,
                       amrex::Array4<amrex::Real const> const& Jx,
                       amrex::Real mu_c2_dt, amrex::Real dtsdx_c2,
                       amrex::Real dtsdy_c2, amrex::Real dtsdz_c2,
                       amrex::Real dxinv, amrex::Real rmin,
                       const long nmodes)
{
#if defined WARPX_DIM_3D
    Ex(i,j,k) += + dtsdy_c2 * (Bz(i,j,k) - Bz(i,j-1,k  ))
                 - dtsdz_c2 * (By(i,j,k) - By(i,j  ,k-1))
                 - mu_c2_dt  * Jx(i,j,k);
#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    // Note that the 2D Cartesian and RZ mode 0 are the same
    Ex(i,j,0,0) += - dtsdz_c2 * (By(i,j,0,0) - By(i,j-1,0,0))
                   - mu_c2_dt  * Jx(i,j,0,0);
#if (defined WARPX_DIM_RZ)
    const amrex::Real r = rmin*dxinv+ i + 0.5;
    for (int imode=1 ; imode < nmodes ; imode++) {
        // Er(i,j,m) = Er(i,j,m) - I*m*dt*Bz(i,j,m)/r - dtsdz*(Bt(i,j,m) - Bt(i,j-1,m)) - mudt*Jr(i,j,m)
        Ex(i,j,0,2*imode-1) += - dtsdz_c2*(By(i,j,0,2*imode-1) - By(i,j-1,0,2*imode-1)) + imode*dtsdx_c2*Bz(i,j,0,2*imode)/r
                             - mu_c2_dt*Jx(i,j,0,2*imode-1);
        Ex(i,j,0,2*imode) += - dtsdz_c2*(By(i,j,0,2*imode) - By(i,j-1,0,2*imode)) - imode*dtsdx_c2*Bz(i,j,0,2*imode-1)/r
                             - mu_c2_dt*Jx(i,j,0,2*imode);
    }
#endif
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ey_yee(int i, int j, int k,
                       amrex::Array4<amrex::Real> const& Ey,
                       amrex::Array4<amrex::Real const> const& Bx,
                       amrex::Array4<amrex::Real const> const& Bz,
                       amrex::Array4<amrex::Real const> const& Jy,
                       amrex::Array4<amrex::Real> const& Ex,
                       amrex::Real mu_c2_dt, amrex::Real dtsdx_c2,
                       amrex::Real dtsdz_c2, amrex::Real rmin,
                       const long nmodes)
{
#if defined WARPX_DIM_3D
    Ey(i,j,k) += - dtsdx_c2 * (Bz(i,j,k) - Bz(i-1,j,k))
                 + dtsdz_c2 * (Bx(i,j,k) - Bx(i,j,k-1))
                 - mu_c2_dt  * Jy(i,j,k);
#elif (defined WARPX_DIM_XZ)
    Ey(i,j,0) += - dtsdx_c2 * (Bz(i,j,0) - Bz(i-1,j,0))
                 + dtsdz_c2 * (Bx(i,j,0) - Bx(i,j-1,0))
                 - mu_c2_dt  * Jy(i,j,0);
#elif (defined WARPX_DIM_RZ)
    if (i != 0 || rmin != 0.) {
        Ey(i,j,0,0) += - dtsdx_c2 * (Bz(i,j,0,0) - Bz(i-1,j,0,0))
                       + dtsdz_c2 * (Bx(i,j,0,0) - Bx(i,j-1,0,0))
                       - mu_c2_dt  * Jy(i,j,0,0);
    } else {
        Ey(i,j,0,0) = 0.;
    }
    for (int imode=1 ; imode < nmodes ; imode++) {
        if (i == 0 && rmin == 0) {
            if (imode == 1) {
                // The bulk equation could in principle be used here since it does not diverge
                // on axis. However, it typically gives poor results e.g. for the propagation
                // of a laser pulse (the field is spuriously reduced on axis). For this reason
                // a modified on-axis condition is used here: we use the fact that
                // Etheta(r=0,m=1) should equal -iEr(r=0,m=1), for the fields Er and Et to be
                // independent of theta at r=0. Now with linear interpolation:
                // Er(r=0,m=1) = 0.5*[Er(r=dr/2,m=1) + Er(r=-dr/2,m=1)]
                // And using the rule applying for the guards cells
                // Er(r=-dr/2,m=1) = Er(r=dr/2,m=1). Thus: Et(i,j,m) = -i*Er(i,j,m)
                Ey(i,j,0,2*imode-1) = Ex(i,j,0,2*imode);
                Ey(i,j,0,2*imode) = -Ex(i,j,0,2*imode-1);
            } else {
                // Etheta should remain 0 on axis, for modes different than m=1
                Ey(i,j,0,2*imode-1) = 0.;
                Ey(i,j,0,2*imode) = 0.;
            }
        } else {
            // Et(i,j,m) = Et(i,j,m) - dtsdr*(Bz(i,j,m) - Bz(i-1,j,m)) + dtsdz*(Br(i,j,m) - Br(i,j-1,m)) - mudt*Jt(i,j,m)
            Ey(i,j,0,2*imode-1) += - dtsdx_c2 * (Bz(i,j,0,2*imode-1) - Bz(i-1,j,0,2*imode-1))
                           + dtsdz_c2 * (Bx(i,j,0,2*imode-1) - Bx(i,j-1,0,2*imode-1))
                           - mu_c2_dt  * Jy(i,j,0,2*imode-1);
            Ey(i,j,0,2*imode) += - dtsdx_c2 * (Bz(i,j,0,2*imode) - Bz(i-1,j,0,2*imode))
                           + dtsdz_c2 * (Bx(i,j,0,2*imode) - Bx(i,j-1,0,2*imode))
                           - mu_c2_dt  * Jy(i,j,0,2*imode);
        }
    }
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ez_yee(int i, int j, int k,
                       amrex::Array4<amrex::Real> const& Ez,
                       amrex::Array4<amrex::Real const> const& Bx,
                       amrex::Array4<amrex::Real const> const& By,
                       amrex::Array4<amrex::Real const> const& Jz,
                       amrex::Real mu_c2_dt,
                       amrex::Real dtsdx_c2, amrex::Real dtsdy_c2,
                       amrex::Real dxinv, amrex::Real rmin,
                       const long nmodes)
{
#if defined WARPX_DIM_3D
    Ez(i,j,k) += + dtsdx_c2 * (By(i,j,k) - By(i-1,j  ,k))
                 - dtsdy_c2 * (Bx(i,j,k) - Bx(i  ,j-1,k))
                 - mu_c2_dt  * Jz(i,j,k);
#elif defined WARPX_DIM_XZ
    Ez(i,j,0) += + dtsdx_c2 * (By(i,j,0) - By(i-1,j,0))
                 - mu_c2_dt  * Jz(i,j,0);
#elif defined WARPX_DIM_RZ
    if (i != 0 || rmin != 0.) {
        const amrex::Real ru = 1. + 0.5/(rmin*dxinv + i);
        const amrex::Real rd = 1. - 0.5/(rmin*dxinv + i);
        Ez(i,j,0,0) += + dtsdx_c2 * (ru*By(i,j,0,0) - rd*By(i-1,j,0,0))
                       - mu_c2_dt  * Jz(i,j,0,0);
    } else {
        Ez(i,j,0,0) += + 4.*dtsdx_c2 * By(i,j,0,0)
                       - mu_c2_dt  * Jz(i,j,0,0);
    }
    for (int imode=1 ; imode < nmodes ; imode++) {
        if (i == 0 && rmin == 0) {
            Ez(i,j,0,2*imode-1) = 0.;
            Ez(i,j,0,2*imode) = 0.;
        } else {
            const amrex::Real r = rmin*dxinv + i + 0.5;
            const amrex::Real ru = 1. + 0.5/(rmin*dxinv + i);
            const amrex::Real rd = 1. - 0.5/(rmin*dxinv + i);
            // Ez(i,j,m) = Ez(i,j,m) + dtsdr*(ru*Bt(i,j,m) - rd*Bt(i-1,j,m)) + I*m*dt*Br(i,j,m)/r - mudt*Jz(i,j,m)
            Ez(i,j,0,2*imode-1) += + dtsdx_c2 * (ru*By(i,j,0,2*imode-1) - rd*By(i-1,j,0,2*imode-1))
                                 - imode*dtsdx_c2*Bx(i,j,0,2*imode)/r
                                 - mu_c2_dt  * Jz(i,j,0,2*imode-1);
            Ez(i,j,0,2*imode) += + dtsdx_c2 * (ru*By(i,j,0,2*imode) - rd*By(i-1,j,0,2*imode))
                                 + imode*dtsdx_c2*Bx(i,j,0,2*imode-1)/r
                                 - mu_c2_dt  * Jz(i,j,0,2*imode);
        }
    }
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ex_f_yee(int j, int k, int l,
                         amrex::Array4<amrex::Real> const& Ex,
                         amrex::Array4<amrex::Real const> const& F,
                         amrex::Real dtsdx_c2)
{
#if defined WARPX_DIM_3D
    Ex(j,k,l) += + dtsdx_c2 * (F(j+1,k,l) - F(j,k,l));
#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    Ex(j,k,0) += + dtsdx_c2 * (F(j+1,k,0) - F(j  ,k,0));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ey_f_yee(int j, int k, int l,
                         amrex::Array4<amrex::Real> const& Ey,
                         amrex::Array4<amrex::Real const> const& F,
                         amrex::Real dtsdy_c2)
{
#if defined WARPX_DIM_3D
    Ey(j,k,l) += + dtsdy_c2 * (F(j,k+1,l) - F(j,k,l));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ez_f_yee(int j, int k, int l,
                         amrex::Array4<amrex::Real> const& Ez,
                         amrex::Array4<amrex::Real const> const& F,
                         amrex::Real dtsdz_c2)
{
#if defined WARPX_DIM_3D
    Ez(j,k,l) += + dtsdz_c2 * (F(j,k,l+1) - F(j,k,l));
#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    Ez(j,k,0) += + dtsdz_c2 * (F(j,k+1,0) - F(j,k,0));
#endif
}

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
void warpx_push_ex_f_ckc(int j, int k, int l,
                         amrex::Array4<amrex::Real> const& Ex,
                         amrex::Array4<amrex::Real const> const& F,
                         amrex::Real betaxy, amrex::Real betaxz, amrex::Real betayx,
                         amrex::Real betayz, amrex::Real betazx, amrex::Real betazy,
                         amrex::Real gammax, amrex::Real gammay, amrex::Real gammaz,
                         amrex::Real alphax, amrex::Real alphay, amrex::Real alphaz)
{
#if defined WARPX_DIM_3D
    Ex(j,k,l) += + alphax * (F(j+1,k  ,l  ) - F(j,  k,  l  ))
                 + betaxy * (F(j+1,k+1,l  ) - F(j,k+1,l  )
                          +  F(j+1,k-1,l  ) - F(j,k-1,l  ))
                 + betaxz * (F(j+1,k  ,l+1) - F(j,k  ,l+1)
                          +  F(j+1,k  ,l-1) - F(j,k  ,l-1))
                 + gammax * (F(j+1,k+1,l+1) - F(j,k+1,l+1)
                          +  F(j+1,k-1,l+1) - F(j,k-1,l+1)
                          +  F(j+1,k+1,l-1) - F(j,k+1,l-1)
                          +  F(j+1,k-1,l-1) - F(j,k-1,l-1));
#elif defined WARPX_DIM_XZ
    Ex(j,k,0) += + alphax * (F(j+1,k  ,0) - F(j,k  ,0))
                 + betaxz * (F(j+1,k+1,0) - F(j,k+1,0)
                          +  F(j+1,k-1,0) - F(j,k-1,0));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ey_f_ckc(int j, int k, int l,
                         amrex::Array4<amrex::Real> const& Ey,
                         amrex::Array4<amrex::Real const> const& F,
                         amrex::Real betaxy, amrex::Real betaxz, amrex::Real betayx,
                         amrex::Real betayz, amrex::Real betazx, amrex::Real betazy,
                         amrex::Real gammax, amrex::Real gammay, amrex::Real gammaz,
                         amrex::Real alphax, amrex::Real alphay, amrex::Real alphaz)
{
#if defined WARPX_DIM_3D
    Ey(j,k,l) += + alphay * (F(j  ,k+1,l  ) - F(j  ,k,l  ))
                 + betayx * (F(j+1,k+1,l  ) - F(j+1,k,l  )
                          +  F(j-1,k+1,l  ) - F(j-1,k,l  ))
                 + betayz * (F(j  ,k+1,l+1) - F(j  ,k,l+1)
                          +  F(j  ,k+1,l-1) - F(j  ,k,l-1))
                 + gammay * (F(j+1,k+1,l+1) - F(j+1,k,l+1)
                          +  F(j-1,k+1,l+1) - F(j-1,k,l+1)
                          +  F(j+1,k+1,l-1) - F(j+1,k,l-1)
                          +  F(j-1,k+1,l-1) - F(j-1,k,l-1));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ez_f_ckc(int j, int k, int l,
                         amrex::Array4<amrex::Real> const& Ez,
                         amrex::Array4<amrex::Real const> const& F,
                         amrex::Real betaxy, amrex::Real betaxz, amrex::Real betayx,
                         amrex::Real betayz, amrex::Real betazx, amrex::Real betazy,
                         amrex::Real gammax, amrex::Real gammay, amrex::Real gammaz,
                         amrex::Real alphax, amrex::Real alphay, amrex::Real alphaz)
{
#if defined WARPX_DIM_3D
    Ez(j,k,l) += + alphaz * (F(j  ,k  ,l+1) - F(j,  k,  l))
                 + betazx * (F(j+1,k  ,l+1) - F(j+1,k  ,l)
                          +  F(j-1,k  ,l+1) - F(j-1,k  ,l))
                 + betazy * (F(j  ,k+1,l+1) - F(j  ,k+1,l)
                          +  F(j  ,k-1,l+1) - F(j  ,k-1,l))
                 + gammaz * (F(j+1,k+1,l+1) - F(j+1,k+1,l)
                          +  F(j-1,k+1,l+1) - F(j-1,k+1,l)
                          +  F(j+1,k-1,l+1) - F(j+1,k-1,l)
                          +  F(j-1,k-1,l+1) - F(j-1,k-1,l));
#elif defined WARPX_DIM_XZ
    Ez(j,k,0) += + alphaz * (F(j  ,k+1,0) - F(j  ,k,0))
                 + betazx * (F(j+1,k+1,0) - F(j+1,k,0)
                          +  F(j-1,k+1,0) - F(j-1,k,0));
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_computedive(int i, int j, int k, int dcomp,
                       amrex::Array4<amrex::Real> const& divE,
                       amrex::Array4<amrex::Real const> const& Ex,
                       amrex::Array4<amrex::Real const> const& Ey,
                       amrex::Array4<amrex::Real const> const& Ez,
                       amrex::Real dxinv, amrex::Real dyinv,
                       amrex::Real dzinv
#ifdef WARPX_DIM_RZ
                       , amrex::Real const rmin
#endif
                       )
{
#if defined WARPX_DIM_3D
    divE(i,j,k,dcomp) = (Ex(i,j,k) - Ex(i-1,j,k))*dxinv
               +        (Ey(i,j,k) - Ey(i,j-1,k))*dyinv
               +        (Ez(i,j,k) - Ez(i,j,k-1))*dzinv;
#elif defined WARPX_DIM_XZ
    divE(i,j,0,dcomp) = (Ex(i,j,0) - Ex(i-1,j,0))*dxinv
                 +      (Ez(i,j,0) - Ez(i,j-1,0))*dzinv;
#elif defined WARPX_DIM_RZ
    if (i == 0 && rmin == 0.) {
        // the bulk equation diverges on axis
        // (due to the 1/r terms). The following expressions regularize
        // these divergences.
        divE(i,j,0,dcomp) = 4.*Ex(i,j,0)*dxinv
                            + (Ez(i,j,0) - Ez(i,j-1,0))*dzinv;
    } else {
        const amrex::Real ru = 1. + 0.5/(rmin*dxinv + i);
        const amrex::Real rd = 1. - 0.5/(rmin*dxinv + i);
        divE(i,j,0,dcomp) = (ru*Ex(i,j,0) - rd*Ex(i-1,j,0))*dxinv
                           + (Ez(i,j,0) - Ez(i,j-1,0))*dzinv;
    }
#endif
}

#endif
