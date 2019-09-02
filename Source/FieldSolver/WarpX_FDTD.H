#ifndef WARPX_FDTD_H_
#define WARPX_FDTD_H_

#include <AMReX_FArrayBox.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_bx_yee(int i, int j, int k, Array4<Real> const& Bx,
                       Array4<Real const> const& Ey, Array4<Real const> const& Ez,
                       Real dtsdy, Real dtsdz)
{
#if defined WARPX_DIM_3D
    Bx(i,j,k) += - dtsdy * (Ez(i,j+1,k  ) - Ez(i,j,k))
                 + dtsdz * (Ey(i,j  ,k+1) - Ey(i,j,k));
#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    // Note that the 2D Cartesian and RZ versions are the same
    Bx(i,j,0) += + dtsdz * (Ey(i,j+1,0) - Ey(i,j,0));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_by_yee(int i, int j, int k, Array4<Real> const& By,
                       Array4<Real const> const& Ex, Array4<Real const> const& Ez,
                       Real dtsdx, Real dtsdz)
{
#if defined WARPX_DIM_3D
    By(i,j,k) += + dtsdx * (Ez(i+1,j,k  ) - Ez(i,j,k))
                 - dtsdz * (Ex(i  ,j,k+1) - Ex(i,j,k));
#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    // Note that the 2D Cartesian and RZ versions are the same
    By(i,j,0) += + dtsdx * (Ez(i+1,j  ,0) - Ez(i,j,0))
                 - dtsdz * (Ex(i  ,j+1,0) - Ex(i,j,0));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_bz_yee(int i, int j, int k, Array4<Real> const& Bz,
                       Array4<Real const> const& Ex, Array4<Real const> const& Ey,
                       Real dtsdx, Real dtsdy, Real dxinv, Real rmin)
{
#if defined WARPX_DIM_3D
    Bz(i,j,k) += - dtsdx * (Ey(i+1,j  ,k) - Ey(i,j,k))
                 + dtsdy * (Ex(i  ,j+1,k) - Ex(i,j,k));
#elif defined WARPX_DIM_XZ
    Bz(i,j,0) += - dtsdx * (Ey(i+1,j,0) - Ey(i,j,0));
#elif defined WARPX_DIM_RZ
    const Real ru = 1. + 0.5/(rmin*dxinv + i + 0.5);
    const Real rd = 1. - 0.5/(rmin*dxinv + i + 0.5);
    Bz(i,j,0) += - dtsdx * (ru*Ey(i+1,j,0) - rd*Ey(i,j,0));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ex_yee(int i, int j, int k, Array4<Real> const& Ex,
                       Array4<Real const> const& By, Array4<Real const> const& Bz, Array4<Real const> const& Jx,
                       Real mu_c2_dt, Real dtsdy_c2, Real dtsdz_c2)
{
#if defined WARPX_DIM_3D
    Ex(i,j,k) += + dtsdy_c2 * (Bz(i,j,k) - Bz(i,j-1,k  ))
                 - dtsdz_c2 * (By(i,j,k) - By(i,j  ,k-1))
                 - mu_c2_dt  * Jx(i,j,k);
#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    // Note that the 2D Cartesian and RZ versions are the same
    Ex(i,j,0) += - dtsdz_c2 * (By(i,j,0) - By(i,j-1,0))
                 - mu_c2_dt  * Jx(i,j,0);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ey_yee(int i, int j, int k, Array4<Real> const& Ey,
                       Array4<Real const> const& Bx, Array4<Real const> const& Bz, Array4<Real const> const& Jy,
                       Real mu_c2_dt, Real dtsdx_c2, Real dtsdz_c2, Real rmin)
{
#if defined WARPX_DIM_3D
    Ey(i,j,k) += - dtsdx_c2 * (Bz(i,j,k) - Bz(i-1,j,k))
                 + dtsdz_c2 * (Bx(i,j,k) - Bx(i,j,k-1))
                 - mu_c2_dt  * Jy(i,j,k);
#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    // 2D Cartesian and RZ are the same, except that the axis is skipped with RZ
#ifdef WARPX_DIM_RZ
    if (i != 0 || rmin != 0.)
#endif
    {
    Ey(i,j,0) += - dtsdx_c2 * (Bz(i,j,0) - Bz(i-1,j,0))
                 + dtsdz_c2 * (Bx(i,j,0) - Bx(i,j-1,0))
                 - mu_c2_dt  * Jy(i,j,0);
    }
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ez_yee(int i, int j, int k, Array4<Real> const& Ez,
                       Array4<Real const> const& Bx, Array4<Real const> const& By, Array4<Real const> const& Jz,
                       Real mu_c2_dt, Real dtsdx_c2, Real dtsdy_c2, Real dxinv, Real rmin)
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
        const Real ru = 1. + 0.5/(rmin*dxinv + i);
        const Real rd = 1. - 0.5/(rmin*dxinv + i);
        Ez(i,j,0) += + dtsdx_c2 * (ru*By(i,j,0) - rd*By(i-1,j,0))
                     - mu_c2_dt  * Jz(i,j,0);
    } else {
        Ez(i,j,0) += + 4.*dtsdx_c2 * By(i,j,0)
                       - mu_c2_dt  * Jz(i,j,0);
    }
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ex_f_yee(int j, int k, int l, Array4<Real> const& Ex,
                         Array4<Real const> const& F, Real dtsdx_c2)
{
#if defined WARPX_DIM_3D
    Ex(j,k,l) += + dtsdx_c2 * (F(j+1,k,l) - F(j,k,l));
#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    Ex(j,k,0) += + dtsdx_c2 * (F(j+1,k,0) - F(j  ,k,0));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ey_f_yee(int j, int k, int l, Array4<Real> const& Ey,
                         Array4<Real const> const& F, Real dtsdy_c2)
{
#if defined WARPX_DIM_3D
    Ey(j,k,l) += + dtsdy_c2 * (F(j,k+1,l) - F(j,k,l));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ez_f_yee(int j, int k, int l, Array4<Real> const& Ez,
                         Array4<Real const> const& F, Real dtsdz_c2)
{
#if defined WARPX_DIM_3D
    Ez(j,k,l) += + dtsdz_c2 * (F(j,k,l+1) - F(j,k,l));
#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    Ez(j,k,0) += + dtsdz_c2 * (F(j,k+1,0) - F(j,k,0));
#endif
}

static void warpx_calculate_ckc_coefficients(Real dtsdx, Real dtsdy, Real dtsdz,
                                             Real &betaxy, Real &betaxz, Real &betayx, Real &betayz, Real &betazx, Real &betazy,
                                             Real &gammax, Real &gammay, Real &gammaz,
                                             Real &alphax, Real &alphay, Real &alphaz)
{
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
void warpx_push_bx_ckc(int j, int k, int l, Array4<Real> const& Bx,
                       Array4<Real const> const& Ey, Array4<Real const> const& Ez,
                       Real betaxy, Real betaxz, Real betayx, Real betayz, Real betazx, Real betazy,
                       Real gammax, Real gammay, Real gammaz,
                       Real alphax, Real alphay, Real alphaz)
{
#if defined WARPX_DIM_3D
    Bx(j,k,l) += - alphay * (Ez(j  ,k+1,l  ) - Ez(j,  k  ,l  ))
                 - betayx * (Ez(j+1,k+1,l  ) - Ez(j+1,k  ,l  )
                          +  Ez(j-1,k+1,l  ) - Ez(j-1,k  ,l  ))
                 - betayz * (Ez(j  ,k+1,l+1) - Ez(j  ,k  ,l+1)
                          +  Ez(j  ,k+1,l-1) - Ez(j  ,k  ,l-1))
                 - gammay * (Ez(j+1,k+1,l+1) - Ez(j+1,k  ,l+1)
                          +  Ez(j-1,k+1,l+1) - Ez(j-1,k  ,l+1)
                          +  Ez(j+1,k+1,l-1) - Ez(j+1,k  ,l-1)
                          +  Ez(j-1,k+1,l-1) - Ez(j-1,k  ,l-1))
                 + alphaz * (Ey(j  ,k  ,l+1) - Ey(j,  k,  l  ))
                 + betazx * (Ey(j+1,k  ,l+1) - Ey(j+1,k  ,l  )
                          +  Ey(j-1,k  ,l+1) - Ey(j-1,k  ,l  ))
                 + betazy * (Ey(j  ,k+1,l+1) - Ey(j  ,k+1,l  )
                          +  Ey(j  ,k-1,l+1) - Ey(j  ,k-1,l  ))
                 + gammaz * (Ey(j+1,k+1,l+1) - Ey(j+1,k+1,l  )
                          +  Ey(j-1,k+1,l+1) - Ey(j-1,k+1,l  )
                          +  Ey(j+1,k-1,l+1) - Ey(j+1,k-1,l  )
                          +  Ey(j-1,k-1,l+1) - Ey(j-1,k-1,l  ));
#elif defined WARPX_DIM_XZ
    Bx(j,k,0) += + alphaz * (Ey(j  ,k+1,0) - Ey(j,  k,0))
                 + betazx * (Ey(j+1,k+1,0) - Ey(j+1,k,0)
                          +  Ey(j-1,k+1,0) - Ey(j-1,k,0));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_by_ckc(int j, int k, int l, Array4<Real> const& By,
                       Array4<Real const> const& Ex, Array4<Real const> const& Ez,
                       Real betaxy, Real betaxz, Real betayx, Real betayz, Real betazx, Real betazy,
                       Real gammax, Real gammay, Real gammaz,
                       Real alphax, Real alphay, Real alphaz)
{
#if defined WARPX_DIM_3D
    By(j,k,l) += + alphax * (Ez(j+1,k  ,l  ) - Ez(j,  k,  l  ))
                 + betaxy * (Ez(j+1,k+1,l  ) - Ez(j  ,k+1,l  )
                          +  Ez(j+1,k-1,l  ) - Ez(j  ,k-1,l  ))
                 + betaxz * (Ez(j+1,k  ,l+1) - Ez(j  ,k  ,l+1)
                          +  Ez(j+1,k  ,l-1) - Ez(j  ,k  ,l-1))
                 + gammax * (Ez(j+1,k+1,l+1) - Ez(j  ,k+1,l+1)
                          +  Ez(j+1,k-1,l+1) - Ez(j  ,k-1,l+1)
                          +  Ez(j+1,k+1,l-1) - Ez(j  ,k+1,l-1)
                          +  Ez(j+1,k-1,l-1) - Ez(j  ,k-1,l-1))
                 - alphaz * (Ex(j  ,k  ,l+1) - Ex(j  ,k  ,l  ))
                 - betazx * (Ex(j+1,k  ,l+1) - Ex(j+1,k  ,l  )
                          +  Ex(j-1,k  ,l+1) - Ex(j-1,k  ,l  ))
                 - betazy * (Ex(j  ,k+1,l+1) - Ex(j  ,k+1,l  )
                          +  Ex(j  ,k-1,l+1) - Ex(j  ,k-1,l  ))
                 - gammaz * (Ex(j+1,k+1,l+1) - Ex(j+1,k+1,l  )
                          +  Ex(j-1,k+1,l+1) - Ex(j-1,k+1,l  )
                          +  Ex(j+1,k-1,l+1) - Ex(j+1,k-1,l  )
                          +  Ex(j-1,k-1,l+1) - Ex(j-1,k-1,l  ));
#elif defined WARPX_DIM_XZ
    By(j,k,0) += + alphax * (Ez(j+1,k  ,0) - Ez(j,k  ,0))
                 + betaxz * (Ez(j+1,k+1,0) - Ez(j,k+1,0)
                          +  Ez(j+1,k-1,0) - Ez(j,k-1,0))
                 - alphaz * (Ex(j  ,k+1,0) - Ex(j,k  ,0))
                 - betazx * (Ex(j+1,k+1,0) - Ex(j+1,k,0)
                          +  Ex(j-1,k+1,0) - Ex(j-1,k,0));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_bz_ckc(int j, int k, int l, Array4<Real> const& Bz,
                       Array4<Real const> const& Ex, Array4<Real const> const& Ey,
                       Real betaxy, Real betaxz, Real betayx, Real betayz, Real betazx, Real betazy,
                       Real gammax, Real gammay, Real gammaz,
                       Real alphax, Real alphay, Real alphaz)
{
#if defined WARPX_DIM_3D
    Bz(j,k,l) += - alphax * (Ey(j+1,k  ,l  ) - Ey(j  ,k  ,l  ))
                 - betaxy * (Ey(j+1,k+1,l  ) - Ey(j  ,k+1,l  )
                          +  Ey(j+1,k-1,l  ) - Ey(j  ,k-1,l  ))
                 - betaxz * (Ey(j+1,k  ,l+1) - Ey(j  ,k  ,l+1)
                          +  Ey(j+1,k  ,l-1) - Ey(j  ,k  ,l-1))
                 - gammax * (Ey(j+1,k+1,l+1) - Ey(j  ,k+1,l+1)
                          +  Ey(j+1,k-1,l+1) - Ey(j  ,k-1,l+1)
                          +  Ey(j+1,k+1,l-1) - Ey(j  ,k+1,l-1)
                          +  Ey(j+1,k-1,l-1) - Ey(j  ,k-1,l-1))
                 + alphay * (Ex(j  ,k+1,l  ) - Ex(j  ,k  ,l  ))
                 + betayx * (Ex(j+1,k+1,l  ) - Ex(j+1,k  ,l  )
                          +  Ex(j-1,k+1,l  ) - Ex(j-1,k  ,l  ))
                 + betayz * (Ex(j  ,k+1,l+1) - Ex(j  ,k  ,l+1)
                          +  Ex(j  ,k+1,l-1) - Ex(j  ,k  ,l-1))
                 + gammay * (Ex(j+1,k+1,l+1) - Ex(j+1,k  ,l+1)
                          +  Ex(j-1,k+1,l+1) - Ex(j-1,k  ,l+1)
                          +  Ex(j+1,k+1,l-1) - Ex(j+1,k  ,l-1)
                          +  Ex(j-1,k+1,l-1) - Ex(j-1,k  ,l-1));
#elif defined WARPX_DIM_XZ
    Bz(j,k,0) += - alphax * (Ey(j+1,k  ,0) - Ey(j,k  ,0))
                 - betaxz * (Ey(j+1,k+1,0) - Ey(j,k+1,0)
                          +  Ey(j+1,k-1,0) - Ey(j,k-1,0));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_ex_f_ckc(int j, int k, int l, Array4<Real> const& Ex,
                         Array4<Real const> const& F,
                         Real betaxy, Real betaxz, Real betayx, Real betayz, Real betazx, Real betazy,
                         Real gammax, Real gammay, Real gammaz,
                         Real alphax, Real alphay, Real alphaz)
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
void warpx_push_ey_f_ckc(int j, int k, int l, Array4<Real> const& Ey,
                         Array4<Real const> const& F,
                         Real betaxy, Real betaxz, Real betayx, Real betayz, Real betazx, Real betazy,
                         Real gammax, Real gammay, Real gammaz,
                         Real alphax, Real alphay, Real alphaz)
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
void warpx_push_ez_f_ckc(int j, int k, int l, Array4<Real> const& Ez,
                         Array4<Real const> const& F,
                         Real betaxy, Real betaxz, Real betayx, Real betayz, Real betazx, Real betazy,
                         Real gammax, Real gammay, Real gammaz,
                         Real alphax, Real alphay, Real alphaz)
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
void warpx_computedivb(int i, int j, int k, int dcomp, Array4<Real> const& divB,
                       Array4<Real const> const& Bx, Array4<Real const> const& By, Array4<Real const> const& Bz,
                       Real dxinv, Real dyinv, Real dzinv
#ifdef WARPX_DIM_RZ
                       ,const Real rmin
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
    const Real ru = 1. + 0.5/(rmin*dxinv + i + 0.5);
    const Real rd = 1. - 0.5/(rmin*dxinv + i + 0.5);
    divB(i,j,0,dcomp) = (ru*Bx(i+1,j,0) - rd*Bx(i,j,0))*dxinv
                       + (Bz(i,j+1,0) - Bz(i,j,0))*dzinv;
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_computedive(int i, int j, int k, int dcomp, Array4<Real> const& divE,
                       Array4<Real const> const& Ex, Array4<Real const> const& Ey, Array4<Real const> const& Ez,
                       Real dxinv, Real dyinv, Real dzinv
#ifdef WARPX_DIM_RZ
                       ,const Real rmin
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
        const Real ru = 1. + 0.5/(rmin*dxinv + i);
        const Real rd = 1. - 0.5/(rmin*dxinv + i);
        divE(i,j,0,dcomp) = (ru*Ex(i,j,0) - rd*Ex(i-1,j,0))*dxinv
                           + (Ez(i,j,0) - Ez(i,j-1,0))*dzinv;
    }
#endif
}

#endif
