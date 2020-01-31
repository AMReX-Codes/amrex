/* Copyright 2019 Remi Lehe, Revathi Jambunathan, Revathi Jambunathan
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PML_KERNELS_H_
#define WARPX_PML_KERNELS_H_

#include <AMReX_FArrayBox.H>

using namespace amrex;


// PML BX YEE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_bx_yee (int i, int j, int k, Array4<Real> const& Bx,
                           Array4<Real const> const& Ey, Array4<Real const> const& Ez,
                           Real dtsdy, Real dtsdz)
{
#if (AMREX_SPACEDIM == 3)
    Bx(i,j,k,0) += - dtsdy * (Ez(i,j+1,k,0) + Ez(i,j+1,k,1) + Ez(i,j+1,k,2)
                   - Ez(i,j,k,0) - Ez(i,j,k,1) - Ez(i,j,k,2)) ;
    Bx(i,j,k,1) += dtsdz * (Ey(i,j,k+1,0) + Ey(i,j,k+1,1) + Ey(i,j,k+1,2)
                   - Ey(i,j,k,0) - Ey(i,j,k,1) - Ey(i,j,k,2)) ;
#else
    Bx(i,j,k,1) += dtsdz * (Ey(i,j+1,k,0) + Ey(i,j+1,k,1) + Ey(i,j+1,k,2)
                   - Ey(i,j,k,0) - Ey(i,j,k,1) - Ey(i,j,k,2)) ;
#endif

}

//PML BY YEE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_by_yee (int i, int j, int k, Array4<Real> const& By,
                           Array4<Real const> const& Ex, Array4<Real const> const& Ez,
                           Real dtsdx, Real dtsdz)
{
#if (AMREX_SPACEDIM == 3)
    By(i,j,k,0) += -dtsdz * (Ex(i,j,k+1,0) + Ex(i,j,k+1,1) + Ex(i,j,k+1,2)
                   - Ex(i,j,k,0) - Ex(i,j,k,1) - Ex(i,j,k,2));
#else
    By(i,j,k,0) += -dtsdz * (Ex(i,j+1,k,0) + Ex(i,j+1,k,1) + Ex(i,j+1,k,2)
                   - Ex(i,j,k,0) - Ex(i,j,k,1) - Ex(i,j,k,2));
#endif
    By(i,j,k,1) += dtsdx * (Ez(i+1,j,k,0) + Ez(i+1,j,k,1) + Ez(i+1,j,k,2)
                   - Ez(i,j,k,0) - Ez(i,j,k,1) - Ez(i,j,k,2));
}


//PML BZ YEE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_bz_yee (int i, int j, int k, Array4<Real> const& Bz,
                           Array4<Real const> const& Ex, Array4<Real const> const& Ey,
                           Real dtsdx, Real dtsdy)
{
    Bz(i,j,k,0) += -dtsdx * (Ey(i+1,j,k,0) + Ey(i+1,j,k,1) + Ey(i+1,j,k,2)
                   -Ey(i,j,k,0) - Ey(i,j,k,1) - Ey(i,j,k,2) );
#if (AMREX_SPACEDIM == 3)
    Bz(i,j,k,1) += dtsdy * (Ex(i,j+1,k,0) + Ex(i,j+1,k,1) + Ex(i,j+1,k,2)
                   -Ex(i,j,k,0) - Ex(i,j,k,1) - Ex(i,j,k,2) );
#endif
}


//PML BX CKC
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_bx_ckc(int i, int j, int k, Array4<Real> const&Bx,
                           Array4<Real const> const& Ey,
                           Array4<Real const> const& Ez,
                           Real betaxy, Real betaxz, Real betayx, Real betayz,
                           Real betazx, Real betazy, Real gammax, Real gammay,
                           Real gammaz, Real alphax, Real alphay, Real alphaz)
{
#if (AMREX_SPACEDIM == 3)
    Bx(i,j,k,0) += -alphay*(Ez(i  ,j+1,k  ,0) + Ez(i  ,j+1,k  ,1) + Ez(i  ,j+1,k  ,2)
                          - Ez(i  ,j  ,k  ,0) - Ez(i  ,j  ,k  ,1) - Ez(i  ,j  ,k  ,2))
                   -betayx*(Ez(i+1,j+1,k  ,0) + Ez(i+1,j+1,k  ,1) + Ez(i+1,j+1,k  ,2)
                          - Ez(i+1,j  ,k  ,0) - Ez(i+1,j  ,k  ,1) - Ez(i+1,j  ,k  ,2)
                          + Ez(i-1,j+1,k  ,0) + Ez(i-1,j+1,k  ,1) + Ez(i-1,j+1,k  ,2)
                          - Ez(i-1,j  ,k  ,0) - Ez(i-1,j  ,k  ,1) - Ez(i-1,j  ,k  ,2))
                   -betayz*(Ez(i  ,j+1,k+1,0) + Ez(i  ,j+1,k+1,1) + Ez(i  ,j+1,k+1,2)
                          - Ez(i  ,j  ,k+1,0) - Ez(i  ,j  ,k+1,1) - Ez(i  ,j  ,k+1,2)
                          + Ez(i  ,j+1,k-1,0) + Ez(i  ,j+1,k-1,1) + Ez(i  ,j+1,k-1,2)
                          - Ez(i  ,j  ,k-1,0) - Ez(i  ,j  ,k-1,1) - Ez(i  ,j  ,k-1,2))
                   -gammay*(Ez(i+1,j+1,k+1,0) + Ez(i+1,j+1,k+1,1) + Ez(i+1,j+1,k+1,2)
                          - Ez(i+1,j  ,k+1,0) - Ez(i+1,j  ,k+1,1) - Ez(i+1,j  ,k+1,2)
                          + Ez(i-1,j+1,k+1,0) + Ez(i-1,j+1,k+1,1) + Ez(i-1,j+1,k+1,2)
                          - Ez(i-1,j  ,k+1,0) - Ez(i-1,j  ,k+1,1) - Ez(i-1,j  ,k+1,2)
                          + Ez(i+1,j+1,k-1,0) + Ez(i+1,j+1,k-1,1) + Ez(i+1,j+1,k-1,2)
                          - Ez(i+1,j  ,k-1,0) - Ez(i+1,j  ,k-1,1) - Ez(i+1,j  ,k-1,2)
                          + Ez(i-1,j+1,k-1,0) + Ez(i-1,j+1,k-1,1) + Ez(i-1,j+1,k-1,2)
                          - Ez(i-1,j  ,k-1,0) - Ez(i-1,j  ,k-1,1) - Ez(i-1,j  ,k-1,2));

   Bx(i,j,k,1) +=   alphaz*(Ey(i  ,j  ,k+1,0) + Ey(i  ,j  ,k+1,1) + Ey(i  ,j  ,k+1,2)
                          - Ey(i  ,j  ,k  ,0) - Ey(i  ,j  ,k  ,1) - Ey(i  ,j  ,k  ,2))
                   +betazx*(Ey(i+1,j  ,k+1,0) + Ey(i+1,j  ,k+1,1) + Ey(i+1,j  ,k+1,2)
                          - Ey(i+1,j  ,k  ,0) - Ey(i+1,j  ,k  ,1) - Ey(i+1,j  ,k  ,2)
                          + Ey(i-1,j  ,k+1,0) + Ey(i-1,j  ,k+1,1) + Ey(i-1,j  ,k+1,2)
                          - Ey(i-1,j  ,k  ,0) - Ey(i-1,j  ,k  ,1) - Ey(i-1,j  ,k  ,2))
                   +betazy*(Ey(i  ,j+1,k+1,0) + Ey(i  ,j+1,k+1,1) + Ey(i  ,j+1,k+1,2)
                          - Ey(i  ,j+1,k  ,0) - Ey(i  ,j+1,k  ,1) - Ey(i  ,j+1,k  ,2)
                          + Ey(i,  j-1,k+1,0) + Ey(i  ,j-1,k+1,1) + Ey(i  ,j-1,k+1,2)
                          - Ey(i  ,j-1,k  ,0) - Ey(i  ,j-1,k  ,1) - Ey(i  ,j-1,k  ,2))
                   +gammaz*(Ey(i+1,j+1,k+1,0) + Ey(i+1,j+1,k+1,1) + Ey(i+1,j+1,k+1,2)
                          - Ey(i+1,j+1,k  ,0) - Ey(i+1,j+1,k  ,1) - Ey(i+1,j+1,k  ,2)
                          + Ey(i-1,j+1,k+1,0) + Ey(i-1,j+1,k+1,1) + Ey(i-1,j+1,k+1,2)
                          - Ey(i-1,j+1,k  ,0) - Ey(i-1,j+1,k  ,1) - Ey(i-1,j+1,k  ,2)
                          + Ey(i+1,j-1,k+1,0) + Ey(i+1,j-1,k+1,1) + Ey(i+1,j-1,k+1,2)
                          - Ey(i+1,j-1,k  ,0) - Ey(i+1,j-1,k  ,1) - Ey(i+1,j-1,k  ,2)
                          + Ey(i-1,j-1,k+1,0) + Ey(i-1,j-1,k+1,1) + Ey(i-1,j-1,k+1,2)
                          - Ey(i-1,j-1,k  ,0) - Ey(i-1,j-1,k  ,1) - Ey(i-1,j-1,k  ,2));

#else

   Bx(i,j,k,1) +=   alphaz*(Ey(i  ,j+1,k  ,0) + Ey(i  ,j+1,k  ,1) + Ey(i  ,j+1,k  ,2)
                          - Ey(i  ,j  ,k  ,0) - Ey(i  ,j  ,k  ,1) - Ey(i  ,j  ,k  ,2))
                   +betazx*(Ey(i+1,j+1,k  ,0) + Ey(i+1,j+1,k  ,1) + Ey(i+1,j+1,k  ,2)
                          - Ey(i+1,j  ,k  ,0) - Ey(i+1,j  ,k  ,1) - Ey(i+1,j  ,k  ,2)
                          + Ey(i-1,j+1,k  ,0) + Ey(i-1,j+1,k  ,1) + Ey(i-1,j+1,k  ,2)
                          - Ey(i-1,j  ,k  ,0) - Ey(i-1,j  ,k  ,1) - Ey(i-1,j  ,k  ,2));


#endif
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_by_ckc(int i, int j, int k, Array4<Real> const&By,
                           Array4<Real const> const& Ex,
                           Array4<Real const> const& Ez,
                           Real betaxy, Real betaxz, Real betayx, Real betayz,
                           Real betazx, Real betazy, Real gammax, Real gammay,
                           Real gammaz, Real alphax, Real alphay, Real alphaz)
{
#if (AMREX_SPACEDIM==3)
   By(i,j,k,0) += -alphaz*(Ex(i  ,j  ,k+1,0) + Ex(i  ,j  ,k+1,1) + Ex(i  ,j  ,k+1,2)
                         - Ex(i  ,j  ,k  ,0) - Ex(i  ,j  ,k  ,1) - Ex(i  ,j  ,k  ,2))
                  -betazx*(Ex(i+1,j  ,k+1,0) + Ex(i+1,j  ,k+1,1) + Ex(i+1,j  ,k+1,2)
                         - Ex(i+1,j  ,k  ,0) - Ex(i+1,j  ,k  ,1) - Ex(i+1,j  ,k  ,2)
                         + Ex(i-1,j  ,k+1,0) + Ex(i-1,j  ,k+1,1) + Ex(i-1,j  ,k+1,2)
                         - Ex(i-1,j  ,k  ,0) - Ex(i-1,j  ,k  ,1) - Ex(i-1,j  ,k  ,2))
                  -betazy*(Ex(i  ,j+1,k+1,0) + Ex(i  ,j+1,k+1,1) + Ex(i  ,j+1,k+1,2)
                         - Ex(i  ,j+1,k  ,0) - Ex(i  ,j+1,k  ,1) - Ex(i  ,j+1,k  ,2)
                         + Ex(i  ,j-1,k+1,0) + Ex(i  ,j-1,k+1,1) + Ex(i  ,j-1,k+1,2)
                         - Ex(i  ,j-1,k  ,0) - Ex(i  ,j-1,k  ,1) - Ex(i  ,j-1,k  ,2))
                  -gammaz*(Ex(i+1,j+1,k+1,0) + Ex(i+1,j+1,k+1,1) + Ex(i+1,j+1,k+1,2)
                         - Ex(i+1,j+1,k  ,0) - Ex(i+1,j+1,k  ,1) - Ex(i+1,j+1,k  ,2)
                         + Ex(i-1,j+1,k+1,0) + Ex(i-1,j+1,k+1,1) + Ex(i-1,j+1,k+1,2)
                         - Ex(i-1,j+1,k  ,0) - Ex(i-1,j+1,k  ,1) - Ex(i-1,j+1,k  ,2)
                         + Ex(i+1,j-1,k+1,0) + Ex(i+1,j-1,k+1,1) + Ex(i+1,j-1,k+1,2)
                         - Ex(i+1,j-1,k  ,0) - Ex(i+1,j-1,k  ,1) - Ex(i+1,j-1,k  ,2)
                         + Ex(i-1,j-1,k+1,0) + Ex(i-1,j-1,k+1,1) + Ex(i-1,j-1,k+1,2)
                         - Ex(i-1,j-1,k  ,0) - Ex(i-1,j-1,k  ,1) - Ex(i-1,j-1,k  ,2));

   By(i,j,k,1) +=  alphax*(Ez(i+1,j  ,k  ,0) + Ez(i+1,j  ,k  ,1) + Ez(i+1,j  ,k  ,2)
                         - Ez(i  ,j  ,k  ,0) - Ez(i  ,j  ,k  ,1) - Ez(i  ,j  ,k  ,2))
                  +betaxy*(Ez(i+1,j+1,k  ,0) + Ez(i+1,j+1,k  ,1) + Ez(i+1,j+1,k  ,2)
                         - Ez(i  ,j+1,k  ,0) - Ez(i  ,j+1,k  ,1) - Ez(i  ,j+1,k  ,2)
                         + Ez(i+1,j-1,k  ,0) + Ez(i+1,j-1,k  ,1) + Ez(i+1,j-1,k  ,2)
                         - Ez(i  ,j-1,k  ,0) - Ez(i  ,j-1,k  ,1) - Ez(i  ,j-1,k  ,2))
                  +betaxz*(Ez(i+1,j  ,k+1,0) + Ez(i+1,j  ,k+1,1) + Ez(i+1,j  ,k+1,2)
                         - Ez(i  ,j  ,k+1,0) - Ez(i  ,j  ,k+1,1) - Ez(i  ,j  ,k+1,2)
                         + Ez(i+1,j  ,k-1,0) + Ez(i+1,j  ,k-1,1) + Ez(i+1,j  ,k-1,2)
                         - Ez(i  ,j  ,k-1,0) - Ez(i  ,j  ,k-1,1) - Ez(i  ,j  ,k-1,2))
                  +gammax*(Ez(i+1,j+1,k+1,0) + Ez(i+1,j+1,k+1,1) + Ez(i+1,j+1,k+1,2)
                         - Ez(i  ,j+1,k+1,0) - Ez(i  ,j+1,k+1,1) - Ez(i  ,j+1,k+1,2)
                         + Ez(i+1,j-1,k+1,0) + Ez(i+1,j-1,k+1,1) + Ez(i+1,j-1,k+1,2)
                         - Ez(i  ,j-1,k+1,0) - Ez(i  ,j-1,k+1,1) - Ez(i  ,j-1,k+1,2)
                         + Ez(i+1,j+1,k-1,0) + Ez(i+1,j+1,k-1,1) + Ez(i+1,j+1,k-1,2)
                         - Ez(i  ,j+1,k-1,0) - Ez(i  ,j+1,k-1,1) - Ez(i  ,j+1,k-1,2)
                         + Ez(i+1,j-1,k-1,0) + Ez(i+1,j-1,k-1,1) + Ez(i+1,j-1,k-1,2)
                         - Ez(i  ,j-1,k-1,0) - Ez(i  ,j-1,k-1,1) - Ez(i  ,j-1,k-1,2));
#else
   By(i,j,k,0) += -alphaz*(Ex(i  ,j+1,k  ,0) + Ex(i  ,j+1,k  ,1) + Ex(i  ,j+1,k  ,2)
                         - Ex(i  ,j  ,k  ,0) - Ex(i  ,j  ,k  ,1) - Ex(i  ,j  ,k  ,2))
                  -betazx*(Ex(i+1,j+1,k  ,0) + Ex(i+1,j+1,k  ,1) + Ex(i+1,j+1,k  ,2)
                         - Ex(i+1,j  ,k  ,0) - Ex(i+1,j  ,k  ,1) - Ex(i+1,j  ,k  ,2)
                         + Ex(i-1,j+1,k  ,0) + Ex(i-1,j+1,k  ,1) + Ex(i-1,j+1,k  ,2)
                         - Ex(i-1,j  ,k  ,0) - Ex(i-1,j  ,k  ,1) - Ex(i-1,j  ,k  ,2));

   By(i,j,k,1) +=  alphax*(Ez(i+1,j  ,k  ,0) + Ez(i+1,j  ,k  ,1) + Ez(i+1,j  ,k  ,2)
                         - Ez(i  ,j  ,k  ,0) - Ez(i  ,j  ,k  ,1) - Ez(i  ,j  ,k  ,2))
                  +betaxz*(Ez(i+1,j+1,k  ,0) + Ez(i+1,j+1,k  ,1) + Ez(i+1,j+1,k  ,2)
                         - Ez(i  ,j+1,k  ,0) - Ez(i  ,j+1,k  ,1) - Ez(i  ,j+1,k  ,2)
                         + Ez(i+1,j-1,k  ,0) + Ez(i+1,j-1,k  ,1) + Ez(i+1,j-1,k  ,2)
                         - Ez(i  ,j-1,k  ,0) - Ez(i  ,j-1,k  ,1) - Ez(i  ,j-1,k  ,2));


#endif

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_bz_ckc(int i, int j, int k, Array4<Real> const&Bz,
                           Array4<Real const> const& Ex,
                           Array4<Real const> const& Ey,
                           Real betaxy, Real betaxz, Real betayx, Real betayz,
                           Real betazx, Real betazy, Real gammax, Real gammay,
                           Real gammaz, Real alphax, Real alphay, Real alphaz)
{
#if (AMREX_SPACEDIM==3)
   Bz(i,j,k,0) += -alphax*(Ey(i+1,j  ,k  ,0) + Ey(i+1,j  ,k  ,1) + Ey(i+1,j  ,k  ,2)
                         - Ey(i  ,j  ,k  ,0) - Ey(i  ,j  ,k  ,1) - Ey(i  ,j  ,k  ,2))
                  -betaxy*(Ey(i+1,j+1,k  ,0) + Ey(i+1,j+1,k  ,1) + Ey(i+1,j+1,k  ,2)
                         - Ey(i  ,j+1,k  ,0) - Ey(i  ,j+1,k  ,1) - Ey(i  ,j+1,k  ,2)
                         + Ey(i+1,j-1,k  ,0) + Ey(i+1,j-1,k  ,1) + Ey(i+1,j-1,k  ,2)
                         - Ey(i  ,j-1,k  ,0) - Ey(i  ,j-1,k  ,1) - Ey(i  ,j-1,k  ,2))
                  -betaxz*(Ey(i+1,j  ,k+1,0) + Ey(i+1,j  ,k+1,1) + Ey(i+1,j  ,k+1,2)
                         - Ey(i  ,j  ,k+1,0) - Ey(i  ,j  ,k+1,1) - Ey(i  ,j  ,k+1,2)
                         + Ey(i+1,j  ,k-1,0) + Ey(i+1,j  ,k-1,1) + Ey(i+1,j  ,k-1,2)
                         - Ey(i  ,j  ,k-1,0) - Ey(i  ,j  ,k-1,1) - Ey(i  ,j  ,k-1,2))
                  -gammax*(Ey(i+1,j+1,k+1,0) + Ey(i+1,j+1,k+1,1) + Ey(i+1,j+1,k+1,2)
                         - Ey(i  ,j+1,k+1,0) - Ey(i  ,j+1,k+1,1) - Ey(i  ,j+1,k+1,2)
                         + Ey(i+1,j-1,k+1,0) + Ey(i+1,j-1,k+1,1) + Ey(i+1,j-1,k+1,2)
                         - Ey(i  ,j-1,k+1,0) - Ey(i  ,j-1,k+1,1) - Ey(i  ,j-1,k+1,2)
                         + Ey(i+1,j+1,k-1,0) + Ey(i+1,j+1,k-1,1) + Ey(i+1,j+1,k-1,2)
                         - Ey(i  ,j+1,k-1,0) - Ey(i  ,j+1,k-1,1) - Ey(i  ,j+1,k-1,2)
                         + Ey(i+1,j-1,k-1,0) + Ey(i+1,j-1,k-1,1) + Ey(i+1,j-1,k-1,2)
                         - Ey(i  ,j-1,k-1,0) - Ey(i  ,j-1,k-1,1) - Ey(i  ,j-1,k-1,2));

   Bz(i,j,k,1) +=  alphay*(Ex(i  ,j+1,k  ,0) + Ex(i  ,j+1,k  ,1) + Ex(i  ,j+1,k  ,2)
                         - Ex(i  ,j  ,k  ,0) - Ex(i  ,j  ,k  ,1) - Ex(i  ,j  ,k  ,2))
                  +betayx*(Ex(i+1,j+1,k  ,0) + Ex(i+1,j+1,k  ,1) + Ex(i+1,j+1,k  ,2)
                         - Ex(i+1,j  ,k  ,0) - Ex(i+1,j  ,k  ,1) - Ex(i+1,j  ,k  ,2)
                         + Ex(i-1,j+1,k  ,0) + Ex(i-1,j+1,k  ,1) + Ex(i-1,j+1,k  ,2)
                         - Ex(i-1,j  ,k  ,0) - Ex(i-1,j  ,k  ,1) - Ex(i-1,j  ,k  ,2))
                  +betayz*(Ex(i  ,j+1,k+1,0) + Ex(i  ,j+1,k+1,1) + Ex(i  ,j+1,k+1,2)
                         - Ex(i  ,j  ,k+1,0) - Ex(i  ,j  ,k+1,1) - Ex(i  ,j  ,k+1,2)
                         + Ex(i  ,j+1,k-1,0) + Ex(i  ,j+1,k-1,1) + Ex(i  ,j+1,k-1,2)
                         - Ex(i  ,j  ,k-1,0) - Ex(i  ,j  ,k-1,1) - Ex(i  ,j  ,k-1,2))
                  +gammay*(Ex(i+1,j+1,k+1,0) + Ex(i+1,j+1,k+1,1) + Ex(i+1,j+1,k+1,2)
                         - Ex(i+1,j  ,k+1,0) - Ex(i+1,j  ,k+1,1) - Ex(i+1,j  ,k+1,2)
                         + Ex(i-1,j+1,k+1,0) + Ex(i-1,j+1,k+1,1) + Ex(i-1,j+1,k+1,2)
                         - Ex(i-1,j  ,k+1,0) - Ex(i-1,j  ,k+1,1) - Ex(i-1,j  ,k+1,2)
                         + Ex(i+1,j+1,k-1,0) + Ex(i+1,j+1,k-1,1) + Ex(i+1,j+1,k-1,2)
                         - Ex(i+1,j  ,k-1,0) - Ex(i+1,j  ,k-1,1) - Ex(i+1,j  ,k-1,2)
                         + Ex(i-1,j+1,k-1,0) + Ex(i-1,j+1,k-1,1) + Ex(i-1,j+1,k-1,2)
                         - Ex(i-1,j  ,k-1,0) - Ex(i-1,j  ,k-1,1) - Ex(i-1,j  ,k-1,2));


#else
   Bz(i,j,k,0) += -alphax*(Ey(i+1,j  ,k  ,0) + Ey(i+1,j  ,k  ,1) + Ey(i+1,j  ,k  ,2)
                         - Ey(i  ,j  ,k  ,0) - Ey(i  ,j  ,k  ,1) - Ey(i  ,j  ,k  ,2))
                  -betaxz*(Ey(i+1,j+1,k  ,0) + Ey(i+1,j+1,k  ,1) + Ey(i+1,j+1,k  ,2)
                         - Ey(i  ,j+1,k  ,0) - Ey(i  ,j+1,k  ,1) - Ey(i  ,j+1,k  ,2)
                         + Ey(i+1,j-1,k  ,0) + Ey(i+1,j-1,k  ,1) + Ey(i+1,j-1,k  ,2)
                         - Ey(i  ,j-1,k  ,0) - Ey(i  ,j-1,k  ,1) - Ey(i  ,j-1,k  ,2));
#endif

}





//PML EX YEE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_ex_yee (int i, int j, int k, Array4<Real> const& Ex,
                           Array4<Real const> const& By,
                           Array4<Real const> const& Bz,
                           Real dtsdy_c2, Real dtsdz_c2)
{
#if (AMREX_SPACEDIM == 3)
    Ex(i,j,k,0) += dtsdy_c2 * ( Bz(i,j  ,k,0) + Bz(i,j  ,k,1)
                            -   Bz(i,j-1,k,0) - Bz(i,j-1,k,1) );
    Ex(i,j,k,1) += -dtsdz_c2 * (By(i,j,k  ,0) + By(i,j,k  ,1)
                             -  By(i,j,k-1,0) - By(i,j,k-1,1) );

#else
    Ex(i,j,k,1) += -dtsdz_c2 * (By(i,j  ,k,0) + By(i,j  ,k,1)
                             -  By(i,j-1,k,0) - By(i,j-1,k,1) );
#endif
}

//PML EY YEE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_ey_yee (int i, int j, int k, Array4<Real> const& Ey,
                            Array4<Real const> const& Bx,
                            Array4<Real const> const& Bz,
                            Real dtsdx_c2, Real dtsdz_c2)
{
#if (AMREX_SPACEDIM == 3)
    Ey(i,j,k,0) += dtsdz_c2 * (Bx(i,j,k  ,0) + Bx(i,j,k  ,1)
                             - Bx(i,j,k-1,0) - Bx(i,j,k-1,1));
#else
    Ey(i,j,k,0) += dtsdz_c2 * (Bx(i,j  ,k,0) + Bx(i,j  ,k,1)
                             - Bx(i,j-1,k,0) - Bx(i,j-1,k,1));
#endif
    Ey(i,j,k,1) += -dtsdx_c2 * (Bz(i  ,j,k,0) + Bz(i  ,j,k,1)
                              - Bz(i-1,j,k,0) - Bz(i-1,j,k,1) );
}

//PML EZ YEE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_ez_yee (int i, int j, int k, Array4<Real> const& Ez,
                           Array4<Real const> const& Bx,
                           Array4<Real const> const& By,
                           Real dtsdx_c2, Real dtsdy_c2)
{
    Ez(i,j,k,0) += dtsdx_c2 * (By(i,  j,k,0) + By(i,  j,k,1)
                             - By(i-1,j,k,0) - By(i-1,j,k,1) );
#if (AMREX_SPACEDIM==3)
    Ez(i,j,k,1) += - dtsdy_c2 * (Bx(i,j,  k,0) + Bx(i,j,  k,1)
                               - Bx(i,j-1,k,0) - Bx(i,j-1,k,1) );
#endif
}


// PML EX F YEE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_ex_f_yee(int i, int j, int k,
                             Array4<Real> const& Exfab,
                             Array4<Real const> const& F_fab,
                             Real dtsdx_c2)
{
    Exfab(i,j,k,2) += dtsdx_c2 * (F_fab(i+1,j,k,0) - F_fab(i,j,k,0)
                        +         F_fab(i+1,j,k,1) - F_fab(i,j,k,1)
                        +         F_fab(i+1,j,k,2) - F_fab(i,j,k,2));
}

// PML EY F YEE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_ey_f_yee(int i, int j, int k,
                             Array4<Real> const& Eyfab,
                             Array4<Real const> const& F_fab,
                             Real dtsdy_c2)
{
#if (AMREX_SPACEDIM == 3)
    Eyfab(i,j,k,2) += dtsdy_c2 * (F_fab(i,j+1,k,0) - F_fab(i,j,k,0)
                         +        F_fab(i,j+1,k,1) - F_fab(i,j,k,1)
                         +        F_fab(i,j+1,k,2) - F_fab(i,j,k,2));
#endif
}

// PML EZ F YEE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_ez_f_yee(int i, int j, int k,
                             Array4<Real> const& Ezfab,
                             Array4<Real const> const& F_fab,
                             Real dtsdz_c2)
{
#if (AMREX_SPACEDIM == 3)
   Ezfab(i,j,k,2) += dtsdz_c2 * (F_fab(i,j,k+1,0) - F_fab(i,j,k,0)
                       +         F_fab(i,j,k+1,1) - F_fab(i,j,k,1)
                       +         F_fab(i,j,k+1,2) - F_fab(i,j,k,2));
#else
   Ezfab(i,j,k,2) += dtsdz_c2 * (F_fab(i,j+1,k,0) - F_fab(i,j,k,0)
                       +         F_fab(i,j+1,k,1) - F_fab(i,j,k,1)
                       +         F_fab(i,j+1,k,2) - F_fab(i,j,k,2));
#endif
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_ex_f_ckc(int i, int j, int k,
                             Array4<Real> const&  Exfab,
                             Array4<Real const> const& F,
                             Real alphax, Real betaxy, Real betaxz,Real gammax)
{
#if (AMREX_SPACEDIM==3)
   Exfab(i,j,k,2) +=  alphax*(
                        F(i+1,j  ,k  ,0) + F(i+1,j  ,k  ,1) + F(i+1,j  ,k  ,2)
                      - F(i  ,j  ,k  ,0) - F(i  ,j  ,k  ,1) - F(i  ,j  ,k  ,2))
                     +betaxy*(
                        F(i+1,j+1,k  ,0) + F(i+1,j+1,k  ,1) + F(i+1,j+1,k  ,2)
                      - F(i  ,j+1,k  ,0) - F(i  ,j+1,k  ,1) - F(i  ,j+1,k  ,2)
                      + F(i+1,j-1,k  ,0) + F(i+1,j-1,k  ,1) + F(i+1,j-1,k  ,2)
                      - F(i  ,j-1,k  ,0) - F(i  ,j-1,k  ,1) - F(i  ,j-1,k  ,2))
                     +betaxz*(
                        F(i+1,j  ,k+1,0) + F(i+1,j  ,k+1,1) + F(i+1,j  ,k+1,2)
                      - F(i  ,j  ,k+1,0) - F(i  ,j  ,k+1,1) - F(i  ,j  ,k+1,2)
                      + F(i+1,j  ,k-1,0) + F(i+1,j  ,k-1,1) + F(i+1,j  ,k-1,2)
                      - F(i  ,j  ,k-1,0) - F(i  ,j  ,k-1,1) - F(i  ,j  ,k-1,2))
                     +gammax*(
                        F(i+1,j+1,k+1,0) + F(i+1,j+1,k+1,1) + F(i+1,j+1,k+1,2)
                      - F(i  ,j+1,k+1,0) - F(i  ,j+1,k+1,1) - F(i  ,j+1,k+1,2)
                      + F(i+1,j-1,k+1,0) + F(i+1,j-1,k+1,1) + F(i+1,j-1,k+1,2)
                      - F(i  ,j-1,k+1,0) - F(i  ,j-1,k+1,1) - F(i  ,j-1,k+1,2)
                      + F(i+1,j+1,k-1,0) + F(i+1,j+1,k-1,1) + F(i+1,j+1,k-1,2)
                      - F(i  ,j+1,k-1,0) - F(i  ,j+1,k-1,1) - F(i  ,j+1,k-1,2)
                      + F(i+1,j-1,k-1,0) + F(i+1,j-1,k-1,1) + F(i+1,j-1,k-1,2)
                      - F(i  ,j-1,k-1,0) - F(i  ,j-1,k-1,1) - F(i  ,j-1,k-1,2));
#else
   Exfab(i,j,k,2) +=  alphax*(
                        F(i+1,j  ,k  ,0) + F(i+1,j  ,k  ,1) + F(i+1,j  ,k  ,2)
                      - F(i  ,j  ,k  ,0) - F(i  ,j  ,k  ,1) - F(i  ,j  ,k  ,2))
                     +betaxz*(
                        F(i+1,j+1,k  ,0) + F(i+1,j+1,k  ,1) + F(i+1,j+1,k  ,2)
                      - F(i  ,j+1,k  ,0) - F(i  ,j+1,k  ,1) - F(i  ,j+1,k  ,2)
                      + F(i+1,j-1,k  ,0) + F(i+1,j-1,k  ,1) + F(i+1,j-1,k  ,2)
                      - F(i  ,j-1,k  ,0) - F(i  ,j-1,k  ,1) - F(i  ,j-1,k  ,2));
#endif

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_ey_f_ckc(int i, int j, int k,
                             Array4<Real> const&  Eyfab,
                             Array4<Real const> const& F,
                             Real alphay, Real betayx, Real betayz,Real gammay)
{
#if (AMREX_SPACEDIM==3)
   Eyfab(i,j,k,2) +=  alphay*(
                        F(i  ,j+1,k  ,0) + F(i  ,j+1,k  ,1) + F(i  ,j+1,k  ,2)
                      - F(i  ,j  ,k  ,0) - F(i  ,j  ,k  ,1) - F(i  ,j  ,k  ,2))
                     +betayx*(
                        F(i+1,j+1,k  ,0) + F(i+1,j+1,k  ,1) + F(i+1,j+1,k  ,2)
                      - F(i+1,j  ,k  ,0) - F(i+1,j  ,k  ,1) - F(i+1,j  ,k  ,2)
                      + F(i-1,j+1,k  ,0) + F(i-1,j+1,k  ,1) + F(i-1,j+1,k  ,2)
                      - F(i-1,j  ,k  ,0) - F(i-1,j  ,k  ,1) - F(i-1,j  ,k  ,2))
                     +betayz*(
                        F(i  ,j+1,k+1,0) + F(i  ,j+1,k+1,1) + F(i  ,j+1,k+1,2)
                      - F(i  ,j  ,k+1,0) - F(i  ,j  ,k+1,1) - F(i  ,j  ,k+1,2)
                      + F(i  ,j+1,k-1,0) + F(i  ,j+1,k-1,1) + F(i  ,j+1,k-1,2)
                      - F(i  ,j  ,k-1,0) - F(i  ,j  ,k-1,1) - F(i  ,j  ,k-1,2))
                     +gammay*(
                        F(i+1,j+1,k+1,0) + F(i+1,j+1,k+1,1) + F(i+1,j+1,k+1,2)
                      - F(i+1,j  ,k+1,0) - F(i+1,j  ,k+1,1) - F(i+1,j  ,k+1,2)
                      + F(i-1,j+1,k+1,0) + F(i-1,j+1,k+1,1) + F(i-1,j+1,k+1,2)
                      - F(i-1,j  ,k+1,0) - F(i-1,j  ,k+1,1) - F(i-1,j  ,k+1,2)
                      + F(i+1,j+1,k-1,0) + F(i+1,j+1,k-1,1) + F(i+1,j+1,k-1,2)
                      - F(i+1,j  ,k-1,0) - F(i+1,j  ,k-1,1) - F(i+1,j  ,k-1,2)
                      + F(i-1,j+1,k-1,0) + F(i-1,j+1,k-1,1) + F(i-1,j+1,k-1,2)
                      - F(i-1,j  ,k-1,0) - F(i-1,j  ,k-1,1) - F(i-1,j  ,k-1,2));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_ez_f_ckc(int i, int j, int k,
                             Array4<Real> const&  Ezfab,
                             Array4<Real const> const& F,
                             Real alphaz, Real betazx, Real betazy,Real gammaz)
{
#if (AMREX_SPACEDIM==3)
   Ezfab(i,j,k,2) +=  alphaz*(
                        F(i  ,j  ,k+1,0) + F(i  ,j  ,k+1,1) + F(i  ,j  ,k+1,2)
                      - F(i  ,j  ,k  ,0) - F(i  ,j  ,k  ,1) - F(i  ,j  ,k  ,2))
                     +betazx*(
                        F(i+1,j  ,k+1,0) + F(i+1,j  ,k+1,1) + F(i+1,j  ,k+1,2)
                      - F(i+1,j  ,k  ,0) - F(i+1,j  ,k  ,1) - F(i+1,j  ,k  ,2)
                      + F(i-1,j  ,k+1,0) + F(i-1,j  ,k+1,1) + F(i-1,j  ,k+1,2)
                      - F(i-1,j  ,k  ,0) - F(i-1,j  ,k  ,1) - F(i-1,j  ,k  ,2))
                     +betazy*(
                        F(i  ,j+1,k+1,0) + F(i  ,j+1,k+1,1) + F(i  ,j+1,k+1,2)
                      - F(i  ,j+1,k  ,0) - F(i  ,j+1,k  ,1) - F(i  ,j+1,k  ,2)
                      + F(i  ,j-1,k+1,0) + F(i  ,j-1,k+1,1) + F(i  ,j-1,k+1,2)
                      - F(i  ,j-1,k  ,0) - F(i  ,j-1,k  ,1) - F(i  ,j-1,k  ,2))
                     +gammaz*(
                        F(i+1,j+1,k+1,0) + F(i+1,j+1,k+1,1) + F(i+1,j+1,k+1,2)
                      - F(i+1,j+1,k  ,0) - F(i+1,j+1,k  ,1) - F(i+1,j+1,k  ,2)
                      + F(i-1,j+1,k+1,0) + F(i-1,j+1,k+1,1) + F(i-1,j+1,k+1,2)
                      - F(i-1,j+1,k  ,0) - F(i-1,j+1,k  ,1) - F(i-1,j+1,k  ,2)
                      + F(i+1,j-1,k+1,0) + F(i+1,j-1,k+1,1) + F(i+1,j-1,k+1,2)
                      - F(i+1,j-1,k  ,0) - F(i+1,j-1,k  ,1) - F(i+1,j-1,k  ,2)
                      + F(i-1,j-1,k+1,0) + F(i-1,j-1,k+1,1) + F(i-1,j-1,k+1,2)
                      - F(i-1,j-1,k  ,0) - F(i-1,j-1,k  ,1) - F(i-1,j-1,k  ,2));
#else
   Ezfab(i,j,k,2) +=  alphaz*(
                        F(i  ,j+1,k  ,0) + F(i  ,j+1,k  ,1) + F(i  ,j+1,k  ,2)
                      - F(i  ,j  ,k  ,0) - F(i  ,j  ,k  ,1) - F(i  ,j  ,k  ,2))
                     +betazx*(
                        F(i+1,j+1,k  ,0) + F(i+1,j+1,k  ,1) + F(i+1,j+1,k  ,2)
                      - F(i+1,j  ,k  ,0) - F(i+1,j  ,k  ,1) - F(i+1,j  ,k  ,2)
                      + F(i-1,j+1,k  ,0) + F(i-1,j+1,k  ,1) + F(i-1,j+1,k  ,2)
                      - F(i-1,j  ,k  ,0) - F(i-1,j  ,k  ,1) - F(i-1,j  ,k  ,2));
#endif

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_push_pml_F(int i, int j, int k, Array4<Real> const& F_fab,
                      Array4<Real const> const& Exfab,
                      Array4<Real const> const& Eyfab,
                      Array4<Real const> const& Ezfab,
                      Real dtsdx, Real dtsdy, Real dtsdz)
{
   F_fab(i,j,k,0) += dtsdx * (Exfab(i,j,k,0) - Exfab(i-1,j,k,0)
                            + Exfab(i,j,k,1) - Exfab(i-1,j,k,1)
                            + Exfab(i,j,k,2) - Exfab(i-1,j,k,2));
#if (AMREX_SPACEDIM==3)
   F_fab(i,j,k,1) += dtsdy * (Eyfab(i,j,k,0) - Eyfab(i,j-1,k,0)
                            + Eyfab(i,j,k,1) - Eyfab(i,j-1,k,1)
                            + Eyfab(i,j,k,2) - Eyfab(i,j-1,k,2));
   F_fab(i,j,k,2) += dtsdz * (Ezfab(i,j,k,0) - Ezfab(i,j,k-1,0)
                            + Ezfab(i,j,k,1) - Ezfab(i,j,k-1,1)
                            + Ezfab(i,j,k,2) - Ezfab(i,j,k-1,2));
#else
   F_fab(i,j,k,2) += dtsdz * (Ezfab(i,j,k,0) - Ezfab(i,j-1,k,0)
                            + Ezfab(i,j,k,1) - Ezfab(i,j-1,k,1)
                            + Ezfab(i,j,k,2) - Ezfab(i,j-1,k,2));
#endif
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_ex (int i, int j, int k, Array4<Real> const& Ex,
                        const Real* const sigma_fac_y,
                        const Real* const sigma_fac_z,
                        const Real* const sigma_star_fac_x,
                        int xlo,int ylo, int zlo)
{
#if (AMREX_SPACEDIM == 3)
    Ex(i,j,k,0) *= sigma_fac_y[j-ylo];
    Ex(i,j,k,1) *= sigma_fac_z[k-zlo];
#else
    Ex(i,j,k,1) *= sigma_fac_z[j-zlo];
#endif
    Ex(i,j,k,2) *= sigma_star_fac_x[i-xlo];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_ey (int i, int j, int k, Array4<Real> const& Ey,
                        const Real* const sigma_fac_z,
                        const Real* const sigma_fac_x,
                        const Real* const sigma_star_fac_y,
                        int xlo,int ylo, int zlo)
{
#if (AMREX_SPACEDIM == 3)
    Ey(i,j,k,0) *= sigma_fac_z[k-zlo];
    Ey(i,j,k,2) *= sigma_star_fac_y[j-ylo];
#else
    Ey(i,j,k,0) *= sigma_fac_z[j-zlo];
#endif
    Ey(i,j,k,1) *= sigma_fac_x[i-xlo];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_ez (int i, int j, int k, Array4<Real> const& Ez,
                        const Real* const sigma_fac_x,
                        const Real* const sigma_fac_y,
                        const Real* const sigma_star_fac_z,
                        int xlo,int ylo, int zlo)
{
    Ez(i,j,k,0) *= sigma_fac_x[i-xlo];
#if (AMREX_SPACEDIM == 3)
    Ez(i,j,k,1) *= sigma_fac_y[j-ylo];
    Ez(i,j,k,2) *= sigma_star_fac_z[k-zlo];
#else
    Ez(i,j,k,2) *= sigma_star_fac_z[j-zlo];
#endif
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_bx (int i, int j, int k, Array4<Real> const& Bx,
                        const Real* const sigma_star_fac_y,
                        const Real* const sigma_star_fac_z,
                        int ylo, int zlo)
{
#if (AMREX_SPACEDIM == 3)
   Bx(i,j,k,0) *= sigma_star_fac_y[j-ylo];
   Bx(i,j,k,1) *= sigma_star_fac_z[k-zlo];
#else
   Bx(i,j,k,1) *= sigma_star_fac_z[j-zlo];
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_by (int i, int j, int k, Array4<Real> const& By,
                        const Real* const sigma_star_fac_z,
                        const Real* const sigma_star_fac_x,
                        int zlo, int xlo)
{
#if (AMREX_SPACEDIM == 3)
   By(i,j,k,0) *= sigma_star_fac_z[k-zlo];
#else
   By(i,j,k,0) *= sigma_star_fac_z[j-zlo];
#endif
   By(i,j,k,1) *= sigma_star_fac_x[i-xlo];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_bz (int i, int j, int k, Array4<Real> const& Bz,
                        const Real* const sigma_star_fac_x,
                        const Real* const sigma_star_fac_y,
                        int xlo, int ylo)
{
   Bz(i,j,k,0) *= sigma_star_fac_x[i-xlo];
#if (AMREX_SPACEDIM == 3)
   Bz(i,j,k,1) *= sigma_star_fac_y[j-ylo];
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_damp_pml_F (int i, int j, int k, Array4<Real> const& F_fab,
                      const Real* const sigma_fac_x,
                      const Real* const sigma_fac_y,
                      const Real* const sigma_fac_z,
                      int xlo, int ylo, int zlo)
{
   F_fab(i,j,k,0) *= sigma_fac_x[i-xlo];
#if (AMREX_SPACEDIM == 3)
   F_fab(i,j,k,1) *= sigma_fac_y[j-ylo];
   F_fab(i,j,k,2) *= sigma_fac_z[k-zlo];
#else
   F_fab(i,j,k,2) *= sigma_fac_z[j-zlo];
#endif
}

#endif
