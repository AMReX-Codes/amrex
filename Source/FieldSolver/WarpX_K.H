#ifndef WARPX_K_H_
#define WARPX_K_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_push_bx_nodal (int j, int k, int l,
                          amrex::Array4<amrex::Real> const& Bx,
                          amrex::Array4<amrex::Real const> const& Ey,
                          amrex::Array4<amrex::Real const> const& Ez,
                          amrex::Real dtsdy, amrex::Real dtsdz)
{
#if (AMREX_SPACEDIM == 3)
    Bx(j,k,l) = Bx(j,k,l) - 0.5*dtsdy * (Ez(j,k+1,l  ) - Ez(j,k-1,l  ))
                          + 0.5*dtsdz * (Ey(j,k  ,l+1) - Ey(j,k  ,l-1));
#else
    Bx(j,k,0) = Bx(j,k,0) + 0.5*dtsdz * (Ey(j,k+1,0) - Ey(j,k-1,0));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_push_by_nodal (int j, int k, int l,
                          amrex::Array4<amrex::Real> const& By,
                          amrex::Array4<amrex::Real const> const& Ex,
                          amrex::Array4<amrex::Real const> const& Ez,
                          amrex::Real dtsdx,
                          amrex::Real dtsdz)
{
#if (AMREX_SPACEDIM == 3)
    By(j,k,l) = By(j,k,l) + 0.5*dtsdx * (Ez(j+1,k,l  ) - Ez(j-1,k,l  ))
                          - 0.5*dtsdz * (Ex(j  ,k,l+1) - Ex(j  ,k,l-1));
#else
    By(j,k,0) = By(j,k,0) + 0.5*dtsdx * (Ez(j+1,k  ,0) - Ez(j-1,k  ,0))
                          - 0.5*dtsdz * (Ex(j  ,k+1,0) - Ex(j  ,k-1,0));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_push_bz_nodal (int j, int k, int l,
                          amrex::Array4<amrex::Real> const& Bz,
                          amrex::Array4<amrex::Real const> const& Ex,
                          amrex::Array4<amrex::Real const> const& Ey,
                          amrex::Real dtsdx, amrex::Real dtsdy)
{
#if (AMREX_SPACEDIM == 3)
    Bz(j,k,l) = Bz(j,k,l) - 0.5*dtsdx * (Ey(j+1,k  ,l) - Ey(j-1,k  ,l))
                          + 0.5*dtsdy * (Ex(j  ,k+1,l) - Ex(j  ,k-1,l));
#else
    Bz(j,k,0) = Bz(j,k,0) - 0.5*dtsdx * (Ey(j+1,k  ,0) - Ey(j-1,k  ,0));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_push_ex_nodal (int j, int k, int l,
                          amrex::Array4<amrex::Real> const& Ex,
                          amrex::Array4<amrex::Real const> const& By,
                          amrex::Array4<amrex::Real const> const& Bz,
                          amrex::Array4<amrex::Real const> const& jx,
                          amrex::Real mudt, amrex::Real dtsdy, amrex::Real dtsdz)
{
#if (AMREX_SPACEDIM == 3)
    Ex(j,k,l) = Ex(j,k,l) + 0.5*dtsdy * (Bz(j,k+1,l  ) - Bz(j,k-1,l  ))
                          - 0.5*dtsdz * (By(j,k  ,l+1) - By(j,k  ,l-1))
        - mudt  * jx(j,k,l);
#else
    Ex(j,k,0) = Ex(j,k,0) - 0.5*dtsdz * (By(j,k+1,0) - By(j,k-1,0))
        - mudt  * jx(j,k,0);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_push_ey_nodal (int j, int k, int l,
                          amrex::Array4<amrex::Real> const& Ey,
                          amrex::Array4<amrex::Real const> const& Bx,
                          amrex::Array4<amrex::Real const> const& Bz,
                          amrex::Array4<amrex::Real const> const& jy,
                          amrex::Real mudt, amrex::Real dtsdx, amrex::Real dtsdz)
{
#if (AMREX_SPACEDIM == 3)
    Ey(j,k,l) = Ey(j,k,l) - 0.5*dtsdx * (Bz(j+1,k,l  ) - Bz(j-1,k,l  ))
                          + 0.5*dtsdz * (Bx(j  ,k,l+1) - Bx(j  ,k,l-1))
        - mudt  * jy(j,k,l);
#else
    Ey(j,k,0) = Ey(j,k,0) - 0.5*dtsdx * (Bz(j+1,k  ,0) - Bz(j-1,k  ,0))
                          + 0.5*dtsdz * (Bx(j  ,k+1,0) - Bx(j  ,k-1,0))
        - mudt  * jy(j,k,0);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_push_ez_nodal (int j, int k, int l,
                          amrex::Array4<amrex::Real> const& Ez,
                          amrex::Array4<amrex::Real const> const& Bx,
                          amrex::Array4<amrex::Real const> const& By,
                          amrex::Array4<amrex::Real const> const& jz,
                          amrex::Real mudt, amrex::Real dtsdx, amrex::Real dtsdy)
{
#if (AMREX_SPACEDIM == 3)
    Ez(j,k,l) = Ez(j,k,l) + 0.5*dtsdx * (By(j+1,k  ,l) - By(j-1,k  ,l))
                          - 0.5*dtsdy * (Bx(j  ,k+1,l) - Bx(j  ,k-1,l))
        - mudt  * jz(j,k,l);
#else
    Ez(j,k,0) = Ez(j,k,0) + 0.5*dtsdx * (By(j+1,k,0) - By(j-1,k,0))
        - mudt  * jz(j,k,0);
#endif
}

#endif
