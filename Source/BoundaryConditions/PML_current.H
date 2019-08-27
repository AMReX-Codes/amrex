#ifndef PML_CURRENT_H_
#define PML_CURRENT_H_

#include <AMReX_FArrayBox.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_push_ex_pml_current_nodal (int j, int k, int l, Array4<Real> const& Ex,
                                      Array4<Real const> const& jx,
                                      Real const* const& sigjy,
                                      Real const* const& sigjz,
                                      Real mu_c2_dt)
{
#if (AMREX_SPACEDIM == 3)
    Real alpha_xy, alpha_xz;
    if (sigjy[k]+sigjz[l] == 0){
        alpha_xy = 0.5;
        alpha_xz = 0.5;
    }
    else {
        alpha_xy = sigjy[k]/(sigjy[k]+sigjz[l]);
        alpha_xz = sigjz[l]/(sigjy[k]+sigjz[l]);
    }
    Ex(j,k,l,0) = Ex(j,k,l,0) - mu_c2_dt  * alpha_xy * jx(j,k,l);
    Ex(j,k,l,1) = Ex(j,k,l,1) - mu_c2_dt  * alpha_xz * jx(j,k,l);
#else
    Ex(j,k,l,1) = Ex(j,k,l,1) - mu_c2_dt  * jx(j,k,l);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_push_ey_pml_current_nodal (int j, int k, int l, Array4<Real> const& Ey,
                                      Array4<Real const> const& jy,
                                      Real const* const& sigjx,
                                      Real const* const& sigjz,
                                      Real mu_c2_dt)
{
#if (AMREX_SPACEDIM == 3)
Real alpha_yx, alpha_yz;
    if (sigjx[j]+sigjz[l] == 0){
        alpha_yx = 0.5;
        alpha_yz = 0.5;
    }
    else {
        alpha_yx = sigjx[j]/(sigjx[j]+sigjz[l]);
        alpha_yz = sigjz[l]/(sigjx[j]+sigjz[l]);
    }
    Ey(j,k,l,0) = Ey(j,k,l,0) - mu_c2_dt  * alpha_yx * jy(j,k,l);
    Ey(j,k,l,1) = Ey(j,k,l,1) - mu_c2_dt  * alpha_yz * jy(j,k,l);
#else
    Ey(j,k,l,0) = Ey(j,k,l,0) - 0.5 * mu_c2_dt  * jy(j,k,l);
    Ey(j,k,l,1) = Ey(j,k,l,1) - 0.5 * mu_c2_dt  * jy(j,k,l);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_push_ez_pml_current_nodal (int j, int k, int l, Array4<Real> const& Ez,
                                      Array4<Real const> const& jz,
                                      Real const* const& sigjx,
                                      Real const* const& sigjy,
                                      Real mu_c2_dt)
{
#if (AMREX_SPACEDIM == 3)
    Real alpha_zx, alpha_zy;
    if (sigjx[j]+sigjy[k]==0){
        alpha_zx = 0.5;
        alpha_zy = 0.5;
    }
    else {
        alpha_zx = sigjx[j]/(sigjx[j]+sigjy[k]);
        alpha_zy = sigjy[k]/(sigjx[j]+sigjy[k]);
    }
    Ez(j,k,l,0) = Ez(j,k,l,0) - mu_c2_dt  * alpha_zx * jz(j,k,l);
    Ez(j,k,l,1) = Ez(j,k,l,1) - mu_c2_dt  * alpha_zy * jz(j,k,l);
#else
    Ez(j,k,l,0) = Ez(j,k,l,0) - mu_c2_dt  * jz(j,k,l);
#endif
}

#endif
