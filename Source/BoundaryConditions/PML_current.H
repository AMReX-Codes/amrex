/* Copyright 2019 Aurore Blelly, Axel Huebl, Maxence Thevenet
 * Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PML_CURRENT_H_
#define PML_CURRENT_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void push_ex_pml_current (int j, int k, int l,
                          amrex::Array4<amrex::Real> const& Ex,
                          amrex::Array4<amrex::Real const> const& jx,
                          amrex::Real const * const sigjy,
                          amrex::Real const * const sigjz,
                          int ylo, int zlo,
                          amrex::Real mu_c2_dt)
{
#if (AMREX_SPACEDIM == 3)
    amrex::Real alpha_xy, alpha_xz;
    if (sigjy[k-ylo]+sigjz[l-zlo] == 0){
        alpha_xy = 0.5;
        alpha_xz = 0.5;
    }
    else {
        alpha_xy = sigjy[k-ylo]/(sigjy[k-ylo]+sigjz[l-zlo]);
        alpha_xz = sigjz[l-zlo]/(sigjy[k-ylo]+sigjz[l-zlo]);
    }
    Ex(j,k,l,0) = Ex(j,k,l,0) - mu_c2_dt  * alpha_xy * jx(j,k,l);
    Ex(j,k,l,1) = Ex(j,k,l,1) - mu_c2_dt  * alpha_xz * jx(j,k,l);
#else
    Ex(j,k,l,1) = Ex(j,k,l,1) - mu_c2_dt  * jx(j,k,l);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void push_ey_pml_current (int j, int k, int l,
                          amrex::Array4<amrex::Real> const& Ey,
                          amrex::Array4<amrex::Real const> const& jy,
                          amrex::Real const * const sigjx,
                          amrex::Real const * const sigjz,
                          int xlo, int zlo,
                          amrex::Real mu_c2_dt)
{
#if (AMREX_SPACEDIM == 3)
    amrex::Real alpha_yx, alpha_yz;
    if (sigjx[j-xlo]+sigjz[l-zlo] == 0){
        alpha_yx = 0.5;
        alpha_yz = 0.5;
    }
    else {
        alpha_yx = sigjx[j-xlo]/(sigjx[j-xlo]+sigjz[l-zlo]);
        alpha_yz = sigjz[l-zlo]/(sigjx[j-xlo]+sigjz[l-zlo]);
    }
    Ey(j,k,l,0) = Ey(j,k,l,0) - mu_c2_dt  * alpha_yx * jy(j,k,l);
    Ey(j,k,l,1) = Ey(j,k,l,1) - mu_c2_dt  * alpha_yz * jy(j,k,l);
#else
    Ey(j,k,l,0) = Ey(j,k,l,0) - 0.5 * mu_c2_dt  * jy(j,k,l);
    Ey(j,k,l,1) = Ey(j,k,l,1) - 0.5 * mu_c2_dt  * jy(j,k,l);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void push_ez_pml_current (int j, int k, int l,
                          amrex::Array4<amrex::Real> const& Ez,
                          amrex::Array4<amrex::Real const> const& jz,
                          amrex::Real const * const sigjx,
                          amrex::Real const * const sigjy,
                          int xlo, int ylo,
                          amrex::Real mu_c2_dt)
{
#if (AMREX_SPACEDIM == 3)
    amrex::Real alpha_zx, alpha_zy;
    if (sigjx[j-xlo]+sigjy[k-ylo]==0){
        alpha_zx = 0.5;
        alpha_zy = 0.5;
    }
    else {
        alpha_zx = sigjx[j-xlo]/(sigjx[j-xlo]+sigjy[k-ylo]);
        alpha_zy = sigjy[k-ylo]/(sigjx[j-xlo]+sigjy[k-ylo]);
    }
    Ez(j,k,l,0) = Ez(j,k,l,0) - mu_c2_dt  * alpha_zx * jz(j,k,l);
    Ez(j,k,l,1) = Ez(j,k,l,1) - mu_c2_dt  * alpha_zy * jz(j,k,l);
#else
    Ez(j,k,l,0) = Ez(j,k,l,0) - mu_c2_dt  * jz(j,k,l);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void damp_jx_pml (int j, int k, int l,
                  amrex::Array4<amrex::Real> const& jx,
                  amrex::Real const* const sigsjx,
                  amrex::Real const* const sigjy,
                  amrex::Real const* const sigjz,
                  int xlo, int ylo, int zlo)
{
#if (AMREX_SPACEDIM == 3)
    jx(j,k,l) = jx(j,k,l) * sigsjx[j-xlo] * sigjy[k-ylo] * sigjz[l-zlo];
#else
    jx(j,k,l) = jx(j,k,l) * sigsjx[j-xlo] * sigjz[k-zlo];
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void damp_jy_pml (int j, int k, int l,
                  amrex::Array4<amrex::Real> const& jy,
                  amrex::Real const * const sigjx,
                  amrex::Real const * const sigsjy,
                  amrex::Real const * const sigjz,
                  int xlo, int ylo, int zlo)
{
#if (AMREX_SPACEDIM == 3)
    jy(j,k,l) = jy(j,k,l) * sigjx[j-xlo] * sigsjy[k-ylo] * sigjz[l-zlo];
#else
    jy(j,k,l) = jy(j,k,l) * sigjx[j-xlo] * sigjz[k-zlo];
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void damp_jz_pml (int j, int k, int l,
                  amrex::Array4<amrex::Real> const& jz,
                  amrex::Real const * const sigjx,
                  amrex::Real const * const sigjy,
                  amrex::Real const * const sigsjz,
                  int xlo, int ylo, int zlo)
{
#if (AMREX_SPACEDIM == 3)
    jz(j,k,l) = jz(j,k,l) * sigjx[j-xlo] * sigjy[k-ylo] * sigsjz[l-zlo];
#else
    jz(j,k,l) = jz(j,k,l) * sigjx[j-xlo] * sigsjz[k-zlo];
#endif
}

#endif
