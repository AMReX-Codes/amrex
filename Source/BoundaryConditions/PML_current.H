#ifndef PML_CURRENT_H_
#define PML_CURRENT_H_

#include <AMReX_FArrayBox.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_push_ex_pml_current_nodal (int j, int k, int l, Array4<Real> const& Ex,
                                      Array4<Real const> const& jx,
                                      Array4<Real const> const& sigjy,
                                      Array4<Real const> const& sigjz,
                                      Real mudt, IntVect pml_type)
{
#if (AMREX_SPACEDIM == 3)
    Real alpha_xy, alpha_xz;
    if (pml_type[0] == 0){
        if (pml_type[1] == 0){
          alpha_xy = 0.5;
          alpha_xz = 0.5;
        }
        if (pml_type[1] == 1){
          alpha_xy = 1.;
          alpha_xz = 0.;
        }
        if (pml_type[1] == 2){
          alpha_xy = 0.;
          alpha_xz = 1.;
        }
    }
    else {
        if ((sigjy[k] == 0.) && (sigjz[l] == 0.)){
            alpha_xy = 0.5;
            alpha_xz = 0.5;
        }
        else {
            alpha_xy = sigjy[k]/(sigjy[k]+sigjz[l]);
            alpha_xz = sigjz[l]/(sigjy[k]+sigjz[l]);
        }
    }
    Ex(j,k,l,0) = Ex(j,k,l,0) - mudt  * alpha_xy * jx(j,k,l);
    Ex(j,k,l,1) = Ex(j,k,l,1) - mudt  * alpha_xz * jx(j,k,l);
#else

    Ex(j,k,0) = Ex(j,k,0) - mudt  * jx(j,k);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_push_ey_pml_current_nodal (int j, int k, int l, Array4<Real> const& Ey,
                                      Array4<Real const> const& jy,
                                      Array4<Real const> const& sigjx,
                                      Array4<Real const> const& sigjz,
                                      Real mudt, int pml_type)
{
#if (AMREX_SPACEDIM == 3)
    Real alpha_yx, alpha yz;
    if (pml_type[0] == 0){
        if (pml_type[1] == 0){
          alpha_yx = 1.;
          alpha_yz = 0.;
        }
        if (pml_type[1] == 1){
          alpha_yx = 0.5;
          alpha_yz = 0.5;
        }
        if (pml_type[1] == 2){
          alpha_yx = 0.;
          alpha_yz = 1.;
        }
    }
    else {
        if ((sigjy[k] == 0.) && (sigjz[l] == 0.)){
            alpha_yx = 0.5;
            alpha_yz = 0.5;
        }
        else {
            alpha_yx = sigjx[j]/(sigjx[j]+sigjz[l]);
            alpha_yz = sigjz[l]/(sigjx[j]+sigjz[l]);
        }
    }
    Ey(j,k,l,0) = Ey(j,k,l,0) - mudt  * alpha_yx * jy(j,k,l);
    Ey(j,k,l,1) = Ey(j,k,l,1) - mudt  * alpha_yz * jy(j,k,l);
#else
    Ey(j,k,0) = Ey(j,k,0) - 0.5 * mudt  * jy(j,k);
    Ey(j,k,1) = Ey(j,k,0) - 0.5 * mudt  * jy(j,k);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void warpx_push_ez_pml_current_nodal (int j, int k, int l, Array4<Real> const& Ez,
                                      Array4<Real const> const& jz,
                                      Array4<Real const> const& sigjx,
                                      Array4<Real const> const& sigjy,
                                      Real mudt, int pml_type)
{
#if (AMREX_SPACEDIM == 3)
    Real alpha_zx, alpha_zy;
    if (pml_type[0] == 0){
        if (pml_type[1] == 0){
          alpha_zx = 1.;
          alpha_zy = 0.;
        }
        if (pml_type[1] == 1){
          alpha_zx = 0.;
          alpha_zy = 1.;
        }
        if (pml_type[1] == 2){
          alpha_zx = 0.5;
          alpha_zy = 0.5;
        }
    }
    else {
        if ((sigjy[k] == 0.) && (sigjz[l] == 0.)){
            alpha_zx = 0.5;
            alpha_zy = 0.5;
        }
        else {
            alpha_zx = sigjx[j]/(sigjx[j]+sigjy[k]);
            alpha_zy = sigjy[k]/(sigjx[j]+sigjy[k]);
        }
    }
    Ez(j,k,l,0) = Ez(j,k,l,0) - mudt  * alpha_zx * jz(j,k,l);
    Ez(j,k,l,1) = Ez(j,k,l,1) - mudt  * alpha_zy * jz(j,k,l);
#else
    Ez(j,k,0) = Ez(j,k,0) - mudt  * jz(j,k);
#endif
}

#endif
