#ifndef _Adv_3d_H_
#define _Adv_3d_H_

#include <AmrCoreAdv_F.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_INLINE
void conservative(int i, int j, int k, 
                  Array4<Real> const& uin,
                  Array4<Real> const& uout,
                  AMREX_D_DECL(Array4<Real> const& flxx,
                               Array4<Real> const& flxy,
                               Array4<Real> const& flxz),
                  const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    uout(i,j,k) = uin(i,j,k) + 
                  ( (flxx(i,j,k) - flxx(i+1,j,k)) * dtdx[0] 
                  + (flxy(i,j,k) - flxy(i,j+1,k)) * dtdx[1] 
                  + (flxz(i,j,k) - flxz(i,j,k+1)) * dtdx[2] );
}

AMREX_GPU_DEVICE
AMREX_INLINE
void flux_scale_x(int i, int j, int k,
                  Array4<Real> const& flxx,
                  Real dt,
                  const GpuArray<Real, AMREX_SPACEDIM>& dx)
{
   flxx(i,j,k) = flxx(i,j,k) * (dt * dx[1]*dx[2]);
}

AMREX_GPU_DEVICE
AMREX_INLINE
void flux_scale_y(int i, int j, int k,
                  Array4<Real> const& flxy,
                  Real dt,
                  const GpuArray<Real, AMREX_SPACEDIM>& dx)
{
   flxy(i,j,k) = flxy(i,j,k) * (dt * dx[0]*dx[2]);
}

AMREX_GPU_DEVICE
AMREX_INLINE
void flux_scale_z(int i, int j ,int k,
                  Array4<Real> const& flxz,
                  Real dt,
                  const GpuArray<Real, AMREX_SPACEDIM>& dx)
{
    flxz(i,j,k) = flxz(i,j,k) * (dt * dx[0]*dx[1]);
}

#endif
