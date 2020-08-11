#ifndef _compute_flux_2d_H_
#define _compute_flux_2d_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_x(int i, int j, int k,
            Array4<Real> const& phi,
            Array4<Real> const& vx,
            Array4<Real> const& px,
            Array4<Real> const& slope,
            const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    px(i,j,k) = ( (vx(i,j,k) < 0) ? 
                phi(i  ,j,k) - slope(i  ,j,k)*(0.5 + 0.5*dtdx[0]*vx(i,j,k)) : 
                phi(i-1,j,k) + slope(i-1,j,k)*(0.5 - 0.5*dtdx[0]*vx(i,j,k)) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_y(int i, int j, int k,
            Array4<Real> const& phi,
            Array4<Real> const& vy,
            Array4<Real> const& py,
            Array4<Real> const& slope,
            const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    py(i,j,k) = ( (vy(i,j,k) < 0) ? 
                phi(i,j  ,k) - slope(i,j  ,k)*(0.5 + 0.5*dtdx[0]*vy(i,j,k)) : 
                phi(i,j-1,k) + slope(i,j-1,k)*(0.5 - 0.5*dtdx[0]*vy(i,j,k)) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void create_flux_x(int i, int j, int k,
                   AMREX_D_DECL(Array4<Real> const& vx, 
                                Array4<Real> const& vy,
                                Array4<Real> const& vz),
                   Array4<Real> const& px,
                   Array4<Real> const& py,
                   Array4<Real> const& fx,
                   const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    fx(i,j,k) = ( (vx(i,j,k) < 0) ? 
                (px(i,j,k) - 0.5*dtdx[1] * ( 0.5*(vy(i  ,j+1,k  ) + vy(i  ,j,k)) * (py(i  ,j+1,k  )-py(i  ,j,k))))*vx(i,j,k) :
                (px(i,j,k) - 0.5*dtdx[1] * ( 0.5*(vy(i-1,j+1,k  ) + vy(i-1,j,k)) * (py(i-1,j+1,k  )-py(i-1,j,k))))*vx(i,j,k) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void create_flux_y(int i, int j, int k,
                   AMREX_D_DECL(Array4<Real> const& vx, 
                                Array4<Real> const& vy,
                                Array4<Real> const& vz),
                   Array4<Real> const& py,
                   Array4<Real> const& px,
                   Array4<Real> const& fy,
                   const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    fy(i,j,k) = ( (vy(i,j,k) < 0) ? 
                (py(i,j,k) - 0.5*dtdx[0] * ( 0.5*(vx(i+1,j  ,k  ) + vx(i,j  ,k)) * (px(i+1,j  ,k  )-px(i,j  ,k))))*vy(i,j,k) :
                (py(i,j,k) - 0.5*dtdx[0] * ( 0.5*(vx(i+1,j-1,k  ) + vx(i,j-1,k)) * (px(i+1,j-1,k  )-px(i,j-1,k))))*vy(i,j,k) );
}

#endif
