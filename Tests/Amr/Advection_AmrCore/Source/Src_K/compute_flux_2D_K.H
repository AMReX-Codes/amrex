#ifndef _compute_flux_2d_H_
#define _compute_flux_2d_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_x(int i, int j, int k,
            amrex::Array4<amrex::Real> const& px,
            amrex::Array4<amrex::Real const> const& phi,
            amrex::Array4<amrex::Real const> const& vx,
            amrex::Array4<amrex::Real const> const& slope,
            amrex::Real dtdx)
{
    px(i,j,k) = ( (vx(i,j,k) < 0) ?
                phi(i  ,j,k) - slope(i  ,j,k)*(0.5 + 0.5*dtdx*vx(i,j,k)) :
                phi(i-1,j,k) + slope(i-1,j,k)*(0.5 - 0.5*dtdx*vx(i,j,k)) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_y(int i, int j, int k,
            amrex::Array4<amrex::Real> const& py,
            amrex::Array4<amrex::Real const> const& phi,
            amrex::Array4<amrex::Real const> const& vy,
            amrex::Array4<amrex::Real const> const& slope,
            amrex::Real dtdy)
{
    py(i,j,k) = ( (vy(i,j,k) < 0) ?
                phi(i,j  ,k) - slope(i,j  ,k)*(0.5 + 0.5*dtdy*vy(i,j,k)) :
                phi(i,j-1,k) + slope(i,j-1,k)*(0.5 - 0.5*dtdy*vy(i,j,k)) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void create_flux_x(int i, int j, int k,
                   amrex::Array4<amrex::Real> const& fx,
                   amrex::Array4<amrex::Real const> const& vx,
                   amrex::Array4<amrex::Real const> const& vy,
                   amrex::Array4<amrex::Real const> const& px,
                   amrex::Array4<amrex::Real const> const& py,
                   amrex::Real dtdy)
{
    fx(i,j,k) = ( (vx(i,j,k) < 0) ?
                (px(i,j,k) - 0.5*dtdy * ( 0.5*(vy(i  ,j+1,k  ) + vy(i  ,j,k)) * (py(i  ,j+1,k  )-py(i  ,j,k))))*vx(i,j,k) :
                (px(i,j,k) - 0.5*dtdy * ( 0.5*(vy(i-1,j+1,k  ) + vy(i-1,j,k)) * (py(i-1,j+1,k  )-py(i-1,j,k))))*vx(i,j,k) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void create_flux_y(int i, int j, int k,
                   amrex::Array4<amrex::Real> const& fy,
                   amrex::Array4<amrex::Real const> const& vx,
                   amrex::Array4<amrex::Real const> const& vy,
                   amrex::Array4<amrex::Real const> const& py,
                   amrex::Array4<amrex::Real const> const& px,
                   amrex::Real dtdx)
{
    fy(i,j,k) = ( (vy(i,j,k) < 0) ?
                (py(i,j,k) - 0.5*dtdx * ( 0.5*(vx(i+1,j  ,k  ) + vx(i,j  ,k)) * (px(i+1,j  ,k  )-px(i,j  ,k))))*vy(i,j,k) :
                (py(i,j,k) - 0.5*dtdx * ( 0.5*(vx(i+1,j-1,k  ) + vx(i,j-1,k)) * (px(i+1,j-1,k  )-px(i,j-1,k))))*vy(i,j,k) );
}

#endif
