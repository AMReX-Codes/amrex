#ifndef _FACE_VELOCITY_H_
#define _FACE_VELOCITY_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AmrCoreAdv_F.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_INLINE
void get_face_velocity_psi(Box const& bx,
                           const Real time,
                           Array4<Real> const& psi,
                           GeometryData const& geomdata)
{
    const auto lo  = lbound(bx);
    const auto hi  = ubound(bx);

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize(); 

    for     (int j = lo.y; j <= hi.y; ++j) {
        Real y = dx[1]*(0.5+j) + prob_lo[1]; 
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real x = dx[0]*(0.5+i) + prob_lo[0];
            psi(i,j,0) = pow(sin(M_PI*x), 2) * pow(sin(M_PI*y), 2)
                       * cos(M_PI*time/2.0) * 1.0/M_PI; 
        }
    }
}

AMREX_GPU_DEVICE
AMREX_INLINE
void get_face_velocity_x(int i, int j, int k,
                         Array4<Real> const& vx,
                         Array4<Real> const& psi,
                         GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                         GpuArray<Real, AMREX_SPACEDIM> dx)
{
    vx(i,j,k) = -( (psi(i,j+1,0)+psi(i-1,j+1,0)) - (psi(i,j-1,0)+psi(i-1,j-1,0)) ) * (0.25/dx[1]); 
}

AMREX_GPU_DEVICE
AMREX_INLINE
void get_face_velocity_y(int i, int j, int k,
                         Array4<Real> const& vy,
                         Array4<Real> const& psi,
                         GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                         GpuArray<Real, AMREX_SPACEDIM> dx)
{
    vy(i,j,k) =  ( (psi(i+1,j,0)+psi(i+1,j-1,0)) - (psi(i-1,j,0)+psi(i-1,j-1,0)) ) * (0.25/dx[0]); 
}

AMREX_GPU_DEVICE
AMREX_INLINE
void get_face_velocity_z(int i, int j, int k,
                         Array4<Real> const& vz,
                         Array4<Real> const& psi,
                         GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                         GpuArray<Real, AMREX_SPACEDIM> dx)
{
    vz(i,j,k) =  1.0; 
}

#endif
