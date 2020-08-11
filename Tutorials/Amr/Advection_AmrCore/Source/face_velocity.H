#ifndef FACE_VELOCITY_H_
#define FACE_VELOCITY_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_psi(amrex::Box const& bx,
                           const amrex::Real time,
                           amrex::Array4<amrex::Real> const& psi,
                           amrex::GeometryData const& geomdata)
{
    using namespace amrex;
    constexpr Real PI = 3.1415926535897932384626;

    const auto lo  = lbound(bx);
    const auto hi  = ubound(bx);

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize(); 

    for     (int j = lo.y; j <= hi.y; ++j) {
        Real y = dx[1]*(0.5+j) + prob_lo[1]; 
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real x = dx[0]*(0.5+i) + prob_lo[0];
            psi(i,j,0) = std::pow(std::sin(PI*x), 2) * std::pow(std::sin(PI*y), 2)
                       * std::cos(PI*time/2.0) * 1.0/PI;
        }
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_x(int i, int j, int k,
                         amrex::Array4<amrex::Real> const& vx,
                         amrex::Array4<amrex::Real> const& psi,
                         amrex::Real dy)
{
    vx(i,j,k) = -( (psi(i,j+1,0)+psi(i-1,j+1,0)) - (psi(i,j-1,0)+psi(i-1,j-1,0)) ) * (0.25/dy);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_y(int i, int j, int k,
                         amrex::Array4<amrex::Real> const& vy,
                         amrex::Array4<amrex::Real> const& psi,
                         amrex::Real dx)
{
    vy(i,j,k) =  ( (psi(i+1,j,0)+psi(i+1,j-1,0)) - (psi(i-1,j,0)+psi(i-1,j-1,0)) ) * (0.25/dx);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_z(int i, int j, int k,
                         amrex::Array4<amrex::Real> const& vz)
{
    vz(i,j,k) =  0.0; 
}

#endif
