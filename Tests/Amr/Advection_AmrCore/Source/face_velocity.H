#ifndef FACE_VELOCITY_H_
#define FACE_VELOCITY_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <numbers>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_psi(int i, int j,
                           const amrex::Real time,
                           amrex::Array4<amrex::Real> const& psi,
                           amrex::GeometryData const& geomdata)
{
    using namespace amrex;
    constexpr Real PI = std::numbers::pi_v<Real>;

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize();

    Real y = dx[1]*(amrex::Real(0.5)+amrex::Real(j)) + prob_lo[1];
    Real x = dx[0]*(amrex::Real(0.5)+amrex::Real(i)) + prob_lo[0];
    psi(i,j,0) = amrex::Real(std::pow(std::sin(PI*x), 2)) * amrex::Real(std::pow(std::sin(PI*y), 2))
                * std::cos(PI*time/amrex::Real(2.0)) * amrex::Real(1.0)/PI;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_x(int i, int j, int k,
                         amrex::Array4<amrex::Real> const& vx,
                         amrex::Array4<amrex::Real> const& psi,
                         amrex::Real dy)
{
    vx(i,j,k) = -( (psi(i,j+1,0)+psi(i-1,j+1,0)) - (psi(i,j-1,0)+psi(i-1,j-1,0)) ) * (amrex::Real(0.25)/dy);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_y(int i, int j, int k,
                         amrex::Array4<amrex::Real> const& vy,
                         amrex::Array4<amrex::Real> const& psi,
                         amrex::Real dx)
{
    vy(i,j,k) =  ( (psi(i+1,j,0)+psi(i+1,j-1,0)) - (psi(i-1,j,0)+psi(i-1,j-1,0)) ) * (amrex::Real(0.25)/dx);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_z(int i, int j, int k,
                         amrex::Array4<amrex::Real> const& vz)
{
    vz(i,j,k) =  0.0;
}

#endif
