#ifndef PROB_H_
#define PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
initdata(int i, int j, int k, amrex::Array4<amrex::Real> const& phi,
         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    using namespace amrex;

    Real y = prob_lo[1] + (amrex::Real(0.5)+amrex::Real(j)) * dx[1];
    Real x = prob_lo[0] + (amrex::Real(0.5)+amrex::Real(i)) * dx[0];
    Real r2 = (amrex::Real(std::pow(x-amrex::Real(0.5), 2)) + amrex::Real(std::pow((y-amrex::Real(0.75)),2))) / amrex::Real(0.01);
    phi(i,j,k) = amrex::Real(1.0) + std::exp(-r2);
}

#endif
