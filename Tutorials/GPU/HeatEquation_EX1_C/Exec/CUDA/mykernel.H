#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_INLINE
void init_phi (int i, int j, int k,
               amrex::Array4<amrex::Real> const& phi,
               GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo)
{
    amrex::Real x = prob_lo[0] + (i+0.5) * dx[0];
    amrex::Real y = prob_lo[1] + (j+0.5) * dx[1];
    amrex::Real z = prob_lo[2] + (k+0.5) * dx[2];
    amrex::Real r2 = ((x-0.25)*(x-0.25)+(y-0.25)*(y-0.25)+(z-0.25)*(z-0.25))/0.01;
    phi(i,j,k) = 1. + std::exp(-r2);
}


AMREX_GPU_DEVICE AMREX_INLINE
void compute_flux_x (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxx,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dxinv)
{
    fluxx(i,j,k) = (phi(i,j,k)-phi(i-1,j,k)) * dxinv;
}


AMREX_GPU_DEVICE AMREX_INLINE
void compute_flux_y (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxy,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dyinv)
{
    fluxy(i,j,k) = (phi(i,j,k)-phi(i,j-1,k)) * dyinv;
}


AMREX_GPU_DEVICE AMREX_INLINE
void compute_flux_z (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxz,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dzinv)
{
    fluxz(i,j,k) = (phi(i,j,k)-phi(i,j,k-1)) * dzinv;
}


AMREX_GPU_DEVICE AMREX_INLINE
void update_phi (int i, int j, int k,
                 amrex::Array4<amrex::Real const> const& phiold,
                 amrex::Array4<amrex::Real      > const& phinew,
                 amrex::Array4<amrex::Real const> const& fluxx,
                 amrex::Array4<amrex::Real const> const& fluxy,
                 amrex::Array4<amrex::Real const> const& fluxz,
                 amrex::Real dt,
                 amrex::Real dxinv, amrex::Real dyinv, amrex::Real dzinv)
{
    phinew(i,j,k) = phiold(i,j,k) + dt * dxinv * (fluxx(i+1,j  ,k  ) - fluxx(i,j,k)) + dt * dyinv * (fluxy(i  ,j+1,k  ) - fluxy(i,j,k)) + dt * dzinv * (fluxz(i  ,j  ,k+1) - fluxz(i,j,k));
}

#endif
