#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void init_phi (int i, int j, int k,
               amrex::Array4<amrex::Real> const& phi,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo)
{
    using namespace amrex::literals;

    amrex::Real x = prob_lo[0] + (i+0.5_rt) * dx[0];
    amrex::Real y = prob_lo[1] + (j+0.5_rt) * dx[1];
#if (AMREX_SPACEDIM > 2)
    amrex::Real z = prob_lo[2] + (k+0.5_rt) * dx[2];
#else
    amrex::Real z = 0._rt;
#endif
    amrex::Real r2 = ((x-0.25_rt)*(x-0.25_rt)+(y-0.25_rt)*(y-0.25_rt)+(z-0.25_rt)*(z-0.25_rt))/0.01_rt;
#ifdef AMREX_USE_HIP
    phi(i,j,k) = 1._rt + ::exp(-r2); // xxxxx HIP todo
#else
    phi(i,j,k) = 1._rt + std::exp(-r2);
#endif
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_x (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxx,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dxinv)
{
    fluxx(i,j,k) = (phi(i,j,k)-phi(i-1,j,k)) * dxinv;
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_y (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxy,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dyinv)
{
    fluxy(i,j,k) = (phi(i,j,k)-phi(i,j-1,k)) * dyinv;
}


#if (AMREX_SPACEDIM > 2)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_z (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxz,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dzinv)
{
    fluxz(i,j,k) = (phi(i,j,k)-phi(i,j,k-1)) * dzinv;
}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void update_phi (int i, int j, int k,
                 amrex::Array4<amrex::Real const> const& phiold,
                 amrex::Array4<amrex::Real      > const& phinew,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fluxx,
                              amrex::Array4<amrex::Real const> const& fluxy,
                              amrex::Array4<amrex::Real const> const& fluxz),
                 amrex::Real dt,
                 AMREX_D_DECL(amrex::Real dxinv,
                              amrex::Real dyinv,
                              amrex::Real dzinv))
{
    phinew(i,j,k) = phiold(i,j,k)
        + dt * dxinv * (fluxx(i+1,j  ,k  ) - fluxx(i,j,k))
        + dt * dyinv * (fluxy(i  ,j+1,k  ) - fluxy(i,j,k))
#if (AMREX_SPACEDIM > 2)
        + dt * dzinv * (fluxz(i  ,j  ,k+1) - fluxz(i,j,k));
#else
        ;
#endif
}

#endif
