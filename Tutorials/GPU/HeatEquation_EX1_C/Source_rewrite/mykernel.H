#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_BLFort.H>

AMREX_GPU_DEVICE AMREX_INLINE
void init_phi (int i, int j, int k,
               amrex::Array4<amrex::Real>& phi,
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
                     amrex::Array4<amrex::Real>& fluxx,
                     amrex::Array4<amrex::Real> const& phi, amrex::Real dxinv)
{
    fluxx(i,j,k) = (phi(i,j,k)-phi(i-1,j,k)) * dxinv;
}


AMREX_GPU_DEVICE AMREX_INLINE
void compute_flux_y (int i, int j, int k,
                     amrex::Array4<amrex::Real>& fluxy,
                     amrex::Array4<amrex::Real> const& phi, amrex::Real dyinv)
{
    fluxy(i,j,k) = (phi(i,j,k)-phi(i,j-1,k)) * dyinv;
}


AMREX_GPU_DEVICE AMREX_INLINE
void compute_flux_z (int i, int j, int k,
                     amrex::Array4<amrex::Real>& fluxz,
                     amrex::Array4<amrex::Real> const& phi, amrex::Real dzinv)
{
    fluxz(i,j,k) = (phi(i,j,k)-phi(i,j,k-1)) * dzinv;
}


AMREX_GPU_DEVICE AMREX_INLINE
void update_phi (int i, int j, int k,
                 amrex::Array4<amrex::Real> const& phiold,
                 amrex::Array4<amrex::Real>      & phinew,
                 amrex::Array4<amrex::Real> const& fluxx,
                 amrex::Array4<amrex::Real> const& fluxy,
                 amrex::Array4<amrex::Real> const& fluxz,
                 amrex::Real dt,
                 amrex::Real dxinv, amrex::Real dyinv, amrex::Real dzinv)
{
    phinew(i,j,k) = phiold(i,j,k)
        + dt * dxinv * (fluxx(i+1,j  ,k  ) - fluxx(i,j,k)),
        + dt * dyinv * (fluxy(i  ,j+1,k  ) - fluxy(i,j,k)),
        + dt * dzinv * (fluxz(i  ,j  ,k+1) - fluxz(i,j,k));
}

#ifdef __cplusplus
extern "C"{
#endif

#ifdef AMREX_USE_ACC
void init_phi_acc(const int *lo, const int *hi,
                  amrex_real *dat, const int *dlo, const int *dhi,
                  const amrex_real *dx,
                  const amrex_real *plo);

void compute_flux_x_acc(const int *lo, const int *hi,
                    amrex_real *fluxx, const int *f_lo, const int *f_hi,
                    amrex_real *phi, const int *p_lo, const int *p_hi,
                    amrex_real dxinv);

void compute_flux_y_acc(const int *lo, const int *hi,
                    amrex_real *fluxy, const int *f_lo, const int *f_hi,
                    amrex_real *phi, const int *p_lo, const int *p_hi,
                    amrex_real dyinv);

void compute_flux_z_acc(const int *lo, const int *hi,
                    amrex_real *fluxz, const int *f_lo, const int *f_hi,
                    amrex_real *phi, const int *p_lo, const int *p_hi,
                    amrex_real dzinv);

void update_phi_acc(const int *lo, const int *hi,
                const amrex_real *fluxx, const int *fxlo, const int *fxhi,
                const amrex_real *fluxy, const int *fylo, const int *fyhi,
                const amrex_real *fluxz, const int *fzlo, const int *fzhi,
                const amrex_real *phi_old, const int *polo, const int *pohi,
                amrex_real *phi_new, const int *pnlo, const int *pnhi,
                const amrex_real dt, const amrex_real dxinv, const amrex_real dyinv, const amrex_real dzinv);
#endif

#ifdef AMREX_OMP_OFFLOAD
void init_phi_omp(const int *lo, const int *hi,
                  amrex_real *dat, const int *dlo, const int *dhi,
                  const  amrex_real *dx, const amrex_real *prob_lo);

void compute_flux_x_omp(const int *lo, const int *hi,
                    amrex_real *fluxx, const int *f_lo, const int *f_hi,
                    amrex_real *phi, const int *p_lo, const int *p_hi,
                    amrex_real dxinv);

void compute_flux_y_omp(const int *lo, const int *hi,
                    amrex_real *fluxy, const int *f_lo, const int *f_hi,
                    amrex_real *phi, const int *p_lo, const int *p_hi,
                    amrex_real dyinv);

void compute_flux_z_omp(const int *lo, const int *hi,
                    amrex_real *fluxz, const int *f_lo, const int *f_hi,
                    amrex_real *phi, const int *p_lo, const int *p_hi,
                    amrex_real dzinv);

void update_phi_omp(const int *lo, const int *hi,
                const amrex_real *fluxx, const int *fxlo, const int *fxhi,
                const amrex_real *fluxy, const int *fylo, const int *fyhi,
                const amrex_real *fluxz, const int *fzlo, const int *fzhi,
                const amrex_real *phi_old, const int *polo, const int *pohi,
                amrex_real *phi_new, const int *pnlo, const int *pnhi,
                const amrex_real dt, const amrex_real dxinv, const amrex_real dyinv, const amrex_real dzinv);
#endif

#ifdef __cplusplus
}
#endif

#endif
