#ifndef MYFUNC_F_H__ 
#define MUFUNC_F_H__ 

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

#if defined (AMREX_USE_ACC) && !defined (AMREX_OMP_OFFLOAD)

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

#elif !defined (AMREX_USE_ACC) && defined (AMREX_OMP_OFFLOAD)

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

#else
    AMREX_GPU_DEVICE
    void init_phi(const int* lo, const int* hi,
                  amrex_real* data, const int* dlo, const int* dhi,
                  const amrex_real* dx, const amrex_real* prob_lo);
    
    AMREX_GPU_DEVICE
    void compute_flux(const int* lo, const int* hi,
                      const int* domlo, const int* domhi,
                      const amrex_real* phi, const int* philo, const int* phihi,
                      amrex_real* fluxx, const int* fxlo, const int* fxhi,
                      amrex_real* fluxy, const int* fylo, const int* fyhi,
                      amrex_real* fluxz, const int* fzlo, const int* fzhi,
                      const amrex_real* dx);

    AMREX_GPU_DEVICE
    void update_phi(const int* lo, const int* hi,
                    const amrex_real* pold, const int* polo, const int* pohi,
                    amrex_real* pnew, const int* pnlo, const int* pnhi,
                    const amrex_real* fluxx, const int* fxlo, const int* fxhi,
                    const amrex_real* fluxy, const int* fylo, const int* fyhi,
                    const amrex_real* fluxz, const int* fzlo, const int* fzhi,
                    const amrex_real* dx, const amrex_real dt);

#endif 

#ifdef __cplusplus
}
#endif

#endif
