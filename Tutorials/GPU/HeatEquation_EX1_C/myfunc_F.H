#ifndef MYFUNC_F_H__ 
#define MUFUNC_F_H__ 

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void init_phi(const int* lo, const int* hi,
                  amrex_real* data, const int* dlo, const int* dhi,
                  const amrex_real* dx, const amrex_real* prob_lo, const amrex_real* prob_hi);
    
    AMREX_FORT_LAUNCHABLE(void compute_flux)(const int* lo, const int* hi,
                                             amrex_real* phi, const int* philo, const int* phihi,
                                             amrex_real* flux, const int* flo, const int* fhi,
                                             const amrex_real* dx, const int idir);

    AMREX_FORT_LAUNCHABLE(void update_phi)(const int* lo, const int* hi,
                                           amrex_real* pold, const int* polo, const int* pohi,
                                           amrex_real* pnew, const int* pnlo, const int* pnhi,
                                           amrex_real* fluxx, const int* fxlo, const int* fxhi,
                                           amrex_real* fluxy, const int* fylo, const int* fyhi,
#if (BL_SPACEDIM == 3)
                                           amrex_real* fluxz, const int* fzlo, const int* fzhi,
#endif
                                           const amrex_real* dx, const amrex_real dt);

#ifdef __cplusplus
}
#endif

#endif
