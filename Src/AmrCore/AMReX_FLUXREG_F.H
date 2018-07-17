#ifndef AMREX_FLUXREG_F_H_
#define AMREX_FLUXREG_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void amrex_frfineadd (amrex_real* reg,       AMREX_ARLIM_P(rlo), AMREX_ARLIM_P(rhi),
                         const amrex_real* flx, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
                         const int* numcomp, const int* dir,
                         const int* ratio, const amrex_real* mult);

    void amrex_frfaadd (amrex_real* reg,        AMREX_ARLIM_P(rlo), AMREX_ARLIM_P(rhi),
                       const amrex_real* flx,  AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
                       const amrex_real* area, AMREX_ARLIM_P(alo), AMREX_ARLIM_P(ahi),
                       const int* numcomp, const int* dir,
                       const int* ratio, const amrex_real* mult);

    void amrex_frreflux (const int* lo, const int* hi,
			      amrex_real* s, const int* slo, const int*shi,
			const amrex_real* f, const int* flo, const int*fhi,
			const amrex_real* v, const int* vlo, const int*vhi,
			const int* n, const amrex_real* mult, 
			const int* dir, const int* islo);

    void amrex_froverwrite_cfb (const int* lo, const int* hi,
                                amrex_real* dst, const int* dlo, const int* dhi,
                                const amrex_real* src, const int* slo, const int* shi,
                                const int* mask, const int* mlo, const int* mhi,
                                const int* ncomp, const int* idim, const amrex_real* scale);

#ifdef __cplusplus
}
#endif

#endif
