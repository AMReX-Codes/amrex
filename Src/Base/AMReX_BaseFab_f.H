#ifndef AMREX_BASEFAB_F_H_
#define AMREX_BASEFAB_F_H_

#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    FORT_LAUNCHABLE(void fort_fab_copy)(const int* lo, const int* hi,
                                        amrex_real* dst, const int* dlo, const int* dhi,
                                        const amrex_real* src, const int* slo, const int* shi, const int* sblo,
                                        const int ncomp);

    FORT_LAUNCHABLE(void fort_fab_copytomem)(const int* lo, const int* hi,
                                             amrex_real* dst,
                                             const amrex_real* src, const int* slo, const int* shi,
                                             const int ncomp);

    FORT_LAUNCHABLE(void fort_fab_copyfrommem)(const int* lo, const int* hi,
                                               const amrex_real* dst, const int* dlo, const int* dhi, const int ncomp,
                                               const amrex_real* src);

    FORT_LAUNCHABLE(void fort_fab_setval)(const int* lo, const int* hi,
                                          const amrex_real* dst, const int* dlo, const int* dhi, const int ncomp,
                                          const amrex_real val);

    FORT_LAUNCHABLE(void fort_fab_norm)(const int* lo, const int* hi,
                                        const amrex_real* src, const int* slo, const int* shi, const int ncomp,
                                        const int p, amrex_real* nrm);

    FORT_LAUNCHABLE(void fort_fab_sum)(const int* lo, const int* hi,
                                       const amrex_real* src, const int* slo, const int* shi,
                                       const int ncomp, amrex_real* sm);

    FORT_LAUNCHABLE(void fort_fab_plus)(const int* lo, const int* hi,
                                        amrex_real* dst, const int* dlo, const int* dhi,
                                        const amrex_real* src, const int* slo, const int* shi, const int* sblo,
                                        const int ncomp);

    FORT_LAUNCHABLE(void fort_fab_minus)(const int* lo, const int* hi,
                                         amrex_real* dst, const int* dlo, const int* dhi,
                                         const amrex_real* src, const int* slo, const int* shi, const int* sblo,
                                         const int ncomp);

    FORT_LAUNCHABLE(void fort_fab_mult)(const int* lo, const int* hi,
                                        amrex_real* dst, const int* dlo, const int* dhi,
                                        const amrex_real* src, const int* slo, const int* shi, const int* sblo,
                                        const int ncomp);

    FORT_LAUNCHABLE(void fort_fab_divide)(const int* lo, const int* hi,
                                          amrex_real* dst, const int* dlo, const int* dhi,
                                          const amrex_real* src, const int* slo, const int* shi, const int* sblo,
                                          const int ncomp);

    FORT_LAUNCHABLE(void fort_fab_protdivide)(const int* lo, const int* hi,
                                              amrex_real* dst, const int* dlo, const int* dhi,
                                              const amrex_real* src, const int* slo, const int* shi, const int* sblo,
                                              const int ncomp);

    FORT_LAUNCHABLE(void fort_fab_invert)(const int* lo, const int* hi,
                                          amrex_real* dst, const int* dlo, const int* dhi, const int* ncomp,
                                          const amrex_real a);

    FORT_LAUNCHABLE(void fort_fab_saxpy)(const int* lo, const int* hi,
                                         amrex_real* dst, const int* dlo, const int* dhi,
                                         const amrex_real a,
                                         const amrex_real* src, const int* slo, const int* shi, const int* sblo,
                                         const int ncomp);

    FORT_LAUNCHABLE(void fort_fab_xpay)(const int* lo, const int* hi,
                                        amrex_real* dst, const int* dlo, const int* dhi,
                                        const amrex_real a,
                                        const amrex_real* src, const int* slo, const int* shi, const int* sblo,
                                        const int ncomp);

    FORT_LAUNCHABLE(void fort_fab_lincomb)(const int* lo, const int* hi,
                                           amrex_real* dst, const int* dlo, const int* dhi,
                                           const amrex_real a, const amrex_real* x, const int* xlo, const int* xhi, const int* xblo,
                                           const amrex_real b, const amrex_real* y, const int* ylo, const int* yhi, const int* yblo,
                                           const int ncomp);

    FORT_LAUNCHABLE(void fort_fab_addproduct)(const int* lo, const int* hi,
                                              amrex_real* dst, const int* dlo, const int* dhi,
                                              const amrex_real* src1, const int* s1lo, const int* s1hi,
                                              const amrex_real* src2, const int* s2lo, const int* s2hi,
                                              const int ncomp);

    FORT_LAUNCHABLE(void fort_fab_dot)(const int* lo, const int* hi,
                                       const amrex_real* x, const int* xlo, const int* xhi,
                                       const amrex_real* y, const int* ylo, const int* yhi, const int* yblo,
                                       const int ncomp, amrex_real* dp);

    amrex_real fort_fab_norminfmask (const int* lo, const int* hi,
                                     const int* msk, const int* mlo, const int* mhi,
                                     const amrex_real* src, const int* slo, const int* shi, const int* ncomp);

    amrex_real fort_fab_dot_mask (const int* lo, const int* hi,
		       const amrex_real* x, const int* xlo, const int* xhi,
		       const amrex_real* y, const int* ylo, const int* yhi, const int* yblo,
		       const int* m, const int* mlo, const int* mhi,
		       const int* ncomp);

    void fort_ifab_copy (const int* lo, const int* hi,
                         int* dst, const int* dlo, const int* dhi,
                         const int* src, const int* slo, const int* shi, const int* sblo,
                         const int* ncomp);

    long fort_ifab_copytomem (const int* lo, const int* hi,
                              int* dst,
                              const int* src, const int* slo, const int* shi,
                              const int* ncomp);

    long fort_ifab_copyfrommem (const int* lo, const int* hi,
                                const int* dst, const int* dlo, const int* dhi, const int* ncomp,
                                const int* src);
    
    void fort_ifab_setval (const int* lo, const int* hi, 
                           const int* dst, const int* dlo, const int* dhi, const int* ncomp,
                           const int* val);

    void fort_ifab_plus (const int* lo, const int* hi,
                         int* dst, const int* dlo, const int* dhi,
                         const int* src, const int* slo, const int* shi, const int* sblo,
                         const int* ncomp);    

    void fort_ifab_minus (const int* lo, const int* hi,
                          int* dst, const int* dlo, const int* dhi,
                          const int* src, const int* slo, const int* shi, const int* sblo,
                          const int* ncomp);    

    void amrex_fab_setval_ifnot (const int* lo, const int* hi,
                                 amrex_real* dst, const int* dlo, const int* dhi, const int* ncomp,
                                 const int* msk, const int* mlo, const int* mhi,
                                 const int val);

#ifdef __cplusplus
}
#endif

#endif
