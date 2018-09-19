#ifndef AMREX_BASEFAB_C_H_
#define AMREX_BASEFAB_C_H_

#include <AMReX_REAL.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_c_fab_copy (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                           amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                           const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT sblo,
                           const int* AMREX_RESTRICT ncomp);

    long amrex_c_fab_copytomem (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                                amrex_real* AMREX_RESTRICT dst,
                                const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                                const int* AMREX_RESTRICT ncomp);

    long amrex_c_fab_copyfrommem (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                                  amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi, const int* AMREX_RESTRICT ncomp,
                                  const amrex_real* AMREX_RESTRICT src);

    void amrex_c_fab_setval (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi, 
                             amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi, const int* AMREX_RESTRICT ncomp,
                             const amrex_real* AMREX_RESTRICT val);

    void amrex_c_fab_setval_ifnot (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                                   amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi, const int* AMREX_RESTRICT ncomp,
                                   const int* AMREX_RESTRICT msk, const int* AMREX_RESTRICT mlo, const int* AMREX_RESTRICT mhi,
                                   const amrex_real* AMREX_RESTRICT val);

    amrex_real amrex_c_fab_norminfmask (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                                        const int* AMREX_RESTRICT msk, const int* AMREX_RESTRICT mlo, const int* AMREX_RESTRICT mhi,
                                        const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT ncomp);

    amrex_real amrex_c_fab_norm (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                                 const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT ncomp,
                                 const int* AMREX_RESTRICT p);

    amrex_real amrex_c_fab_sum (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                                const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT ncomp);
    
    void amrex_c_fab_plus (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                           amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                           const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT sblo,
                           const int* AMREX_RESTRICT ncomp);    

    void amrex_c_fab_minus (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                            amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                            const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT sblo,
                            const int* AMREX_RESTRICT ncomp);    
    
    void amrex_c_fab_mult (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                           amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                           const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT sblo,
                           const int* AMREX_RESTRICT ncomp);    

    void amrex_c_fab_divide (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                             amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                             const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT sblo,
                             const int* AMREX_RESTRICT ncomp);    

    void amrex_c_fab_protdivide (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                                 amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                                 const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT sblo,
                                 const int* AMREX_RESTRICT ncomp);        
    
    void amrex_c_fab_invert (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                             amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi, const int* AMREX_RESTRICT ncomp,
                             const amrex_real* AMREX_RESTRICT a);

    void amrex_c_fab_saxpy (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                            amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                            const amrex_real* AMREX_RESTRICT a,
                            const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT sblo,
                            const int* AMREX_RESTRICT ncomp);
    
    void amrex_c_fab_xpay (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                           amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                           const amrex_real* AMREX_RESTRICT a,
                           const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT sblo,
                           const int* AMREX_RESTRICT ncomp);
    
    void amrex_c_fab_lincomb (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                              amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                              const amrex_real* AMREX_RESTRICT a, const amrex_real* AMREX_RESTRICT x, const int* AMREX_RESTRICT xlo, const int* AMREX_RESTRICT xhi, const int* AMREX_RESTRICT xblo,
                              const amrex_real* AMREX_RESTRICT b, const amrex_real* AMREX_RESTRICT y, const int* AMREX_RESTRICT ylo, const int* AMREX_RESTRICT yhi, const int* AMREX_RESTRICT yblo,
                              const int* AMREX_RESTRICT ncomp);
    
    void amrex_c_fab_addproduct (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                                 amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                                 const amrex_real* AMREX_RESTRICT src1, const int* AMREX_RESTRICT s1lo, const int* AMREX_RESTRICT s1hi,
                                 const amrex_real* AMREX_RESTRICT src2, const int* AMREX_RESTRICT s2lo, const int* AMREX_RESTRICT s2hi,
                                 const int* AMREX_RESTRICT ncomp);

    amrex_real amrex_c_fab_dot (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                                const amrex_real* AMREX_RESTRICT x, const int* AMREX_RESTRICT xlo, const int* AMREX_RESTRICT xhi,
                                const amrex_real* AMREX_RESTRICT y, const int* AMREX_RESTRICT ylo, const int* AMREX_RESTRICT yhi, const int* AMREX_RESTRICT yblo,
                                const int* AMREX_RESTRICT ncomp);
    
    amrex_real amrex_c_fab_dot_mask (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                                     const amrex_real* AMREX_RESTRICT x, const int* AMREX_RESTRICT xlo, const int* AMREX_RESTRICT xhi,
                                     const amrex_real* AMREX_RESTRICT y, const int* AMREX_RESTRICT ylo, const int* AMREX_RESTRICT yhi, const int* AMREX_RESTRICT yblo,
                                     const int* AMREX_RESTRICT m, const int* AMREX_RESTRICT mlo, const int* AMREX_RESTRICT mhi,
                                     const int* AMREX_RESTRICT ncomp);

    void amrex_c_ifab_copy (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                            int* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                            const int* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT sblo,
                            const int* AMREX_RESTRICT ncomp);

    long amrex_c_ifab_copytomem (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                                 int* AMREX_RESTRICT dst,
                                 const int* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                                 const int* AMREX_RESTRICT ncomp);

    long amrex_c_ifab_copyfrommem (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                                   int* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi, const int* AMREX_RESTRICT ncomp,
                                   const int* AMREX_RESTRICT src);
    
    void amrex_c_ifab_setval (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi, 
                              int* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi, const int* AMREX_RESTRICT ncomp,
                              const int* AMREX_RESTRICT val);

    void amrex_c_ifab_plus (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                            int* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                            const int* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT sblo,
                            const int* AMREX_RESTRICT ncomp);    
    
    void amrex_c_ifab_minus (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                             int* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                             const int* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT sblo,
                             const int* AMREX_RESTRICT ncomp);    
    
#ifdef __cplusplus
}
#endif

#endif
