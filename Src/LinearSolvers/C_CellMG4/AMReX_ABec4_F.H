#ifndef AMREX_AABBEC4_F_H_
#define AMREX_AABBEC4_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void amrex_ab4_applybc4(
        const int *flagden, const int *flagbc, const int *maxorder,
        amrex_real *phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const int *cdr,
        const int *bct,
        const amrex_real *bcl,
        const amrex_real *bcval, ARLIM_P(bcval_lo), ARLIM_P(bcval_hi),
        const int *mask,   ARLIM_P(mask_lo),  ARLIM_P(mask_hi), 
        amrex_real *den,         ARLIM_P(den_lo),   ARLIM_P(den_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );

    void amrex_ab4_applybc4_touchup(
        amrex_real *phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const int *lo, const int *hi, const int *nc);

    void amrex_ab4_adotx(
        amrex_real *y      , ARLIM_P(y_lo), ARLIM_P(y_hi),
        const amrex_real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a, ARLIM_P(a_lo), ARLIM_P(a_hi),
        const amrex_real* b, ARLIM_P(b_lo), ARLIM_P(b_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );
    
    void amrex_ab4_flux(
        const amrex_real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a, ARLIM_P(a_lo), ARLIM_P(a_hi),
        const amrex_real* b, ARLIM_P(b_lo), ARLIM_P(b_hi),
        const int *nc,
        const amrex_real *h,
        const int *xlo, const int *xhi,
        const amrex_real* xflux, ARLIM_P(xflux_lo), ARLIM_P(xflux_hi)
#if (BL_SPACEDIM >= 2)
        ,const int *ylo, const int *yhi
        ,const amrex_real* yflux, ARLIM_P(yflux_lo), ARLIM_P(yflux_hi)
#endif
#if (BL_SPACEDIM == 3)
        ,const int *zlo, const int *zhi
        ,const amrex_real* zflux, ARLIM_P(zflux_lo), ARLIM_P(zflux_hi)
#endif
        );

      void amrex_ab4_ca2cc(const int* lo, const int* hi,
                      const amrex_real* ca, ARLIM_P(ca_lo), ARLIM_P(ca_hi),
                      amrex_real*       cc, ARLIM_P(cc_lo), ARLIM_P(cc_hi),
                      const int* nc);
      void amrex_ab4_cc2ca(const int* lo, const int* hi,
                      const amrex_real* cc, ARLIM_P(cc_lo), ARLIM_P(cc_hi),
                      amrex_real*       ca, ARLIM_P(ca_lo), ARLIM_P(ca_hi),
                      const int* nc);
      void amrex_ab4_lo_cc2ec(const int* lo, const int* hi,
                         const amrex_real* cc, ARLIM_P(cc_lo), ARLIM_P(cc_hi),
                         amrex_real*       ca, ARLIM_P(ca_lo), ARLIM_P(ca_hi),
                         const int* nc, const int* dir, const int* isharm);
#ifdef __cplusplus
}
#endif

#endif /*_AABBEC4_F_H_*/
