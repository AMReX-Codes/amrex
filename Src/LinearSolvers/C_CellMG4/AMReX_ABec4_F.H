#ifndef _AABBEC4_F_H_
#define _AABBEC4_F_H_

#include <AMReX_REAL.H>

#if        defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 1)
#define FORT_ADOTX         adotx1daabbec4
#define FORT_FLUX          flux1daabbec4
#define FORT_APPLYBC4      applybc1daabbec4
#define FORT_CA2CC         ca2cc1daabbec4
#define FORT_CC2CA         cc2ca1daabbec4
#define FORT_LO_CC2EC      cc2ec1daabbec4
#endif

#if (BL_SPACEDIM == 2)
#define FORT_ADOTX         adotx2daabbec4
#define FORT_FLUX          flux2daabbec4
#define FORT_APPLYBC4      applybc2daabbec4
#define FORT_APPLYBC4_TOUCHUP abc4tu
#define FORT_CA2CC         ca2cc2daabbec4
#define FORT_CC2CA         cc2ca2daabbec4
#define FORT_LO_CC2EC      cc2ec2daabbec4
#endif

#if (BL_SPACEDIM == 3)
#define FORT_ADOTX         adotx3daabbec4
#define FORT_FLUX          flux3daabbec4
#define FORT_APPLYBC4      applybc3daabbec4
#define FORT_APPLYBC4_TOUCHUP abc4tu
#define FORT_CA2CC         ca2cc3daabbec4
#define FORT_CC2CA         cc2ca3daabbec4
#define FORT_LO_CC2EC      cc2ec3daabbec4
#endif

#else

#if (BL_SPACEDIM == 1)
 
#if  defined(BL_FORT_USE_UPPERCASE)
#define FORT_ADOTX         ADOTX1DAABBEC4
#define FORT_FLUX          FLUX1DAABBEC4
#define FORT_APPLYBC4      APPLYBC1DAABBEC4
#define FORT_APPLYBC4_TOUCHUP ABC4TU
#define FORT_CA2CC         CA2CC1DAABBEC4
#define FORT_CC2CA         CC2CA1DAABBEC4
#define FORT_LO_CC2EC      CC2EC1DAABBEC4
#elif defined(BL_FORT_USE_LOWERCASE)
#define FORT_ADOTX         adotx1daabbec4
#define FORT_FLUX          flux1daabbec4
#define FORT_APPLYBC4      applybc1daabbec4
#define FORT_APPLYBC4_TOUCHUP abc4tu
#define FORT_CA2CC         ca2cc1daabbec4
#define FORT_CC2CA         cc2ca1daabbec4
#define FORT_LO_CC2EC      cc2ec1daabbec4
#elif defined(BL_FORT_USE_UNDERSCORE)
#define FORT_ADOTX         adotx1daabbec4_
#define FORT_FLUX          flux1daabbec4_
#define FORT_APPLYBC4      applybc1daabbec4_
#define FORT_APPLYBC4_TOUCHUP abc4tu_
#define FORT_CA2CC         ca2cc1daabbec4_
#define FORT_CC2CA         cc2ca1daabbec4_
#define FORT_LO_CC2EC      cc2ec1daabbec4_
#endif
 
#endif

#if (BL_SPACEDIM == 2)

#if  defined(BL_FORT_USE_UPPERCASE)
#define FORT_ADOTX    ADOTX2DAABBEC4
#define FORT_FLUX     FLUX2DAABBEC4
#define FORT_APPLYBC4 APPLYBC2DAABBEC4
#define FORT_APPLYBC4_TOUCHUP ABC4TU
#define FORT_CA2CC    CA2CC2DAABBEC4
#define FORT_CC2CA    CC2CA2DAABBEC4
#define FORT_LO_CC2EC  CC2EC2DAABBEC4
#elif defined(BL_FORT_USE_LOWERCASE)
#define FORT_ADOTX    adotx2daabbec4
#define FORT_FLUX     flux2daabbec4
#define FORT_APPLYBC4 applybc2daabbec4
#define FORT_APPLYBC4_TOUCHUP abc4tu
#define FORT_CA2CC    ca2cc2daabbec4
#define FORT_CC2CA    cc2ca2daabbec4
#define FORT_LO_CC2EC  cc2ec2daabbec4
#elif defined(BL_FORT_USE_UNDERSCORE)
#define FORT_ADOTX    adotx2daabbec4_
#define FORT_FLUX     flux2daabbec4_
#define FORT_APPLYBC4 applybc2daabbec4_
#define FORT_APPLYBC4_TOUCHUP abc4tu_
#define FORT_CA2CC    ca2cc2daabbec4_
#define FORT_CC2CA    cc2ca2daabbec4_
#define FORT_LO_CC2EC  cc2ec2daabbec4_
#endif

#endif


#if (BL_SPACEDIM == 3)

#if   defined(BL_FORT_USE_UPPERCASE)
#define FORT_ADOTX    ADOTX3DAABBEC4
#define FORT_FLUX     FLUX3DAABBEC4
#define FORT_APPLYBC4 APPLYBC3DAABBEC4
#define FORT_APPLYBC4_TOUCHUP ABC4TU
#define FORT_CA2CC    CA2CC3DAABBEC4
#define FORT_CC2CA    CC2CA3DAABBEC4
#define FORT_LO_CC2EC  CC2EC3DAABBEC4
#elif defined(BL_FORT_USE_LOWERCASE)
#define FORT_ADOTX    adotx3daabbec4
#define FORT_APPLYBC4 applybc3daabbec4
#define FORT_APPLYBC4_TOUCHUP abc4tu
#define FORT_FLUX     flux3daabbec4
#define FORT_CA2CC    ca2cc3daabbec4
#define FORT_CC2CA    cc2ca3daabbec4
#define FORT_LO_CC2EC  cc2ec3daabbec4
#elif defined(BL_FORT_USE_UNDERSCORE)
#define FORT_ADOTX    adotx3daabbec4_
#define FORT_FLUX     flux3daabbec4_
#define FORT_APPLYBC4 applybc3daabbec4_
#define FORT_APPLYBC4_TOUCHUP abc4tu_
#define FORT_CA2CC    ca2cc3daabbec4_
#define FORT_CC2CA    cc2ca3daabbec4_
#define FORT_LO_CC2EC  cc2ec3daabbec4_
#endif

#endif

#include <AMReX_ArrayLim.H>

extern "C"
{
    void FORT_APPLYBC4(
        const int *flagden, const int *flagbc, const int *maxorder,
        amrex::Real *phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const int *cdr,
        const int *bct,
        const amrex::Real *bcl,
        const amrex::Real *bcval, ARLIM_P(bcval_lo), ARLIM_P(bcval_hi),
        const int *mask,   ARLIM_P(mask_lo),  ARLIM_P(mask_hi), 
        amrex::Real *den,         ARLIM_P(den_lo),   ARLIM_P(den_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex::Real *h
        );

    void FORT_APPLYBC4_TOUCHUP(
        amrex::Real *phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const int *lo, const int *hi, const int *nc);

    void FORT_ADOTX(
        amrex::Real *y      , ARLIM_P(y_lo), ARLIM_P(y_hi),
        const amrex::Real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const amrex::Real* alpha, const amrex::Real* beta,
        const amrex::Real* a, ARLIM_P(a_lo), ARLIM_P(a_hi),
        const amrex::Real* b, ARLIM_P(b_lo), ARLIM_P(b_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex::Real *h
        );
    
    void FORT_FLUX(
        const amrex::Real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const amrex::Real* alpha, const amrex::Real* beta,
        const amrex::Real* a, ARLIM_P(a_lo), ARLIM_P(a_hi),
        const amrex::Real* b, ARLIM_P(b_lo), ARLIM_P(b_hi),
        const int *nc,
        const amrex::Real *h,
        const int *xlo, const int *xhi,
        const amrex::Real* xflux, ARLIM_P(xflux_lo), ARLIM_P(xflux_hi)
#if (BL_SPACEDIM >= 2)
        ,const int *ylo, const int *yhi
        ,const amrex::Real* yflux, ARLIM_P(yflux_lo), ARLIM_P(yflux_hi)
#endif
#if (BL_SPACEDIM == 3)
        ,const int *zlo, const int *zhi
        ,const amrex::Real* zflux, ARLIM_P(zflux_lo), ARLIM_P(zflux_hi)
#endif
        );

      void FORT_CA2CC(const int* lo, const int* hi,
                      const amrex::Real* ca, ARLIM_P(ca_lo), ARLIM_P(ca_hi),
                      amrex::Real*       cc, ARLIM_P(cc_lo), ARLIM_P(cc_hi),
                      const int* nc);
      void FORT_CC2CA(const int* lo, const int* hi,
                      const amrex::Real* cc, ARLIM_P(cc_lo), ARLIM_P(cc_hi),
                      amrex::Real*       ca, ARLIM_P(ca_lo), ARLIM_P(ca_hi),
                      const int* nc);
      void FORT_LO_CC2EC(const int* lo, const int* hi,
                         const amrex::Real* cc, ARLIM_P(cc_lo), ARLIM_P(cc_hi),
                         amrex::Real*       ca, ARLIM_P(ca_lo), ARLIM_P(ca_hi),
                         const int* nc, const int* dir, const int* isharm);
}
#endif

#endif /*_AABBEC4_F_H_*/
