
#ifndef _MACOPERATOR_F_H_
#define _MACOPERATOR_F_H_

#ifdef BL_LANG_FORT
#    define FORT_MACCOEF       maccoef
#    define FORT_MACRHS        macrhs
#    define FORT_MACUPDATE     macupdate
#    define FORT_MACSYNCRHS    macsyncrhs
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_MACCOEF       MACCOEF
#    define FORT_MACRHS        MACRHS
#    define FORT_MACUPDATE     MACUPDATE
#    define FORT_MACSYNCRHS    MACSYNCRHS
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_MACCOEF       maccoef
#    define FORT_MACRHS        macrhs
#    define FORT_MACUPDATE     macupdate
#    define FORT_MACSYNCRHS    macsyncrhs
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_MACCOEF       maccoef_
#    define FORT_MACRHS        macrhs_
#    define FORT_MACUPDATE     macupdate_
#    define FORT_MACSYNCRHS    macsyncrhs_
#endif

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C"
{
#endif
#if (BL_SPACEDIM == 2)
    void FORT_MACCOEF (amrex_real* cx, AMREX_ARLIM_P(cxlo), AMREX_ARLIM_P(cxhi),
                       amrex_real* cy, AMREX_ARLIM_P(cylo), AMREX_ARLIM_P(cyhi),
                       const amrex_real* ax, AMREX_ARLIM_P(axlo), AMREX_ARLIM_P(axhi),
                       const amrex_real* ay, AMREX_ARLIM_P(aylo), AMREX_ARLIM_P(ayhi),
                       const amrex_real* rho, AMREX_ARLIM_P(rlo), AMREX_ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const amrex_real* dx);

    void FORT_MACRHS  (const amrex_real* ux, AMREX_ARLIM_P(uxlo), AMREX_ARLIM_P(uxhi),
                       const amrex_real* uy, AMREX_ARLIM_P(uylo), AMREX_ARLIM_P(uyhi),
                       const amrex_real* ax, AMREX_ARLIM_P(axlo), AMREX_ARLIM_P(axhi),
                       const amrex_real* ay, AMREX_ARLIM_P(aylo), AMREX_ARLIM_P(ayhi),
                       const amrex_real* vol, AMREX_ARLIM_P(vlo), AMREX_ARLIM_P(vhi),
                       amrex_real* rhs, AMREX_ARLIM_P(rlo), AMREX_ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const amrex_real* scale);

    void FORT_MACUPDATE (const int *init,
                         amrex_real* ux, AMREX_ARLIM_P(uxlo), AMREX_ARLIM_P(uxhi),
                         amrex_real* uy, AMREX_ARLIM_P(uylo), AMREX_ARLIM_P(uyhi),
                         const amrex_real* phi, AMREX_ARLIM_P(p_lo), AMREX_ARLIM_P(p_hi),
                         const amrex_real* rho, AMREX_ARLIM_P(rlo), AMREX_ARLIM_P(rhi),
                         const int* lo, const int* hi,
                         const amrex_real* dx, const amrex_real* scale);
#endif

#if (BL_SPACEDIM == 3)
    void FORT_MACCOEF (amrex_real* cx, AMREX_ARLIM_P(cxlo), AMREX_ARLIM_P(cxhi),
                       amrex_real* cy, AMREX_ARLIM_P(cylo), AMREX_ARLIM_P(cyhi),
                       amrex_real* cz, AMREX_ARLIM_P(czlo), AMREX_ARLIM_P(czhi),
                       const amrex_real* ax, AMREX_ARLIM_P(axlo), AMREX_ARLIM_P(axhi),
                       const amrex_real* ay, AMREX_ARLIM_P(aylo), AMREX_ARLIM_P(ayhi),
                       const amrex_real* az, AMREX_ARLIM_P(azlo), AMREX_ARLIM_P(azhi),
                       const amrex_real* rho, AMREX_ARLIM_P(rlo), AMREX_ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const amrex_real* dx);

    void FORT_MACRHS  (const amrex_real* ux, AMREX_ARLIM_P(uxlo), AMREX_ARLIM_P(uxhi),
                       const amrex_real* uy, AMREX_ARLIM_P(uylo), AMREX_ARLIM_P(uyhi),
                       const amrex_real* uz, AMREX_ARLIM_P(uzlo), AMREX_ARLIM_P(uzhi),
                       const amrex_real* ax, AMREX_ARLIM_P(axlo), AMREX_ARLIM_P(axhi),
                       const amrex_real* ay, AMREX_ARLIM_P(aylo), AMREX_ARLIM_P(ayhi),
                       const amrex_real* az, AMREX_ARLIM_P(azlo), AMREX_ARLIM_P(azhi),
                       const amrex_real* vol, AMREX_ARLIM_P(vlo), AMREX_ARLIM_P(vhi),
                       amrex_real* rhs, AMREX_ARLIM_P(rlo), AMREX_ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const amrex_real* scale);

    void FORT_MACUPDATE (const int *init,
                         amrex_real* ux, AMREX_ARLIM_P(uxlo), AMREX_ARLIM_P(uxhi),
                         amrex_real* uy, AMREX_ARLIM_P(uylo), AMREX_ARLIM_P(uyhi),
                         amrex_real* uz, AMREX_ARLIM_P(uzlo), AMREX_ARLIM_P(uzhi),
                         const amrex_real* phi, AMREX_ARLIM_P(p_lo), AMREX_ARLIM_P(p_hi),
                         const amrex_real* rho, AMREX_ARLIM_P(rlo), AMREX_ARLIM_P(rhi),
                         const int* lo, const int* hi,
                         const amrex_real* dx, const amrex_real* scale); 

#endif

    void FORT_MACSYNCRHS (amrex_real* rhs, AMREX_ARLIM_P(rlo), AMREX_ARLIM_P(rhi),
                          const int* lo, const int* hi,
                          const amrex_real* vol, AMREX_ARLIM_P(vlo), AMREX_ARLIM_P(vhi),
                          const amrex_real* rhs_scale);
#ifdef __cplusplus
}
#endif
#endif

#endif /*_MACOPERATOR_F_H_*/
