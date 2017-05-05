
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
    void FORT_MACCOEF (amrex_real* cx, ARLIM_P(cxlo), ARLIM_P(cxhi),
                       amrex_real* cy, ARLIM_P(cylo), ARLIM_P(cyhi),
                       const amrex_real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const amrex_real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const amrex_real* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const amrex_real* dx);

    void FORT_MACRHS  (const amrex_real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                       const amrex_real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                       const amrex_real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const amrex_real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const amrex_real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       amrex_real* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const amrex_real* scale);

    void FORT_MACUPDATE (const int *init,
                         amrex_real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                         amrex_real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                         const amrex_real* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                         const amrex_real* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                         const int* lo, const int* hi,
                         const amrex_real* dx, const amrex_real* scale);
#endif

#if (BL_SPACEDIM == 3)
    void FORT_MACCOEF (amrex_real* cx, ARLIM_P(cxlo), ARLIM_P(cxhi),
                       amrex_real* cy, ARLIM_P(cylo), ARLIM_P(cyhi),
                       amrex_real* cz, ARLIM_P(czlo), ARLIM_P(czhi),
                       const amrex_real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const amrex_real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const amrex_real* az, ARLIM_P(azlo), ARLIM_P(azhi),
                       const amrex_real* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const amrex_real* dx);

    void FORT_MACRHS  (const amrex_real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                       const amrex_real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                       const amrex_real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
                       const amrex_real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const amrex_real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const amrex_real* az, ARLIM_P(azlo), ARLIM_P(azhi),
                       const amrex_real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       amrex_real* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const amrex_real* scale);

    void FORT_MACUPDATE (const int *init,
                         amrex_real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                         amrex_real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                         amrex_real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
                         const amrex_real* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                         const amrex_real* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                         const int* lo, const int* hi,
                         const amrex_real* dx, const amrex_real* scale); 

#endif

    void FORT_MACSYNCRHS (amrex_real* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
                          const int* lo, const int* hi,
                          const amrex_real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                          const amrex_real* rhs_scale);
#ifdef __cplusplus
}
#endif
#endif

#endif /*_MACOPERATOR_F_H_*/
