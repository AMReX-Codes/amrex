
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

#include <ArrayLim.H>

extern "C"
{
#if (BL_SPACEDIM == 2)
    void FORT_MACCOEF (Real* cx, ARLIM_P(cxlo), ARLIM_P(cxhi),
                       Real* cy, ARLIM_P(cylo), ARLIM_P(cyhi),
                       const Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const Real* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const Real* dx);

    void FORT_MACRHS  (const Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                       const Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                       const Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       Real* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const Real* scale);

    void FORT_MACUPDATE (const int *init,
                         Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                         Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                         const Real* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                         const Real* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                         const int* lo, const int* hi,
                         const Real* dx, const Real* scale);
#endif

#if (BL_SPACEDIM == 3)
    void FORT_MACCOEF (Real* cx, ARLIM_P(cxlo), ARLIM_P(cxhi),
                       Real* cy, ARLIM_P(cylo), ARLIM_P(cyhi),
                       Real* cz, ARLIM_P(czlo), ARLIM_P(czhi),
                       const Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const Real* az, ARLIM_P(azlo), ARLIM_P(azhi),
                       const Real* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const Real* dx);

    void FORT_MACRHS  (const Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                       const Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                       const Real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
                       const Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),
                       const Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
                       const Real* az, ARLIM_P(azlo), ARLIM_P(azhi),
                       const Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       Real* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
                       const int* lo, const int* hi,
                       const Real* scale);

    void FORT_MACUPDATE (const int *init,
                         Real* ux, ARLIM_P(uxlo), ARLIM_P(uxhi),
                         Real* uy, ARLIM_P(uylo), ARLIM_P(uyhi),
                         Real* uz, ARLIM_P(uzlo), ARLIM_P(uzhi),
                         const Real* phi, ARLIM_P(p_lo), ARLIM_P(p_hi),
                         const Real* rho, ARLIM_P(rlo), ARLIM_P(rhi),
                         const int* lo, const int* hi,
                         const Real* dx, const Real* scale); 

#endif

    void FORT_MACSYNCRHS (Real* rhs, ARLIM_P(rlo), ARLIM_P(rhi),
                          const int* lo, const int* hi,
                          const Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                          const Real* rhs_scale);
}
#endif

#endif /*_MACOPERATOR_F_H_*/
