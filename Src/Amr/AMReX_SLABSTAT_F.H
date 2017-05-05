#ifndef _SLABSTAT_F_H_
#define _SLABSTAT_F_H_

#ifdef BL_LANG_FORT
#    define FORT_REYNAVG_FIRST4MOMS    reynavg_first4moms
#    define FORT_FAVREAVG_FIRST4MOMS   favreavg_first4moms
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_REYNAVG_FIRST4MOMS    REYNAVG_FIRST4MOMS
#    define FORT_FAVREAVG_FIRST4MOMS   FAVREAVG_FIRST4MOMS
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_REYNAVG_FIRST4MOMS    reynavg_first4moms
#    define FORT_FAVREAVG_FIRST4MOMS   favreavg_first4moms
#  else
#    define FORT_REYNAVG_FIRST4MOMS    reynavg_first4moms_
#    define FORT_FAVREAVG_FIRST4MOMS   favreavg_first4moms_
#  endif
#  endif

#include <AMReX_REAL.H>

#ifdef __cplusplus
extern "C"
{
#endif

     void FORT_REYNAVG_FIRST4MOMS (const amrex_real* src, 
                                     ARLIM_P(src_lo), ARLIM_P(src_hi),
                                   const int* nsrc,
                                   amrex_real* dst, ARLIM_P(dst_lo), ARLIM_P(dst_hi),
                                   const int* ndst,
                                   const amrex_real* dt,
                                   const amrex_real* dx);

     void FORT_FAVREAVG_FIRST4MOMS (const amrex_real* src, 
                                      ARLIM_P(src_lo), ARLIM_P(src_hi),
                                    const int* nsrc,
                                    amrex_real* dst, ARLIM_P(dst_lo), ARLIM_P(dst_hi),
                                    const int* ndst,
                                    const amrex_real* dt,
                                    const amrex_real* dx);
#ifdef __cplusplus
}
#endif

#endif

#endif /*_SLABSTAT_F_H_*/
