#ifndef _SLABSTAT_F_H_
#define _SLABSTAT_F_H_

/*
** $Id: SLABSTAT_F.H,v 1.3 2002-04-26 22:10:02 lijewski Exp $
*/

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

extern "C"
{

     void FORT_REYNAVG_FIRST4MOMS (const Real* src, 
                                     ARLIM_P(src_lo), ARLIM_P(src_hi),
                                   const int* nsrc,
                                   Real* dst, ARLIM_P(dst_lo), ARLIM_P(dst_hi),
                                   const int* ndst,
                                   const Real* dt,
                                   const Real* dx);

     void FORT_FAVREAVG_FIRST4MOMS (const Real* src, 
                                      ARLIM_P(src_lo), ARLIM_P(src_hi),
                                    const int* nsrc,
                                    Real* dst, ARLIM_P(dst_lo), ARLIM_P(dst_hi),
                                    const int* ndst,
                                    const Real* dt,
                                    const Real* dx);
}

#endif

#endif /*_SLABSTAT_F_H_*/
