#ifndef _PROB_AMR_F_H_
#define _PROB_AMR_F_H_

/*
** $Id: PROB_AMR_F.H,v 1.8 2002-04-16 18:06:28 vince Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_PROBINIT    probinit
#ifdef BL_SYNC_RANTABLES
#    define FORT_SYNC_RANTABLES    syncrantables
#endif
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_PROBINIT    PROBINIT
#ifdef BL_SYNC_RANTABLES
#    define FORT_SYNC_RANTABLES    SYNCRANTABLES
#endif
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_PROBINIT    probinit
#ifdef BL_SYNC_RANTABLES
#    define FORT_SYNC_RANTABLES    syncrantables
#endif
#  else
#    define FORT_PROBINIT    probinit_
#ifdef BL_SYNC_RANTABLES
#    define FORT_SYNC_RANTABLES    syncrantables_
#endif
#  endif
#  endif

extern "C"
{
     void FORT_PROBINIT (const int* init,
                         const int* name,
                         const int* namelen,
                         const Real* problo,
                         const Real* probhi);

#ifdef BL_SYNC_RANTABLES
     void FORT_SYNC_RANTABLES (Real *RanPhaseTable,
                               Real *RanAmplTable,
                               const int* igetput);
#endif

}

#endif

#endif /*_PROB_AMR_F_H_*/
