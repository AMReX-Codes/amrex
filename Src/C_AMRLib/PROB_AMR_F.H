#ifndef _PROB_AMR_F_H_
#define _PROB_AMR_F_H_

/*
** $Id: PROB_AMR_F.H,v 1.5 1998-10-07 21:18:20 vince Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_PROBINIT    probinit
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_PROBINIT    PROBINIT
#  else
#  ifdef BL_AIX
#    define FORT_PROBINIT    probinit
#  else
#    define FORT_PROBINIT    probinit_
#  endif
#  endif

extern "C"
{
     void FORT_PROBINIT (const int* init,
                         const int* name,
                         const int* namelen);
}

#endif

#endif /*_PROB_AMR_F_H_*/
