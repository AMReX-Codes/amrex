#ifndef _PROB_AMR_F_H_
#define _PROB_AMR_F_H_

/*
** $Id: PROB_AMR_F.H,v 1.2 1997-12-14 23:34:29 lijewski Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_PROBINIT    probinit
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_PROBINIT    PROBINIT
#  else
#    define FORT_PROBINIT    probinit_
#  endif

extern "C"
{
     void FORT_PROBINIT (const int* init);
};

#endif

#endif /*_PROB_AMR_F_H_*/
