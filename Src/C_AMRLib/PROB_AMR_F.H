#ifndef _PROB_AMR_F_H_
#define _PROB_AMR_F_H_

/*
** $Id: PROB_AMR_F.H,v 1.4 1998-05-12 22:12:43 lijewski Exp $
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
     void FORT_PROBINIT (const int* init,
                         const int* name,
                         const int* namelen);
}

#endif

#endif /*_PROB_AMR_F_H_*/
