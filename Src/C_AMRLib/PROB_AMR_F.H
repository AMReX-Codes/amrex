#ifndef _PROB_AMR_F_H_
#define _PROB_AMR_F_H_

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

#endif

