#ifndef _PROB_F_H_
#define _PROB_F_H_

#ifdef BL_LANG_FORT
#    define FORT_PROBINIT    probinit
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_PROBINIT    PROBINIT
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_PROBINIT    probinit
#  else
#    define FORT_PROBINIT    probinit_
#  endif
#  endif

#include "ArrayLim.H"

extern "C"
{
     void FORT_PROBINIT (const int* init,
                         const int* name,
                         const int* namelen,
                         const Real* problo,
                         const Real* probhi);

}
#endif

#endif /*_PROB_F_H_*/
