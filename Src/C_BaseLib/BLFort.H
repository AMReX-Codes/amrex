#define _BLFORT_H_
#define _BLFORT_H_ 

#include <ArrayLim.H>

#ifndef BL_LANG_FORT
#if   defined(BL_FORT_USE_UPPERCASE)
#define __BL_FORT_NAME__(A,B) A
#elif defined(BL_FORT_USE_LOWERCASE)
#define __BL_FORT_NAME__(A,B) B
#elif defined(BL_FORT_USE_UNDERSCORE)
#define __BL_FORT_NAME__(A,B) B##_
#else
#error "One of BL_FORT_USE_{UPPERCASE,LOWERCASE,UNDERSCORE} must be defined"
#endif
#define BL_FORT_PROC_DECL(A,B) extern "C" void __BL_FORT_NAME__(A,B)
#define BL_FORT_PROC_CALL(A,B) __BL_FORT_NAME__(A,B)
#define BL_FORT_PROC_NAME(A,B) __BL_FORT_NAME__(A,B)

#define BL_FORT_FAB_ARG(A) Real* A##_fab, ARLIM_P(A##_lo), ARLIM_P(A##_hi)

#endif

#endif
