#ifndef _AMRPROFTEST_F_H 
#define _AMRPROFTEST_F_H
               
#ifdef BL_LANG_FORT
#    define FORT_AMRPROFTEST      amrproftest
#    define FORT_AMRPROFTESTINT    amrproftestint
#else 
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_AMRPROFTEST      AMRPROFTEST
#    define FORT_AMRPROFTESTINT    AMRPROFTESTNT
#  else
#    define FORT_AMRPROFTEST      amrproftest_
#    define FORT_AMRPROFTESTINT    amrproftestint_ 
#  endif


#ifdef __cplusplus
extern "C" {
#endif
    void FORT_AMRPROFTEST();
    void FORT_AMRPROFTESTINT();
#ifdef __cplusplus
}     
#endif
 
#endif
#endif
