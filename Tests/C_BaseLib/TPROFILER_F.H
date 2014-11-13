#ifndef _AVGDOWN_F_H 
#define _AVGDOWN_F_H
               
#ifdef BL_LANG_FORT
#    define FORT_TESTPROFILER         testprofiler
#else 
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_TESTPROFILER         TESTPROFILER
#  else
#    define FORT_TESTPROFILER         testprofiler_ 
#  endif


extern "C" {
      
    void FORT_TESTPROFILER ();

}     
 
#endif
#endif
