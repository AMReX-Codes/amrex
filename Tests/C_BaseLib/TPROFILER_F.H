#ifndef _AVGDOWN_F_H 
#define _AVGDOWN_F_H
               
#ifdef BL_LANG_FORT
#    define FORT_TESTPROFILER         testprofiler
#    define FORT_TESTPROFILERINT         testprofilerint
#else 
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_TESTPROFILER         TESTPROFILER
#    define FORT_TESTPROFILERINT         TESTPROFILERINT
#  else
#    define FORT_TESTPROFILER         testprofiler_ 
#    define FORT_TESTPROFILERINT         testprofilerint_ 
#  endif


extern "C" {
    void FORT_TESTPROFILER ();
    void FORT_TESTPROFILERINT ();
}     
 
#endif
#endif
