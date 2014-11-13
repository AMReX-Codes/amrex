#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "TPROFILER_F.H"


c     ----------------------------------------------------------
      subroutine FORT_TESTPROFILER ()

      integer  i
      
      print *, "_in FORT_TESTPROFILER"
      do i = 1, 1000000
      call bl_proffortfuncstart("testprofiler")
      call bl_proffortfuncstop("testprofiler")
      enddo
      print *, "_out FORT_TESTPROFILER"

      end
