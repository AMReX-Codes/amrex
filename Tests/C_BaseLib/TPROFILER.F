#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "TPROFILER_F.H"


c     ----------------------------------------------------------
      subroutine FORT_TESTPROFILER ()
        integer  i
        print *, "_in FORT_TESTPROFILER"
        do i = 1, 12000
c         call bl_proffortfuncstart("testprofiler")
c         call bl_proffortfuncstop("testprofiler")
          call bl_proffortfuncstart("t")
          call bl_proffortfuncstop("t")
        enddo
        print *, "_out FORT_TESTPROFILER"
      end

c     ----------------------------------------------------------
      subroutine FORT_TESTPROFILERINT ()
        integer  i
        print *, "_in FORT_TESTPROFILERINT"
        do i = 1, 12000
          call bl_proffortfuncstart_int(3)
          call bl_proffortfuncstop_int(3)
        enddo
        print *, "_out FORT_TESTPROFILERINT"
      end
