#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "PROB_F.H"
#include "ArrayLim.H"

#define SDIM 3


      subroutine FORT_PROBINIT (init,name,namlen,problo,probhi)
      integer init,namlen
      integer name(namlen)
      REAL_T  problo(SDIM), probhi(SDIM)

      return
      end
