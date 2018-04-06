#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include "AMReX_MG_F.H"

    subroutine FORT_AVERAGE ( &
           c, c_l1,c_h1, &
           f, f_l1,f_h1, &
           lo, hi, nc)

      integer nc
      integer f_l1,f_h1
      integer c_l1,c_h1
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(f_l1:f_h1,nc)
      REAL_T c(c_l1:c_h1,nc)

      integer i,n

      do n = 1, nc
         do i = lo(1), hi(1)
            c(i,n) =  half * ( f(2*i+1,n) + f(2*i,n) )
         end do
      end do

    end subroutine FORT_AVERAGE

    subroutine FORT_INTERP ( &
           f, f_l1,f_h1, &
           c, c_l1,c_h1, &
           lo, hi, nc)

      integer nc
      integer f_l1,f_h1
      integer c_l1,c_h1
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(f_l1:f_h1,nc)
      REAL_T c(c_l1:c_h1,nc)

      integer i,n

      do n = 1, nc
         do i = lo(1), hi(1)
            f(2*i+1,n) = c(i,n) + f(2*i+1,n)
            f(2*i  ,n) = c(i,n) + f(2*i  ,n)
         end do
      end do

    end subroutine FORT_INTERP
