
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include "AMReX_MG_F.H"
#include "AMReX_ArrayLim.H"

    subroutine FORT_AVERAGE ( &
           c, c_l1,c_l2,c_h1,c_h2, &
           f, f_l1,f_l2,f_h1,f_h2, &
           lo, hi, nc)

      implicit none

      integer nc
      integer f_l1,f_l2,f_h1,f_h2
      integer c_l1,c_l2,c_h1,c_h2
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)

      integer i
      integer j
      integer n
      REAL_T denom
      parameter(denom=fourth)

      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               c(i,j,n) =  ( &
                    f(2*i+1,2*j+1,n) + f(2*i  ,2*j+1,n) &
                    + f(2*i+1,2*j,n ) + f(2*i  ,2*j ,n))*denom
            end do
         end do
      end do

    end subroutine FORT_AVERAGE

    subroutine FORT_INTERP ( &
           f, f_l1,f_l2,f_h1,f_h2, &
           c, c_l1,c_l2,c_h1,c_h2, &
           lo, hi, nc)

      implicit none

      integer nc
      integer f_l1,f_l2,f_h1,f_h2
      integer c_l1,c_l2,c_h1,c_h2
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)

      integer i, j, n, twoi, twoj, twoip1, twojp1

!     MultiGrid::relax(...) does only V-cycles (not F-cycles), and for V-cycles, 
!     piecewise-constant interpolation performs better than linear interpolation,
!     as measured both by run-time and number of V-cycles for convergence.

      do n = 1, nc
         do j = lo(2),hi(2)
            twoj   = 2*j
            twojp1 = twoj+1

            do i = lo(1),hi(1)

               twoi   = 2*i
               twoip1 = twoi+1

               f(twoi,   twoj  ,n) = f(twoi,   twoj  ,n) + c(i,j,n)
               f(twoip1, twoj  ,n) = f(twoip1, twoj  ,n) + c(i,j,n)
               f(twoi,   twojp1,n) = f(twoi,   twojp1,n) + c(i,j,n)
               f(twoip1, twojp1,n) = f(twoip1, twojp1,n) + c(i,j,n)

            end do
         end do
      end do

    end subroutine FORT_INTERP
