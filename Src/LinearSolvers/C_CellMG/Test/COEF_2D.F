#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include "COEF_F.H"
#include "ArrayLim.H"

c-----------------------------------------------------------------------
      subroutine FORT_SET_CC_COEF (
     $     coef, DIMS(coef), lo, hi, dx)

      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(coef)
      REAL_T coef(DIMV(coef))
      REAL_T  dx(BL_SPACEDIM)
c
      integer i,j
      REAL_T  x,y,r
      REAL_T  xcent,ycent

      xcent = 0.25d0
      ycent = 0.25d0
c
      do j = lo(2)-1, hi(2)+1
         y = (dble(j)-lo(2)+0.5d0)*dx(2) - ycent
         do i = lo(1)-1, hi(1)+1
            x = (dble(i)-lo(1)+0.5d0)*dx(1) - xcent
            r = sqrt(x**2 + y**2) / 0.1d0
            coef(i,j) = 1.0d0 + 1000.d0 * (0.5d0 * (1.0d0 + tanh((2.d0-r))))
         end do
      end do
c
      end

c-----------------------------------------------------------------------
      subroutine FORT_COEF_TO_EDGES (dir, edge, DIMS(edge),
     $                               coef, DIMS(coef), lo, hi)

      integer dir
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(edge)
      REAL_T edge(DIMV(edge))
      integer DIMDEC(coef)
      REAL_T coef(DIMV(coef))
c
      integer i,j
c
      if (dir .eq. 0) then
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
            edge(i,j) = 0.5d0 * (coef(i,j) + coef(i-1,j))
         end do
        end do
      else if (dir .eq. 1) then
        do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
            edge(i,j) = 0.5d0 * (coef(i,j) + coef(i,j-1))
         end do
        end do
      end if
c
      end

