#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include "AMReX_ArrayLim.H"
#include "COEF_F.H"

c-----------------------------------------------------------------------
      subroutine FORT_SET_CC_COEF (
     $     coef, DIMS(coef), lo, hi, dx, prob_lo, prob_hi)

      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(coef)
      REAL_T coef(DIMV(coef))
      REAL_T dx(BL_SPACEDIM)
      REAL_T prob_lo(BL_SPACEDIM)
      REAL_T prob_hi(BL_SPACEDIM)
c
      integer i,j,k
      REAL_T  x,y,z,r
      REAL_T  xcent,ycent,zcent
      REAL_T  xlen ,ylen ,zlen

      xlen = prob_hi(1)-prob_lo(1)
      ylen = prob_hi(2)-prob_lo(2)
      zlen = prob_hi(3)-prob_lo(3)

      xcent = prob_lo(1) + 0.5d0 * xlen
      ycent = prob_lo(2) + 0.5d0 * ylen
      zcent = prob_lo(3) + 0.5d0 * zlen
c
      do k = lo(3)-1, hi(3)+1
       z = ( (dble(k)+0.5)*dx(3) - zcent ) / zlen
       do j = lo(2)-1, hi(2)+1
         y = ( (dble(j)+0.5)*dx(2) - ycent ) / ylen
         do i = lo(1)-1, hi(1)+1
            x = ( (dble(i)+0.5)*dx(1) - xcent ) / xlen
            r = sqrt(x**2 + y**2 + z**2) / 0.1d0
c           coef(i,j,k) = 1.d0 + 1000.d0 * (0.5d0 * (1.0d0 + tanh((2.d0-r))))
            coef(i,j,k) = 1.d0 + 10.d0 * (0.5d0 * (1.0d0 + tanh((2.d0-r))))
c           coef(i,j,k) = 1.d0
         end do
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
      integer i,j,k
c
      if (dir .eq. 0) then 
        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
          edge(i,j,k) = 0.5d0 * (coef(i,j,k) + coef(i-1,j,k))
        end do
        end do
        end do
      else if (dir .eq. 1) then
        do k = lo(3), hi(3)
        do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
          edge(i,j,k) = 0.5d0 * (coef(i,j,k) + coef(i,j-1,k))
        end do
        end do
        end do
      else if (dir .eq. 2) then
        do k = lo(3), hi(3)+1
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           edge(i,j,k) = 0.5d0 * (coef(i,j,k) + coef(i,j,k-1))
        end do
        end do
        end do
      end if
c
      end

