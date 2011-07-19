#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include "ArrayLim.H"
#include "COEF_F.H"

c-----------------------------------------------------------------------
      subroutine FORT_SET_CC_COEF (
     $     coef, DIMS(coef), lo, hi, dx, prob_lo, prob_hi)

      implicit none

      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(coef)
      REAL_T coef(DIMV(coef))
      REAL_T dx(BL_SPACEDIM)
      REAL_T prob_lo(BL_SPACEDIM)
      REAL_T prob_hi(BL_SPACEDIM)
c
      integer i,j,k
c
      do k = lo(3)-1, hi(3)+1
      do j = lo(2)-1, hi(2)+1
      do i = lo(1)-1, hi(1)+1
         coef(i,j,k) = 1.d0
      end do
      end do
      end do
c
      end

c-----------------------------------------------------------------------
      subroutine FORT_COEF_TO_EDGES (dir, edge, DIMS(edge),
     $                               coef, DIMS(coef), lo, hi)

      implicit none

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

c-----------------------------------------------------------------------
      subroutine FORT_SET_RHS (
     $     rhs, DIMS(rhs), lo, hi, dx, prob_lo, prob_hi)

      implicit none

      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(rhs)
      REAL_T rhs(DIMV(rhs))
      REAL_T dx(BL_SPACEDIM)
      REAL_T prob_lo(BL_SPACEDIM)
      REAL_T prob_hi(BL_SPACEDIM)
c
      integer i,j,k
      REAL_T  x,y,z
      REAL_T  pi, phi, fpi, tpi

      pi = 4.d0 * datan(1.d0)

      phi = .25d0 * pi
      tpi = 2.0d0 * pi
      fpi = 4.0d0 * pi
c
      do k = lo(3), hi(3)
         z = (dble(k)+0.5)*dx(3)
         do j = lo(2), hi(2)
            y = (dble(j)+0.5)*dx(2)
            do i = lo(1), hi(1)
               x = (dble(i)+0.5)*dx(1)

               rhs(i,j,k) =         sin(tpi*x)     * sin(tpi*y)       * sin(tpi*z)  
     $                    + .25d0 * sin(fpi*x+phi) * sin(fpi*y + phi) * sin(fpi*z + phi)

            end do
         end do
      end do
c
      end
