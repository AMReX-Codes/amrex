#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include "AMReX_ArrayLim.H"
#include "COEF_F.H"
#include "AMReX_LO_BCTYPES.H"

c-----------------------------------------------------------------------
      subroutine FORT_COMP_ASOL (
     $     asol, a_l1, a_l2, a_l3, a_h1, a_h2, a_h3,
     $     lo, hi, dx, ibnd)

      implicit none

      integer          :: lo(3), hi(3), ibnd
      integer          :: a_l1, a_l2, a_l3, a_h1, a_h2, a_h3
      double precision :: asol(a_l1:a_h1, a_l2:a_h2, a_l3:a_h3)
      double precision :: dx(3)
c
      integer i,j,k
      REAL_T  x,y,z
      REAL_T  pi, fpi, tpi

      pi = 4.d0 * datan(1.d0)

      tpi = 2.0d0 * pi
      fpi = 4.0d0 * pi
 
      do k = lo(3), hi(3)
         z = (dble(k)+0.5d0)*dx(3)
         do j = lo(2), hi(2)
            y = (dble(j)+0.5d0)*dx(2)
            do i = lo(1), hi(1)
               x = (dble(i)+0.5d0)*dx(1)

               if (ibnd .eq. 0 .or. ibnd.eq. LO_NEUMANN) then
                  asol(i,j,k) = 1.d0 * cos(tpi*x) * cos(tpi*y) * cos(tpi*z)  
     $                       + .25d0 * cos(fpi*x) * cos(fpi*y) * cos(fpi*z)
               else if (ibnd .eq. LO_DIRICHLET) then
                  asol(i,j,k) = 1.d0 * sin(tpi*x) * sin(tpi*y) * sin(tpi*z)  
     $                       + .25d0 * sin(fpi*x) * sin(fpi*y) * sin(fpi*z)
               else
                  print *, 'FORT_COMP_ASOL: unknow boundary type'
                  stop
               endif
            end do
         end do
      end do
 
      end

c-----------------------------------------------------------------------
      subroutine FORT_SET_ALPHA (
     $     alpha, a_l1, a_l2, a_l3, a_h1, a_h2, a_h3, lo, hi, dx)

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: a_l1, a_l2, a_l3, a_h1, a_h2, a_h3
      double precision :: alpha(a_l1:a_h1, a_l2:a_h2, a_l3:a_h3)
      double precision :: dx(3)
 
      integer i,j,k

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         alpha(i,j,k) = 1.d0
      end do
      end do
      end do
c
      end

c-----------------------------------------------------------------------
      subroutine FORT_SET_CC_COEF (coef, c_l1, c_l2, c_l3, c_h1, c_h2, c_h3, 
     $                             lo, hi, dx, sigma, w)

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: c_l1, c_l2, c_l3, c_h1, c_h2, c_h3
      double precision :: coef(c_l1:c_h1, c_l2:c_h2, c_l3:c_h3)
      double precision :: dx(3), sigma, w
 
      integer i,j,k
      double precision theta, x, y, z, r
 
      theta = 0.5d0*log(3.0) / (w + 1.d-50)

      do k = c_l3, c_h3
         z = (dble(k)+0.5d0)*dx(3)
         do j = c_l2, c_h2
            y = (dble(j)+0.5d0)*dx(2)
            do i = c_l1, c_h1
               x = (dble(i)+0.5d0)*dx(1)

               r = sqrt((x-0.5d0)**2 + (y-0.5d0)**2 + (z-0.5d0)**2)

               coef(i,j,k) = (sigma-1.d0)/2.d0*tanh(theta*(r-0.25d0)) 
     $              + (sigma+1.d0)/2.d0
            end do
         end do
      end do
 
      end

c-----------------------------------------------------------------------
      subroutine FORT_COEF_TO_EDGES (dir, edge, e_l1, e_l2, e_l3, e_h1, e_h2, e_h3,
     $                                    coef, c_l1, c_l2, c_l3, c_h1, c_h2, c_h3,
     $                                    lo, hi)

      implicit none

      integer dir, lo(3), hi(3)
      integer e_l1, e_l2, e_l3, e_h1, e_h2, e_h3
      integer c_l1, c_l2, c_l3, c_h1, c_h2, c_h3
      double precision edge(e_l1:e_h1, e_l2:e_h2, e_l3:e_h3)
      double precision coef(c_l1:c_h1, c_l2:c_h2, c_l3:c_h3)
 
      integer i,j,k
 
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
 
      end

