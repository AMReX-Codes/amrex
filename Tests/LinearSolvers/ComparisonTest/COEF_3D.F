#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_LO_BCTYPES.H>

#include "AMReX_ArrayLim.H"
#include "COEF_F.H"

c-----------------------------------------------------------------------
      subroutine FORT_SET_COEF (
     $     exac, alph, beta, rhs, 
     $     lo, hi, prob_lo, prob_hi, 
     $     dx, a, b, sigma, w, ibnd)

      implicit none

      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer ibnd
      double precision a, b, sigma, w
      double precision prob_lo(BL_SPACEDIM)
      double precision prob_hi(BL_SPACEDIM)
      double precision dx(BL_SPACEDIM)
      double precision exac(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
      double precision alph(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
      double precision rhs (lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
      double precision beta 
     $     (lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1)
c
      integer i,j,k
      double precision x, y, z, xc, yc, zc
      double precision r, theta, dbdrfac
      double precision pi, fpi, tpi, fac
c
      pi = 4.d0 * datan(1.d0)
      tpi = 2.0d0 * pi
      fpi = 4.0d0 * pi
      fac = 12.d0 * pi**2

      xc = (prob_hi(1) + prob_lo(1))/2.d0
      yc = (prob_hi(2) + prob_lo(2))/2.d0
      zc = (prob_hi(3) + prob_lo(3))/2.d0

      theta = 0.5d0*log(3.0) / (w + 1.d-50)
      
      do k = lo(3)-1, hi(3)+1
         z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
         do j = lo(2)-1, hi(2)+1
            y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
            do i = lo(1)-1, hi(1)+1
               x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
               
               r = sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
               
               beta(i,j,k) = (sigma-1.d0)/2.d0*tanh(theta*(r-0.25d0))
     $              + (sigma+1.d0)/2.d0
            end do
         end do
      end do
      
      do k = lo(3), hi(3)
         z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
         do j = lo(2), hi(2)
            y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
            do i = lo(1), hi(1)
               x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
               
               r = sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
               
               dbdrfac = (sigma-1.d0)/2.d0/(cosh(theta*(r-0.25d0)))**2 * theta/r
               dbdrfac = dbdrfac * b
               
               alph(i,j,k) = 1.d0

               if (ibnd .eq. 0 .or. ibnd.eq. LO_NEUMANN) then

                  exac(i,j,k) = 1.d0 * cos(tpi*x) * cos(tpi*y) * cos(tpi*z)  
     $                 + .25d0 * cos(fpi*x) * cos(fpi*y) * cos(fpi*z)

                  rhs(i,j,k) = beta(i,j,k)*b*fac*(cos(tpi*x) * cos(tpi*y) * cos(tpi*z)  
     $                                          + cos(fpi*x) * cos(fpi*y) * cos(fpi*z))
     $                 + dbdrfac*((x-xc)*(tpi*sin(tpi*x) * cos(tpi*y) * cos(tpi*z)  
     $                                      + pi*sin(fpi*x) * cos(fpi*y) * cos(fpi*z))
     $                          + (y-yc)*(tpi*cos(tpi*x) * sin(tpi*y) * cos(tpi*z)  
     $                                      + pi*cos(fpi*x) * sin(fpi*y) * cos(fpi*z))
     $                          + (z-zc)*(tpi*cos(tpi*x) * cos(tpi*y) * sin(tpi*z)  
     $                                      + pi*cos(fpi*x) * cos(fpi*y) * sin(fpi*z)))
     $                 + a * (cos(tpi*x) * cos(tpi*y) * cos(tpi*z)  
     $                          + 0.25d0 * cos(fpi*x) * cos(fpi*y) * cos(fpi*z))

               else if (ibnd .eq. LO_DIRICHLET) then

                  exac(i,j,k) = 1.d0 * sin(tpi*x) * sin(tpi*y) * sin(tpi*z)  
     $                 + .25d0 * sin(fpi*x) * sin(fpi*y) * sin(fpi*z)

                  rhs(i,j,k) = beta(i,j,k)*b*fac*(sin(tpi*x) * sin(tpi*y) * sin(tpi*z)  
     $                                          + sin(fpi*x) * sin(fpi*y) * sin(fpi*z))
     $                 + dbdrfac*((x-xc)*(-tpi*cos(tpi*x) * sin(tpi*y) * sin(tpi*z)  
     $                                       - pi*cos(fpi*x) * sin(fpi*y) * sin(fpi*z))
     $                          + (y-yc)*(-tpi*sin(tpi*x) * cos(tpi*y) * sin(tpi*z)  
     $                                       - pi*sin(fpi*x) * cos(fpi*y) * sin(fpi*z))
     $                          + (z-zc)*(-tpi*sin(tpi*x) * sin(tpi*y) * cos(tpi*z)  
     $                                       - pi*sin(fpi*x) * sin(fpi*y) * cos(fpi*z)))
     $                 + a * (sin(tpi*x) * sin(tpi*y) * sin(tpi*z)  
     $                          + 0.25d0 * sin(fpi*x) * sin(fpi*y) * sin(fpi*z))

               else
                  print *, 'FORT_SET_COEF: unknow boundary type'
                  stop
               endif
            end do
         end do
      end do
c
      end

