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
      double precision exac(lo(1):hi(1))
      double precision alph(lo(1):hi(1))
      double precision rhs (lo(1):hi(1))
      double precision beta(lo(1)-1:hi(1)+1)
c
      integer i
      double precision x, xc
      double precision r, theta, dbdrfac
      double precision pi, fpi, tpi, fac
c
      pi = 4.d0 * datan(1.d0)
      tpi = 2.0d0 * pi
      fpi = 4.0d0 * pi
      fac = 4.d0 * pi**2

      xc = (prob_hi(1) + prob_lo(1))/2.d0

      theta = 0.5d0*log(3.0) / (w + 1.d-50)
      
      do i = lo(1)-1, hi(1)+1
         x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
               
         r = abs(x-xc)
               
         beta(i) = (sigma-1.d0)/2.d0*tanh(theta*(r-0.25d0))
     $        + (sigma+1.d0)/2.d0
      end do
      
      do i = lo(1), hi(1)
         x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
               
         r = abs(x-xc)
               
         dbdrfac = (sigma-1.d0)/2.d0/(cosh(theta*(r-0.25d0)))**2 * theta/r
         dbdrfac = dbdrfac * b
               
         alph(i) = 1.d0

         if (ibnd .eq. 0 .or. ibnd.eq. LO_NEUMANN) then

            exac(i) = 1.d0 * cos(tpi*x) + .25d0 * cos(fpi*x) 

            rhs(i) = beta(i)*b*fac*(cos(tpi*x) + cos(fpi*x))
     $           + dbdrfac*((x-0.5d0)*(tpi*sin(tpi*x) + pi*sin(fpi*x)))
     $           + a * (cos(tpi*x)+ 0.25d0 * cos(fpi*x))

         else if (ibnd .eq. LO_DIRICHLET) then

            exac(i) = 1.d0 * sin(tpi*x) + .25d0 * sin(fpi*x) 

            rhs(i) = beta(i)*b*fac*(sin(tpi*x) + sin(fpi*x))
     $           + dbdrfac*((x-0.5d0)*(-tpi*cos(tpi*x) - pi*cos(fpi*x)))
     $           + a * (sin(tpi*x) + 0.25d0 * sin(fpi*x))

         else
            print *, 'FORT_SET_COEF: unknow boundary type'
            stop
         endif
      end do
c
      end

