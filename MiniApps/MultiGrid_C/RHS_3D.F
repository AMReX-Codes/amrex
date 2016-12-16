#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include "AMReX_ArrayLim.H"
#include "RHS_F.H"
#include "AMReX_LO_BCTYPES.H"

c-----------------------------------------------------------------------

      subroutine FORT_SET_RHS (
     $     rhs, r_l1, r_l2, r_l3, r_h1, r_h2, r_h3,
     $     lo, hi, dx, a, b, sigma, w, ibnd)

      implicit none

      integer          :: lo(3), hi(3), ibnd
      integer          :: r_l1, r_l2, r_l3, r_h1, r_h2, r_h3
      double precision :: rhs(r_l1:r_h1, r_l2:r_h2, r_l3:r_h3)
      double precision ::  dx(3), a, b, sigma, w
c
      integer i,j,k
      REAL_T  x,y,z,r, theta, beta, dbdrfac
      REAL_T  pi, fpi, tpi, fac

      pi = 4.d0 * datan(1.d0)
      tpi = 2.0d0 * pi
      fpi = 4.0d0 * pi
      fac = 12.d0 * pi**2
 
      theta = 0.5d0*log(3.0) / (w + 1.d-50)

      do k = lo(3), hi(3)
         z = (dble(k)+0.5d0)*dx(3)
         do j = lo(2), hi(2)
            y = (dble(j)+0.5d0)*dx(2)
            do i = lo(1), hi(1)
               x = (dble(i)+0.5d0)*dx(1)

               r = sqrt((x-0.5d0)**2 + (y-0.5d0)**2 + (z-0.5d0)**2)

               beta = (sigma-1.d0)/2.d0*tanh(theta*(r-0.25d0)) + (sigma+1.d0)/2.d0
               beta = beta * b
               dbdrfac = (sigma-1.d0)/2.d0/(cosh(theta*(r-0.25d0)))**2 * theta/r
               dbdrfac = dbdrfac * b

               if (ibnd .eq. 0 .or. ibnd.eq. LO_NEUMANN) then
                  rhs(i,j,k) = beta*fac*(cos(tpi*x) * cos(tpi*y) * cos(tpi*z)  
     $                                 + cos(fpi*x) * cos(fpi*y) * cos(fpi*z))
     $                 + dbdrfac*((x-0.5d0)*(tpi*sin(tpi*x) * cos(tpi*y) * cos(tpi*z)  
     $                                      + pi*sin(fpi*x) * cos(fpi*y) * cos(fpi*z))
     $                          + (y-0.5d0)*(tpi*cos(tpi*x) * sin(tpi*y) * cos(tpi*z)  
     $                                      + pi*cos(fpi*x) * sin(fpi*y) * cos(fpi*z))
     $                          + (z-0.5d0)*(tpi*cos(tpi*x) * cos(tpi*y) * sin(tpi*z)  
     $                                      + pi*cos(fpi*x) * cos(fpi*y) * sin(fpi*z)))
     $                 + a * (cos(tpi*x) * cos(tpi*y) * cos(tpi*z)  
     $                          + 0.25d0 * cos(fpi*x) * cos(fpi*y) * cos(fpi*z))
               else if (ibnd .eq. LO_DIRICHLET) then
                  rhs(i,j,k) = beta*fac*(sin(tpi*x) * sin(tpi*y) * sin(tpi*z)  
     $                                 + sin(fpi*x) * sin(fpi*y) * sin(fpi*z))
     $                 + dbdrfac*((x-0.5d0)*(-tpi*cos(tpi*x) * sin(tpi*y) * sin(tpi*z)  
     $                                       - pi*cos(fpi*x) * sin(fpi*y) * sin(fpi*z))
     $                          + (y-0.5d0)*(-tpi*sin(tpi*x) * cos(tpi*y) * sin(tpi*z)  
     $                                       - pi*sin(fpi*x) * cos(fpi*y) * sin(fpi*z))
     $                          + (z-0.5d0)*(-tpi*sin(tpi*x) * sin(tpi*y) * cos(tpi*z)  
     $                                       - pi*sin(fpi*x) * sin(fpi*y) * cos(fpi*z)))
     $                 + a * (sin(tpi*x) * sin(tpi*y) * sin(tpi*z)  
     $                          + 0.25d0 * sin(fpi*x) * sin(fpi*y) * sin(fpi*z))
               else
                  print *, 'FORT_SET_RHS: unknow boundary type'
                  stop
               endif
            end do
         end do
      end do
 
      end
