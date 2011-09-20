#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include "ArrayLim.H"
#include "RHS_F.H"
#include "LO_BCTYPES.H"

c-----------------------------------------------------------------------

      subroutine FORT_SET_RHS (
     $     rhs, r_l1, r_l2, r_l3, r_h1, r_h2, r_h3,
     $     lo, hi, dx, ibnd)

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: ibnd
      integer          :: r_l1, r_l2, r_l3, r_h1, r_h2, r_h3
      double precision :: rhs(r_l1:r_h1, r_l2:r_h2, r_l3:r_h3)
      double precision ::  dx(3)
c
      integer i,j,k
      REAL_T  x,y,z
      REAL_T  pi, phi, fpi, tpi

      pi = 4.d0 * datan(1.d0)

      phi = .25d0 * pi
      tpi = 2.0d0 * pi
      fpi = 4.0d0 * pi
 
      do k = lo(3), hi(3)
         z = (dble(k)+0.5d0)*dx(3)
         do j = lo(2), hi(2)
            y = (dble(j)+0.5d0)*dx(2)
            do i = lo(1), hi(1)
               x = (dble(i)+0.5d0)*dx(1)

               if (ibnd .eq. LO_DIRICHLET) then
                  rhs(i,j,k) = sin(tpi*x) * sin(tpi*y) * sin(tpi*z)  
     $                 +       sin(fpi*x) * sin(fpi*y) * sin(fpi*z)
                  rhs(i,j,k) = rhs(i,j,k) * 12.d0*pi**2
               else if (ibnd .eq. LO_NEUMANN) then
                  rhs(i,j,k) = cos(tpi*x) * cos(tpi*y) * cos(tpi*z)  
     $                 +       cos(fpi*x) * cos(fpi*y) * cos(fpi*z)
                  rhs(i,j,k) = rhs(i,j,k) * 12.d0*pi**2                  
               else if (ibnd .eq. 0) then 
c                 periodic
                  rhs(i,j,k) = sin(tpi*x+phi) * sin(tpi*y+phi) * sin(tpi*z+phi)  
     $                 +       sin(fpi*x) * sin(fpi*y) * sin(fpi*z)
                  rhs(i,j,k) = rhs(i,j,k) * 12.d0*pi**2
               endif
            end do
         end do
      end do
 
      end
