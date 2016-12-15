#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <CONSTANTS.H>
#include <REAL.H>

#include "LP_F.H"
#include "ArrayLim.H"

c-----------------------------------------------------------------------
c      
c     Gauss-Seidel Red-Black (GSRB):
c     Apply the GSRB relaxation to the state phi for the equation
c     L(phi) = Div(Grad(phi(x))) = rhs(x) central differenced, according
c     to the arrays of boundary masks (m#) and auxiliary data (f#).
c     
c     In general, if the linear operator L=gamma*y-rho, the GS relaxation
c     is y = (R - rho)/gamma.  Near a boundary, the ghost data is filled
c     using a polynomial interpolant based on the "old" phi values, so
c     L=(gamma-delta)*y - rho + delta*yOld.   The resulting iteration is
c     
c     y = (R - delta*yOld + rho)/(gamma - delta)
c     
c     This expression is valid additionally in the interior provided
c     delta->0 there.  delta is constructed by summing all the
c     contributions to the central stencil element coming from boundary 
c     interpolants.  The f#s contain the corresponding coefficient of 
c     the interpolating polynomial.  The masks are set > 0 if the boundary 
c     value was filled with an interpolant involving the central stencil 
c     element.
c     
c-----------------------------------------------------------------------
      subroutine FORT_GSRB (
     $     phi, DIMS(phi),
     $     rhs, DIMS(rhs),
     $     f0, DIMS(f0), m0, DIMS(m0),
     $     f1, DIMS(f1), m1, DIMS(m1),
     $     f2, DIMS(f2), m2, DIMS(m2),
     $     f3, DIMS(f3), m3, DIMS(m3),
     $     lo, hi, blo, bhi,
     $     nc, h, redblack
     $     )
      implicit none
      integer nc
      integer DIMDEC(phi)
      REAL_T phi(DIMV(phi),nc)
      integer DIMDEC(rhs)
      REAL_T rhs(DIMV(rhs),nc)
      integer  lo(BL_SPACEDIM),  hi(BL_SPACEDIM)
      integer blo(BL_SPACEDIM), bhi(BL_SPACEDIM)
      integer DIMDEC(f0)
      integer DIMDEC(f1)
      integer DIMDEC(f2)
      integer DIMDEC(f3)
      REAL_T f0(DIMV(f0))
      REAL_T f1(DIMV(f1))
      REAL_T f2(DIMV(f2))
      REAL_T f3(DIMV(f3))
      integer DIMDEC(m0)
      integer DIMDEC(m1)
      integer DIMDEC(m2)
      integer DIMDEC(m3)
      integer m0(DIMV(m0))
      integer m1(DIMV(m1))
      integer m2(DIMV(m2))
      integer m3(DIMV(m3))
      integer redblack
      REAL_T  h
c
      integer  i, j, ioff, n
c
      REAL_T cf0, cf1, cf2, cf3
      REAL_T delta, gamma, rho
c
      gamma = 4.0D0
      do n = 1, nc
         do j = lo(2), hi(2)
            ioff = MOD(lo(1) + j +  redblack, 2)
            do i = lo(1) + ioff,hi(1),2
c     
               cf0 = merge(f0(blo(1),j), 0.0D0,
     $              (i .eq. blo(1)) .and. (m0(blo(1)-1,j).gt.0))
               cf1 = merge(f1(i,blo(2)), 0.0D0,
     $              (j .eq. blo(2)) .and. (m1(i,blo(2)-1).gt.0))
               cf2 = merge(f2(bhi(1),j), 0.0D0,
     $              (i .eq. bhi(1)) .and. (m2(bhi(1)+1,j).gt.0))
               cf3 = merge(f3(i,bhi(2)), 0.0D0,
     $              (j .eq. bhi(2)) .and. (m3(i,bhi(2)+1).gt.0))
c     
               delta = cf0 + cf1 + cf2 + cf3
c     
               rho =  phi(i-1,j,n) + phi(i+1,j,n)
     $              + phi(i,j-1,n) + phi(i,j+1,n)
c     
               phi(i,j,n) = (rhs(i,j,n)*h*h - rho + phi(i,j,n)*delta)
     $              /                (delta - gamma)
c     
            end do
         end do
      end do
      
      end
c-----------------------------------------------------------------------
c
c     Fill in a matrix x vector operator here
c
      subroutine FORT_ADOTX(
     $     y, DIMS(y),
     $     x, DIMS(x),
     $     lo, hi, nc,
     $     h
     $     )
      implicit none
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(y)
      REAL_T y(DIMV(y),nc)
      integer DIMDEC(x)
      REAL_T x(DIMV(x),nc)
      REAL_T h
c
      integer i, j, n
      REAL_T scal
c
      scal = 1.0D0/h**2
c
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               y(i,j,n) = scal*
     $              ( x(i-1,j,n) + x(i+1,j,n) 
     $              + x(i,j-1,n) + x(i,j+1,n)
     $              - 4*x(i,j,n) )
            end do
         end do
      end do
c
      end

c-----------------------------------------------------------------------
c
c     Fill in fluxes
c
      subroutine FORT_FLUX(
     $     x,DIMS(x),
     $     xlo,xhi,
     $     ylo,yhi,
     $     nc,
     $     h,
     $     xflux,DIMS(xflux),
     $     yflux,DIMS(yflux)
     $     )
      implicit none
      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM)
      integer ylo(BL_SPACEDIM), yhi(BL_SPACEDIM)
      integer nc
      integer DIMDEC(x)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      REAL_T  x(DIMV(x),nc)
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T h(BL_SPACEDIM)
c
      REAL_T dhx, dhy
      integer i,j,n
c
      dhx = one/h(1)
      dhy = one/h(2)
c
      do n = 1, nc
         do j = xlo(2), xhi(2)
            do i = xlo(1), xhi(1)
               xflux(i,j,n) = - dhx*( x(i,j,n) - x(i-1,j,n) )
            end do
         end do
      end do
      do n = 1, nc
         do j = ylo(2), yhi(2)
            do i = ylo(1), yhi(1)
               yflux(i,j,n) = - dhy*( x(i,j,n) - x(i,j-1,n) )
            end do
         end do
      end do
      end

      
