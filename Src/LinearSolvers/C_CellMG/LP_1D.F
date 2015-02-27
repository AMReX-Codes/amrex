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
      subroutine FORT_LINESOLVE (
     $     phi, DIMS(phi),
     $     rhs, DIMS(rhs),
     $     f0, DIMS(f0), m0, DIMS(m0),
     $     f2, DIMS(f2), m2, DIMS(m2),
     $     lo, hi, nc,
     $     h)
      integer nc
      integer DIMDEC(phi)
      REAL_T phi(DIMV(phi),nc)
      integer DIMDEC(rhs)
      REAL_T rhs(DIMV(rhs),nc)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer DIMDEC(f0)
      integer DIMDEC(f2)
      REAL_T f0(DIMV(f0))
      REAL_T f2(DIMV(f2))
      integer DIMDEC(m0)
      integer DIMDEC(m2)
      integer m0(DIMV(m0))
      integer m2(DIMV(m2))
      REAL_T  h
c
      integer  i, n
c
      REAL_T cf0, cf2
      REAL_T delta, gamma, rho
c
      gamma = 4.0D0
      do n = 1, nc
         do i = lo(1),hi(1)
c     
               cf0 = merge(f0(lo(1)), 0.0D0,
     $              (i .eq. lo(1)) .and. (m0(lo(1)-1).gt.0))
               cf2 = merge(f2(hi(1)), 0.0D0,
     $              (i .eq. hi(1)) .and. (m2(hi(1)+1).gt.0))
c     
               delta = cf0 + cf2
c     
               rho =  phi(i-1,n) + phi(i+1,n)
c     
               phi(i,n) = (rhs(i,n)*h*h - rho + phi(i,n)*delta)
     $              /                (delta - gamma)
c     
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
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(y)
      REAL_T y(DIMV(y),nc)
      integer DIMDEC(x)
      REAL_T x(DIMV(x),nc)
      REAL_T h
c
      integer i, n
      REAL_T scal
c
      scal = 1.0D0/h**2
c
      do n = 1, nc
         do i = lo(1), hi(1)
            y(i,n) = scal*
     $       ( x(i-1,n) + x(i+1,n) - 2.d0*x(i,n) )
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
     $     xlo,xhi,nc,
     $     h,
     $     xflux,DIMS(xflux)
     $     )
      implicit none
      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM), nc
      integer DIMDEC(x)
      integer DIMDEC(xflux)
      REAL_T  x(DIMV(x),nc)
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T h(BL_SPACEDIM)
c
      REAL_T dhx
      integer i,n
c
      dhx = one/h(1)
c
      do n = 1, nc
         do i = xlo(1), xhi(1)
            xflux(i,n) = - dhx*( x(i,n) - x(i-1,n) )
         end do
      end do
      end

      
