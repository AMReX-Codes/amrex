#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_CONSTANTS.H>
#include <AMReX_REAL.H>

#include "AMReX_LP_F.H"
#include "AMReX_ArrayLim.H"

!-----------------------------------------------------------------------
!      
!     Gauss-Seidel Red-Black (GSRB):
!     Apply the GSRB relaxation to the state phi for the equation
!     L(phi) = Div(Grad(phi(x))) = rhs(x) central differenced, according
!     to the arrays of boundary masks (m#) and auxiliary data (f#).
!     
!     In general, if the linear operator L=gamma*y-rho, the GS relaxation
!     is y = (R - rho)/gamma.  Near a boundary, the ghost data is filled
!     using a polynomial interpolant based on the "old" phi values, so
!     L=(gamma-delta)*y - rho + delta*yOld.   The resulting iteration is
!     
!     y = (R - delta*yOld + rho)/(gamma - delta)
!     
!     This expression is valid additionally in the interior provided
!     delta->0 there.  delta is constructed by summing all the
!     contributions to the central stencil element coming from boundary 
!     interpolants.  The f#s contain the corresponding coefficient of 
!     the interpolating polynomial.  The masks are set > 0 if the boundary 
!     value was filled with an interpolant involving the central stencil 
!     element.
!     
!-----------------------------------------------------------------------
    subroutine FORT_LINESOLVE ( &
           phi, DIMS(phi), &
           rhs, DIMS(rhs), &
           f0, DIMS(f0), m0, DIMS(m0), &
           f2, DIMS(f2), m2, DIMS(m2), &
           lo, hi, nc, &
           h)

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

      integer  i, n

      REAL_T cf0, cf2
      REAL_T delta, gamma, rho

      gamma = 4.0D0
      do n = 1, nc
         do i = lo(1),hi(1)
     
               cf0 = merge(f0(lo(1)), 0.0D0, &
                    (i .eq. lo(1)) .and. (m0(lo(1)-1).gt.0))
               cf2 = merge(f2(hi(1)), 0.0D0, &
                    (i .eq. hi(1)) .and. (m2(hi(1)+1).gt.0))

               delta = cf0 + cf2

               rho =  phi(i-1,n) + phi(i+1,n)

               phi(i,n) = (rhs(i,n)*h*h - rho + phi(i,n)*delta) &
                    /                (delta - gamma)

         end do
      end do
      
    end subroutine FORT_LINESOLVE
!-----------------------------------------------------------------------
!
!     Fill in a matrix x vector operator here
!
    subroutine FORT_ADOTX( &
           y, DIMS(y), &
           x, DIMS(x), &
           lo, hi, nc, &
           h &
           )

      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(y)
      REAL_T y(DIMV(y),nc)
      integer DIMDEC(x)
      REAL_T x(DIMV(x),nc)
      REAL_T h

      integer i, n
      REAL_T scal

      scal = 1.0D0/h**2

      do n = 1, nc
         do i = lo(1), hi(1)
            y(i,n) = scal* &
             ( x(i-1,n) + x(i+1,n) - 2.d0*x(i,n) )
         end do
      end do

    end subroutine FORT_ADOTX

!-----------------------------------------------------------------------
!
!     Fill in fluxes
!
    subroutine FORT_FLUX( &
           x,DIMS(x), &
           xlo,xhi,nc, &
           h, &
           xflux,DIMS(xflux) &
           )

      implicit none

      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM), nc
      integer DIMDEC(x)
      integer DIMDEC(xflux)
      REAL_T  x(DIMV(x),nc)
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T h(BL_SPACEDIM)

      REAL_T dhx
      integer i,n

      dhx = one/h(1)

      do n = 1, nc
         do i = xlo(1), xhi(1)
            xflux(i,n) = - dhx*( x(i,n) - x(i-1,n) )
         end do
      end do

    end subroutine FORT_FLUX
