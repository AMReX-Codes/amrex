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
    subroutine FORT_GSRB ( &
           phi, DIMS(phi), &
           rhs, DIMS(rhs), &
           f0, DIMS(f0), m0, DIMS(m0), &
           f1, DIMS(f1), m1, DIMS(m1), &
           f2, DIMS(f2), m2, DIMS(m2), &
           f3, DIMS(f3), m3, DIMS(m3), &
           f4, DIMS(f4), m4, DIMS(m4), &
           f5, DIMS(f5), m5, DIMS(m5), &
           lo, hi, blo, bhi, &
           nc, h, redblack &
           )

      implicit none

      integer nc
      integer DIMDEC(phi)
      REAL_T  phi(DIMV(phi),nc)
      integer DIMDEC(rhs)
      REAL_T  rhs(DIMV(rhs),nc)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer blo(BL_SPACEDIM), bhi(BL_SPACEDIM)
      integer redblack
      integer DIMDEC(f0)
      integer DIMDEC(f1)
      integer DIMDEC(f2)
      integer DIMDEC(f3)
      integer DIMDEC(f4)
      integer DIMDEC(f5)
      REAL_T  f0(DIMV(f0))
      REAL_T  f1(DIMV(f1))
      REAL_T  f2(DIMV(f2))
      REAL_T  f3(DIMV(f3))
      REAL_T  f4(DIMV(f4))
      REAL_T  f5(DIMV(f5))
      integer DIMDEC(m0)
      integer DIMDEC(m1)
      integer DIMDEC(m2)
      integer DIMDEC(m3)
      integer DIMDEC(m4)
      integer DIMDEC(m5)
      integer m0(DIMV(m0))
      integer m1(DIMV(m1))
      integer m2(DIMV(m2))
      integer m3(DIMV(m3))
      integer m4(DIMV(m4))
      integer m5(DIMV(m5))
      REAL_T  h

      integer  i, j, k, ioff, n

      REAL_T cf0, cf1, cf2, cf3, cf4, cf5
      REAL_T delta, gamma, rho

      parameter(gamma = 6.0D0)

      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               ioff = MOD(lo(1) + j + k + redblack,2)
               do i = lo(1) + ioff,hi(1),2

                  cf0 = merge(f0(blo(1),j,k), 0.0D0, &
                       (i .eq. blo(1)) .and. (m0(blo(1)-1,j,k).gt.0))
                  cf1 = merge(f1(i,blo(2),k), 0.0D0, &
                       (j .eq. blo(2)) .and. (m1(i,blo(2)-1,k).gt.0))
                  cf2 = merge(f2(i,j,blo(3)), 0.0D0, &
                       (k .eq. blo(3)) .and. (m2(i,j,blo(3)-1).gt.0))
                  cf3 = merge(f3(bhi(1),j,k), 0.0D0, &
                       (i .eq. bhi(1)) .and. (m3(bhi(1)+1,j,k).gt.0))
                  cf4 = merge(f4(i,bhi(2),k), 0.0D0, &
                       (j .eq. bhi(2)) .and. (m4(i,bhi(2)+1,k).gt.0))
                  cf5 = merge(f5(i,j,bhi(3)), 0.0D0, &
                       (k .eq. bhi(3)) .and. (m5(i,j,bhi(3)+1).gt.0))

                  delta = cf0 + cf1 + cf2 + cf3 + cf4 + cf5

                  rho =  phi(i-1,j,k,n) + phi(i+1,j,k,n) &
                       + phi(i,j-1,k,n) + phi(i,j+1,k,n) &
                       + phi(i,j,k-1,n) + phi(i,j,k+1,n)

                  phi(i,j,k,n) &
                       = (rhs(i,j,k,n)*h*h - rho + phi(i,j,k,n)*delta) &
                       /             (delta - gamma)
                  
               end do
            end do
         end do
      end do

    end subroutine FORT_GSRB
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

      implicit none

      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(y)
      REAL_T  y(DIMV(y),nc)
      integer DIMDEC(x)
      REAL_T  x(DIMV(x),nc)
      REAL_T  h

      integer i, j, k, n
      REAL_T scal

      scal = 1.0D0/h**2

      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  y(i,j,k,n) = scal* &
                       ( x(i-1,j,k,n) + x(i+1,j,k,n) &
                       + x(i,j-1,k,n) + x(i,j+1,k,n) &
                       + x(i,j,k-1,n) + x(i,j,k+1,n) &
                       - 6*x(i,j,k,n) )
               end do
            end do
         end do
      end do

    end subroutine FORT_ADOTX

!-----------------------------------------------------------------------
!
!     Fill in fluxes
!
    subroutine FORT_FLUX( &
           x,DIMS(x), &
           xlo,xhi, &
           ylo,yhi, &
           zlo,zhi, &
           nc, &
           h, &
           xflux,DIMS(xflux), &
           yflux,DIMS(yflux), &
           zflux,DIMS(zflux) &
           )

      implicit none

      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM)
      integer ylo(BL_SPACEDIM), yhi(BL_SPACEDIM)
      integer zlo(BL_SPACEDIM), zhi(BL_SPACEDIM)
      integer nc
      integer DIMDEC(x)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      REAL_T  x(DIMV(x),nc)
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T zflux(DIMV(zflux),nc)
      REAL_T h(BL_SPACEDIM)

      REAL_T dhx, dhy, dhz
      integer i,j,k,n

      dhx = one/h(1)
      dhy = one/h(2)
      dhz = one/h(3)

      do n = 1, nc
         do k = xlo(3), xhi(3)
            do j = xlo(2), xhi(2)
               do i = xlo(1), xhi(1)
                  xflux(i,j,k,n) = - dhx*( x(i,j,k,n) - x(i-1,j,k,n) )
               end do
            end do
         end do

         do k = ylo(3), yhi(3)
            do j = ylo(2), yhi(2)
               do i = ylo(1), yhi(1)
                  yflux(i,j,k,n) = - dhy*( x(i,j,k,n) - x(i,j-1,k,n) )
               end do
            end do
         end do

         do k = zlo(3), zhi(3)
            do j = zlo(2), zhi(2)
               do i = zlo(1), zhi(1)
                  zflux(i,j,k,n) = - dhz*( x(i,j,k,n) - x(i,j,k-1,n) )
               end do
            end do
         end do
      end do

    end subroutine FORT_FLUX
