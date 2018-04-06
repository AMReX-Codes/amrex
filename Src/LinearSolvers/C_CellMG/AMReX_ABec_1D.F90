#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>

#include "AMReX_ABec_F.H"
#include "AMReX_ArrayLim.H"
#include "AMReX_CONSTANTS.H"

!-----------------------------------------------------------------------
!      
!     LINESOLVE 
!     Apply the line solve to the state phi for the equation
!     L(phi) = alpha*a(x)*phi(x) - beta*Div(b(x)Grad(phi(x))) = rhs(x)
!     central differenced, according to the arrays of boundary
!     masks (m#) and auxiliary data (f#).
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
           phi,phi_l1,phi_h1, &
           rhs,rhs_l1,rhs_h1, &
           alpha, beta, &
           a,  a_l1,a_h1, &
           bX, bX_l1,bX_h1, &
           f0, f0_l1,f0_h1, &
           m0, m0_l1,m0_h1, &
           f2, f2_l1,f2_h1, &
           m2, m2_l1,m2_h1, &
           lo,hi,nc, &
           h &
           )

      REAL_T alpha, beta
      integer phi_l1,phi_h1
      integer rhs_l1,rhs_h1
      integer a_l1,a_h1
      integer bX_l1,bX_h1
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer nc
      integer f0_l1,f0_h1
      REAL_T f0(DIMV(f0))
      integer f2_l1,f2_h1
      REAL_T f2(DIMV(f2))
      integer m0_l1,m0_h1
      integer m0(DIMV(m0))
      integer m2_l1,m2_h1
      integer m2(DIMV(m2))
      REAL_T  h(BL_SPACEDIM)
      REAL_T   phi(DIMV(phi),nc)
      REAL_T   rhs(DIMV(rhs),nc)
      REAL_T     a(DIMV(a))
      REAL_T    bX(DIMV(bX))

      integer  i, n

      REAL_T dhx, cf0, cf2
      REAL_T delta, gamma, rho, rho_x

      integer LSDIM
      parameter(LSDIM=127)
      REAL_T a_ls(0:LSDIM)
      REAL_T b_ls(0:LSDIM)
      REAL_T c_ls(0:LSDIM)
      REAL_T r_ls(0:LSDIM)
      REAL_T u_ls(0:LSDIM)

      integer ilen
      ilen = hi(1)-lo(1)+1

      dhx = beta/h(1)**2
      do n = 1, nc
             do i = lo(1), hi(1)
     
               cf0 = merge(f0(lo(1)), 0.0D0, &
                    (i .eq. lo(1)) .and. (m0(lo(1)-1).gt.0))
               cf2 = merge(f2(hi(1)), 0.0D0, &
                    (i .eq. hi(1)) .and. (m2(hi(1)+1).gt.0))
     
               delta = dhx*(bX(i)*cf0 + bX(i+1)*cf2)
     
               gamma = alpha*a(i) + dhx*( bX(i) + bX(i+1) )

               a_ls(i-lo(1)) = -dhx*bX(i)
               b_ls(i-lo(1)) = gamma - delta
               c_ls(i-lo(1)) = -dhx*bX(i+1)
               r_ls(i-lo(1)) = rhs(i,n) - phi(i,n)*delta

               if (i .eq. lo(1)) &
                  r_ls(i-lo(1)) = r_ls(i-lo(1)) + dhx*bX(i)*phi(i-1,n)

               if (i .eq. hi(1)) &
                  r_ls(i-lo(1)) = r_ls(i-lo(1)) + dhx*bX(i+1)*phi(i+1,n)
             end do

             call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,ilen)

             do i = lo(1), hi(1)
               phi(i,n) = u_ls(i-lo(1))
             end do
      end do

    end subroutine FORT_LINESOLVE

!-----------------------------------------------------------------------
!
!     Fill in a matrix x vector operator here
!
    subroutine FORT_ADOTX( &
           y,y_l1,y_h1, &
           x,x_l1,x_h1, &
           alpha, beta, &
           a, a_l1,a_h1, &
           bX, bX_l1,bX_h1, &
           lo,hi,nc, &
           h &
           )

      REAL_T alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer y_l1,y_h1
      integer x_l1,x_h1
      integer a_l1,a_h1
      integer bX_l1,bX_h1
      REAL_T  x(DIMV(x),nc)
      REAL_T  y(DIMV(x),nc)
      REAL_T  a(DIMV(a))
      REAL_T bX(DIMV(bX))
      REAL_T h(BL_SPACEDIM)

      integer i,n
      REAL_T dhx

      dhx = beta/h(1)**2

      do n = 1, nc
         do i = lo(1), hi(1)
            y(i,n) = alpha*a(i)*x(i,n) &
                 - dhx* &
                 (   bX(i+1)*( x(i+1,n) - x(i  ,n) ) &
                 -   bX(i  )*( x(i  ,n) - x(i-1,n) ) )
         end do
      end do

    end subroutine FORT_ADOTX

!-----------------------------------------------------------------------
!
!     Fill in a matrix x vector operator here
!
    subroutine FORT_NORMA( &
           res, &
           alpha, beta, &
           a, a_l1,a_h1, &
           bX,bX_l1,bX_h1, &
           lo,hi,nc, &
           h &
           )

      REAL_T res
      REAL_T alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer a_l1,a_h1
      integer bX_l1,bX_h1
      REAL_T  a(DIMV(a))
      REAL_T bX(DIMV(bX))
      REAL_T h(BL_SPACEDIM)

      integer i,n
      REAL_T dhx

      dhx = beta/h(1)**2

      res = 0.0D0
      do n = 1, nc
         do i = lo(1), hi(1)
            res = max(res, &
                 + abs( alpha*a(i) &
                      + dhx * (bX(i+1) + bX(i)) ) &
                 + abs(-dhx*bX(i+1)) + abs(-dhx*bX(i)) )
         end do
      end do

    end subroutine FORT_NORMA

!-----------------------------------------------------------------------
!
!     Fill in fluxes
!
    subroutine FORT_FLUX( &
           x,x_l1,x_h1, &
           alpha, beta, &
           a, a_l1,a_h1, &
           bX,bX_l1,bX_h1, &
           xlo,xhi,nc, &
           h, &
           xflux,xflux_l1,xflux_h1 &
           )

      implicit none

      REAL_T alpha, beta
      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM), nc
      integer x_l1,x_h1
      integer a_l1,a_h1
      integer bX_l1,bX_h1
      integer xflux_l1,xflux_h1
      REAL_T  x(DIMV(x),nc)
      REAL_T  a(DIMV(a))
      REAL_T bX(DIMV(bX))
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T h(BL_SPACEDIM)

      REAL_T dhx
      integer i,n

      dhx = one/h(1)

      do n = 1, nc
         do i = xlo(1), xhi(1)
            xflux(i,n) = - dhx*bX(i)*( x(i,n) - x(i-1,n) )
         end do
      end do

    end subroutine FORT_FLUX

