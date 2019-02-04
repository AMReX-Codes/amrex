module amrex_lp_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

contains

!-----------------------------------------------------------------------
!>
!>     Gauss-Seidel Red-Black (GSRB):
!>     Apply the GSRB relaxation to the state phi for the equation
!>     ``L(phi) = Div(Grad(phi(x))) = rhs(x)`` central differenced, according
!>     to the arrays of boundary masks (m#) and auxiliary data (f#).
!>
!>     In general, if the linear operator ``L=gamma*y-rho``, the GS relaxation
!>     is ``y = (R - rho)/gamma``.  Near a boundary, the ghost data is filled
!>     using a polynomial interpolant based on the "old" phi values, so
!>     ``L=(gamma-delta)*y - rho + delta*yOld``.   The resulting iteration is
!>
!>     ``y = (R - delta*yOld + rho)/(gamma - delta)``
!>
!>     This expression is valid additionally in the interior provided
!>     delta->0 there.  delta is constructed by summing all the
!>     contributions to the central stencil element coming from boundary
!>     interpolants.  The f#s contain the corresponding coefficient of
!>     the interpolating polynomial.  The masks are set > 0 if the boundary
!>     value was filled with an interpolant involving the central stencil
!>     element.
!>
!-----------------------------------------------------------------------
    subroutine amrex_lp_gsrb ( &
           phi, phi_l1,phi_l2,phi_h1,phi_h2, &
           rhs, rhs_l1,rhs_l2,rhs_h1,rhs_h2, &
           f0, f0_l1,f0_l2,f0_h1,f0_h2, m0, m0_l1,m0_l2,m0_h1,m0_h2, &
           f1, f1_l1,f1_l2,f1_h1,f1_h2, m1, m1_l1,m1_l2,m1_h1,m1_h2, &
           f2, f2_l1,f2_l2,f2_h1,f2_h2, m2, m2_l1,m2_l2,m2_h1,m2_h2, &
           f3, f3_l1,f3_l2,f3_h1,f3_h2, m3, m3_l1,m3_l2,m3_h1,m3_h2, &
           lo, hi, blo, bhi, &
           nc, h, redblack &
           ) bind(c,name='amrex_lp_gsrb')

      implicit none

      integer nc
      integer phi_l1,phi_l2,phi_h1,phi_h2
      real(amrex_real) phi(phi_l1:phi_h1,phi_l2:phi_h2,nc)
      integer rhs_l1,rhs_l2,rhs_h1,rhs_h2
      real(amrex_real) rhs(rhs_l1:rhs_h1,rhs_l2:rhs_h2,nc)
      integer  lo(BL_SPACEDIM),  hi(BL_SPACEDIM)
      integer blo(BL_SPACEDIM), bhi(BL_SPACEDIM)
      integer f0_l1,f0_l2,f0_h1,f0_h2
      integer f1_l1,f1_l2,f1_h1,f1_h2
      integer f2_l1,f2_l2,f2_h1,f2_h2
      integer f3_l1,f3_l2,f3_h1,f3_h2
      real(amrex_real) f0(f0_l1:f0_h1,f0_l2:f0_h2)
      real(amrex_real) f1(f1_l1:f1_h1,f1_l2:f1_h2)
      real(amrex_real) f2(f2_l1:f2_h1,f2_l2:f2_h2)
      real(amrex_real) f3(f3_l1:f3_h1,f3_l2:f3_h2)
      integer m0_l1,m0_l2,m0_h1,m0_h2
      integer m1_l1,m1_l2,m1_h1,m1_h2
      integer m2_l1,m2_l2,m2_h1,m2_h2
      integer m3_l1,m3_l2,m3_h1,m3_h2
      integer m0(m0_l1:m0_h1,m0_l2:m0_h2)
      integer m1(m1_l1:m1_h1,m1_l2:m1_h2)
      integer m2(m2_l1:m2_h1,m2_l2:m2_h2)
      integer m3(m3_l1:m3_h1,m3_l2:m3_h2)
      integer redblack
      real(amrex_real)  h

      integer  i, j, ioff, n

      real(amrex_real) cf0, cf1, cf2, cf3
      real(amrex_real) delta, gamma, rho

      gamma = 4.0D0
      do n = 1, nc
         do j = lo(2), hi(2)
            ioff = MOD(lo(1) + j +  redblack, 2)
            do i = lo(1) + ioff,hi(1),2

               cf0 = merge(f0(blo(1),j), 0.0D0, &
                    (i .eq. blo(1)) .and. (m0(blo(1)-1,j).gt.0))
               cf1 = merge(f1(i,blo(2)), 0.0D0, &
                    (j .eq. blo(2)) .and. (m1(i,blo(2)-1).gt.0))
               cf2 = merge(f2(bhi(1),j), 0.0D0, &
                    (i .eq. bhi(1)) .and. (m2(bhi(1)+1,j).gt.0))
               cf3 = merge(f3(i,bhi(2)), 0.0D0, &
                    (j .eq. bhi(2)) .and. (m3(i,bhi(2)+1).gt.0))

               delta = cf0 + cf1 + cf2 + cf3

               rho =  phi(i-1,j,n) + phi(i+1,j,n) &
                    + phi(i,j-1,n) + phi(i,j+1,n)

               phi(i,j,n) = (rhs(i,j,n)*h*h - rho + phi(i,j,n)*delta) &
                    /                (delta - gamma)

            end do
         end do
      end do

    end subroutine amrex_lp_gsrb
!-----------------------------------------------------------------------
!>
!>     Fill in a matrix x vector operator here
!>
    subroutine amrex_lp_adotx( &
           y, y_l1,y_l2,y_h1,y_h2, &
           x, x_l1,x_l2,x_h1,x_h2, &
           lo, hi, nc, &
           h &
           ) bind(c,name='amrex_lp_adotx')

      implicit none

      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer y_l1,y_l2,y_h1,y_h2
      real(amrex_real) y(y_l1:y_h1,y_l2:y_h2,nc)
      integer x_l1,x_l2,x_h1,x_h2
      real(amrex_real) x(x_l1:x_h1,x_l2:x_h2,nc)
      real(amrex_real) h

      integer i, j, n
      real(amrex_real) scal

      scal = 1.0D0/h**2

      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               y(i,j,n) = scal* &
                    ( x(i-1,j,n) + x(i+1,j,n) &
                    + x(i,j-1,n) + x(i,j+1,n) &
                    - 4*x(i,j,n) )
            end do
         end do
      end do

    end subroutine amrex_lp_adotx

!-----------------------------------------------------------------------
!>
!>     Fill in fluxes
!>
    subroutine amrex_lp_flux( &
           x,x_l1,x_l2,x_h1,x_h2, &
           xlo,xhi, &
           ylo,yhi, &
           nc, &
           h, &
           xflux,xflux_l1,xflux_l2,xflux_h1,xflux_h2, &
           yflux,yflux_l1,yflux_l2,yflux_h1,yflux_h2 &
           ) bind(c,name='amrex_lp_flux')

      implicit none

      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM)
      integer ylo(BL_SPACEDIM), yhi(BL_SPACEDIM)
      integer nc
      integer x_l1,x_l2,x_h1,x_h2
      integer xflux_l1,xflux_l2,xflux_h1,xflux_h2
      integer yflux_l1,yflux_l2,yflux_h1,yflux_h2
      real(amrex_real)  x(x_l1:x_h1,x_l2:x_h2,nc)
      real(amrex_real) xflux(xflux_l1:xflux_h1,xflux_l2:xflux_h2,nc)
      real(amrex_real) yflux(yflux_l1:yflux_h1,yflux_l2:yflux_h2,nc)
      real(amrex_real) h(BL_SPACEDIM)

      real(amrex_real) dhx, dhy
      integer i,j,n

      dhx = one/h(1)
      dhy = one/h(2)

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

    end subroutine amrex_lp_flux

end module amrex_lp_module
