
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
    subroutine amrex_lp_linesolve ( &
           phi, phi_l1,phi_h1, &
           rhs, rhs_l1,rhs_h1, &
           f0, f0_l1,f0_h1, m0, m0_l1,m0_h1, &
           f2, f2_l1,f2_h1, m2, m2_l1,m2_h1, &
           lo, hi, nc, &
           h) bind(c,name='amrex_lp_linesolve')

      integer nc
      integer phi_l1,phi_h1
      real(amrex_real) phi(phi_l1:phi_h1,nc)
      integer rhs_l1,rhs_h1
      real(amrex_real) rhs(rhs_l1:rhs_h1,nc)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer f0_l1,f0_h1
      integer f2_l1,f2_h1
      real(amrex_real) f0(f0_l1:f0_h1)
      real(amrex_real) f2(f2_l1:f2_h1)
      integer m0_l1,m0_h1
      integer m2_l1,m2_h1
      integer m0(m0_l1:m0_h1)
      integer m2(m2_l1:m2_h1)
      real(amrex_real)  h

      integer  i, n

      real(amrex_real) cf0, cf2
      real(amrex_real) delta, gamma, rho

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

    end subroutine amrex_lp_linesolve
!-----------------------------------------------------------------------
!>
!>     Fill in a matrix x vector operator here
!>
    subroutine amrex_lp_adotx( &
           y, y_l1,y_h1, &
           x, x_l1,x_h1, &
           lo, hi, nc, &
           h &
           ) bind(c,name='amrex_lp_adotx')

      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer y_l1,y_h1
      real(amrex_real) y(y_l1:y_h1,nc)
      integer x_l1,x_h1
      real(amrex_real) x(x_l1:x_h1,nc)
      real(amrex_real) h

      integer i, n
      real(amrex_real) scal

      scal = 1.0D0/h**2

      do n = 1, nc
         do i = lo(1), hi(1)
            y(i,n) = scal* &
             ( x(i-1,n) + x(i+1,n) - 2.d0*x(i,n) )
         end do
      end do

    end subroutine amrex_lp_adotx

!-----------------------------------------------------------------------
!>
!>     Fill in fluxes
!>
    subroutine amrex_lp_flux( &
           x,x_l1,x_h1, &
           xlo,xhi,nc, &
           h, &
           xflux,xflux_l1,xflux_h1 &
           ) bind(c,name='amrex_lp_flux')

      implicit none

      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM), nc
      integer x_l1,x_h1
      integer xflux_l1,xflux_h1
      real(amrex_real)  x(x_l1:x_h1,nc)
      real(amrex_real) xflux(xflux_l1:xflux_h1,nc)
      real(amrex_real) h(BL_SPACEDIM)

      real(amrex_real) dhx
      integer i,n

      dhx = one/h(1)

      do n = 1, nc
         do i = xlo(1), xhi(1)
            xflux(i,n) = - dhx*( x(i,n) - x(i-1,n) )
         end do
      end do

    end subroutine amrex_lp_flux

end module amrex_lp_module
