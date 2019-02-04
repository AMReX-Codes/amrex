
module amrex_abec_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

contains

!-----------------------------------------------------------------------
!>
!>     LINESOLVE
!>     Apply the line solve to the state phi for the equation
!>     ``L(phi) = alpha*a(x)*phi(x) - beta*Div(b(x)Grad(phi(x))) = rhs(x)``
!>     central differenced, according to the arrays of boundary
!>     masks (m#) and auxiliary data (f#).
!>
!>     In general, if the linear operator`` L=gamma*y-rho``, the GS relaxation
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
    subroutine amrex_abec_linesolve ( &
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
           ) bind(c,name='amrex_abec_linesolve')

      use amrex_abec_util_module, only : tridiag

      real(amrex_real) alpha, beta
      integer phi_l1,phi_h1
      integer rhs_l1,rhs_h1
      integer a_l1,a_h1
      integer bX_l1,bX_h1
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer nc
      integer f0_l1,f0_h1
      real(amrex_real) f0(f0_l1:f0_h1)
      integer f2_l1,f2_h1
      real(amrex_real) f2(f2_l1:f2_h1)
      integer m0_l1,m0_h1
      integer m0(m0_l1:m0_h1)
      integer m2_l1,m2_h1
      integer m2(m2_l1:m2_h1)
      real(amrex_real)  h(BL_SPACEDIM)
      real(amrex_real)   phi(phi_l1:phi_h1,nc)
      real(amrex_real)   rhs(rhs_l1:rhs_h1,nc)
      real(amrex_real)     a(a_l1:a_h1)
      real(amrex_real)    bX(bX_l1:bX_h1)

      integer  i, n

      real(amrex_real) dhx, cf0, cf2
      real(amrex_real) delta, gamma, rho, rho_x

      integer LSDIM
      parameter(LSDIM=127)
      real(amrex_real) a_ls(0:LSDIM)
      real(amrex_real) b_ls(0:LSDIM)
      real(amrex_real) c_ls(0:LSDIM)
      real(amrex_real) r_ls(0:LSDIM)
      real(amrex_real) u_ls(0:LSDIM)

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

    end subroutine amrex_abec_linesolve

!-----------------------------------------------------------------------
!>
!>     Fill in a matrix x vector operator here
!>
    subroutine amrex_abec_adotx( &
           y,y_l1,y_h1, &
           x,x_l1,x_h1, &
           alpha, beta, &
           a, a_l1,a_h1, &
           bX, bX_l1,bX_h1, &
           lo,hi,nc, &
           h &
           ) bind(c,name='amrex_abec_adotx')

      real(amrex_real) alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer y_l1,y_h1
      integer x_l1,x_h1
      integer a_l1,a_h1
      integer bX_l1,bX_h1
      real(amrex_real)  x(x_l1:x_h1,nc)
      real(amrex_real)  y(x_l1:x_h1,nc)
      real(amrex_real)  a(a_l1:a_h1)
      real(amrex_real) bX(bX_l1:bX_h1)
      real(amrex_real) h(BL_SPACEDIM)

      integer i,n
      real(amrex_real) dhx

      dhx = beta/h(1)**2

      do n = 1, nc
         do i = lo(1), hi(1)
            y(i,n) = alpha*a(i)*x(i,n) &
                 - dhx* &
                 (   bX(i+1)*( x(i+1,n) - x(i  ,n) ) &
                 -   bX(i  )*( x(i  ,n) - x(i-1,n) ) )
         end do
      end do

    end subroutine amrex_abec_adotx

!-----------------------------------------------------------------------
!>
!>     Fill in a matrix x vector operator here
!>
    subroutine amrex_abec_norma( &
           res, &
           alpha, beta, &
           a, a_l1,a_h1, &
           bX,bX_l1,bX_h1, &
           lo,hi,nc, &
           h &
           ) bind(c,name='amrex_abec_norma')

      real(amrex_real) res
      real(amrex_real) alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer a_l1,a_h1
      integer bX_l1,bX_h1
      real(amrex_real)  a(a_l1:a_h1)
      real(amrex_real) bX(bX_l1:bX_h1)
      real(amrex_real) h(BL_SPACEDIM)

      integer i,n
      real(amrex_real) dhx

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

    end subroutine amrex_abec_norma

!-----------------------------------------------------------------------
!>
!>     Fill in fluxes
!>
    subroutine amrex_abec_flux( &
           x,x_l1,x_h1, &
           alpha, beta, &
           a, a_l1,a_h1, &
           bX,bX_l1,bX_h1, &
           xlo,xhi,nc, &
           h, &
           xflux,xflux_l1,xflux_h1 &
           ) bind(c,name='amrex_abec_flux')

      implicit none

      real(amrex_real) alpha, beta
      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM), nc
      integer x_l1,x_h1
      integer a_l1,a_h1
      integer bX_l1,bX_h1
      integer xflux_l1,xflux_h1
      real(amrex_real)  x(x_l1:x_h1,nc)
      real(amrex_real)  a(a_l1:a_h1)
      real(amrex_real) bX(bX_l1:bX_h1)
      real(amrex_real) xflux(xflux_l1:xflux_h1,nc)
      real(amrex_real) h(BL_SPACEDIM)

      real(amrex_real) dhx
      integer i,n

      dhx = one/h(1)

      do n = 1, nc
         do i = xlo(1), xhi(1)
            xflux(i,n) = - dhx*bX(i)*( x(i,n) - x(i-1,n) )
         end do
      end do

    end subroutine amrex_abec_flux

end module amrex_abec_module
