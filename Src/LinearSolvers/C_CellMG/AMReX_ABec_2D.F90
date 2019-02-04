

module amrex_abec_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

contains

!-----------------------------------------------------------------------
!>
!>     Gauss-Seidel Red-Black (GSRB):
!>     Apply the GSRB relaxation to the state phi for the equation
!>     ``L(phi) = alpha*a(x)*phi(x) - beta*Div(b(x)Grad(phi(x))) = rhs(x)``
!>     central differenced, according to the arrays of boundary
!>     masks (m#) and auxiliary data (f#).
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
    subroutine amrex_abec_gsrb ( &
           phi,phi_l1,phi_l2,phi_h1,phi_h2, &
           rhs,rhs_l1,rhs_l2,rhs_h1,rhs_h2, &
           alpha, beta, &
           a,  a_l1,a_l2,a_h1,a_h2, &
           bX, bX_l1,bX_l2,bX_h1,bX_h2, &
           bY, bY_l1,bY_l2,bY_h1,bY_h2, &
           f0, f0_l1,f0_l2,f0_h1,f0_h2, &
           m0, m0_l1,m0_l2,m0_h1,m0_h2, &
           f1, f1_l1,f1_l2,f1_h1,f1_h2, &
           m1, m1_l1,m1_l2,m1_h1,m1_h2, &
           f2, f2_l1,f2_l2,f2_h1,f2_h2, &
           m2, m2_l1,m2_l2,m2_h1,m2_h2, &
           f3, f3_l1,f3_l2,f3_h1,f3_h2, &
           m3, m3_l1,m3_l2,m3_h1,m3_h2, &
           lo,hi,blo,bhi, &
           nc,h,redblack &
           ) bind(c,name='amrex_abec_gsrb')

      use amrex_abec_util_module, only : tridiag
      implicit none

      real(amrex_real) alpha, beta
      integer phi_l1,phi_l2,phi_h1,phi_h2
      integer rhs_l1,rhs_l2,rhs_h1,rhs_h2
      integer a_l1,a_l2,a_h1,a_h2
      integer bX_l1,bX_l2,bX_h1,bX_h2
      integer bY_l1,bY_l2,bY_h1,bY_h2
      integer  lo(BL_SPACEDIM),  hi(BL_SPACEDIM)
      integer blo(BL_SPACEDIM), bhi(BL_SPACEDIM)
      integer nc
      integer redblack
      integer f0_l1,f0_l2,f0_h1,f0_h2
      real(amrex_real) f0(f0_l1:f0_h1,f0_l2:f0_h2)
      integer f1_l1,f1_l2,f1_h1,f1_h2
      real(amrex_real) f1(f1_l1:f1_h1,f1_l2:f1_h2)
      integer f2_l1,f2_l2,f2_h1,f2_h2
      real(amrex_real) f2(f2_l1:f2_h1,f2_l2:f2_h2)
      integer f3_l1,f3_l2,f3_h1,f3_h2
      real(amrex_real) f3(f3_l1:f3_h1,f3_l2:f3_h2)
      integer m0_l1,m0_l2,m0_h1,m0_h2
      integer m0(m0_l1:m0_h1,m0_l2:m0_h2)
      integer m1_l1,m1_l2,m1_h1,m1_h2
      integer m1(m1_l1:m1_h1,m1_l2:m1_h2)
      integer m2_l1,m2_l2,m2_h1,m2_h2
      integer m2(m2_l1:m2_h1,m2_l2:m2_h2)
      integer m3_l1,m3_l2,m3_h1,m3_h2
      integer m3(m3_l1:m3_h1,m3_l2:m3_h2)
      real(amrex_real)  h(BL_SPACEDIM)
      real(amrex_real)   phi(phi_l1:phi_h1,phi_l2:phi_h2,nc)
      real(amrex_real)   rhs(rhs_l1:rhs_h1,rhs_l2:rhs_h2,nc)
      real(amrex_real)     a(a_l1:a_h1,a_l2:a_h2)
      real(amrex_real)    bX(bX_l1:bX_h1,bX_l2:bX_h2)
      real(amrex_real)    bY(bY_l1:bY_h1,bY_l2:bY_h2)

      integer  i, j, ioff, joff, n

      real(amrex_real) dhx, dhy, cf0, cf1, cf2, cf3
      real(amrex_real) delta, gamma, rho, rho_x, rho_y

      integer LSDIM
      parameter(LSDIM=127)
      real(amrex_real) a_ls(0:LSDIM)
      real(amrex_real) b_ls(0:LSDIM)
      real(amrex_real) c_ls(0:LSDIM)
      real(amrex_real) r_ls(0:LSDIM)
      real(amrex_real) u_ls(0:LSDIM)

      integer do_line
      integer ilen,jlen

      if (h(2) .gt. 1.5D0*h(1)) then
        do_line = 1
        ilen = hi(1)-lo(1)+1
        if (ilen .gt. LSDIM) then
#ifdef AMREX_DEBUG
          print *,'TOO BIG FOR LINE SOLVE IN GSRB: ilen = ',ilen
#endif
          call bl_error("stop")
        end if
      else if (h(1) .gt. 1.5D0*h(2)) then
        do_line = 2
        jlen = hi(2)-lo(2)+1
        if (jlen .gt. LSDIM) then
#ifdef AMREX_DEBUG
          print *,'TOO BIG FOR LINE SOLVE IN GSRB: jlen = ',jlen
#endif
          call bl_error("stop")
        end if
      else
        do_line = 0
      end if

      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      do n = 1, nc
       if (do_line .eq. 0) then
         do j = lo(2), hi(2)
            ioff = MOD(lo(1) + j + redblack, 2)
            do i = lo(1) + ioff,hi(1),2

               cf0 = merge(f0(blo(1),j), 0.0D0, &
                    (i .eq. blo(1)) .and. (m0(blo(1)-1,j).gt.0))
               cf1 = merge(f1(i,blo(2)), 0.0D0, &
                    (j .eq. blo(2)) .and. (m1(i,blo(2)-1).gt.0))
               cf2 = merge(f2(bhi(1),j), 0.0D0, &
                    (i .eq. bhi(1)) .and. (m2(bhi(1)+1,j).gt.0))
               cf3 = merge(f3(i,bhi(2)), 0.0D0, &
                    (j .eq. bhi(2)) .and. (m3(i,bhi(2)+1).gt.0))

               delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2) &
                    +  dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)

               gamma = alpha*a(i,j) &
                    +   dhx*( bX(i,j) + bX(i+1,j) ) &
                    +   dhy*( bY(i,j) + bY(i,j+1) )

               rho = dhx*(bX(i,j)*phi(i-1,j,n) + bX(i+1,j)*phi(i+1,j,n)) &
                    +dhy*(bY(i,j)*phi(i,j-1,n) + bY(i,j+1)*phi(i,j+1,n))

               phi(i,j,n) = (rhs(i,j,n) + rho - phi(i,j,n)*delta) &
                    /                (gamma - delta)

            end do
         end do
       else if (do_line .eq. 2) then
         ioff = MOD(lo(1) + redblack, 2)
         do i = lo(1) + ioff,hi(1),2
             do j = lo(2), hi(2)

               cf0 = merge(f0(blo(1),j), 0.0D0, &
                    (i .eq. blo(1)) .and. (m0(blo(1)-1,j).gt.0))
               cf1 = merge(f1(i,blo(2)), 0.0D0, &
                    (j .eq. blo(2)) .and. (m1(i,blo(2)-1).gt.0))
               cf2 = merge(f2(bhi(1),j), 0.0D0, &
                    (i .eq. bhi(1)) .and. (m2(bhi(1)+1,j).gt.0))
               cf3 = merge(f3(i,bhi(2)), 0.0D0, &
                    (j .eq. bhi(2)) .and. (m3(i,bhi(2)+1).gt.0))

               delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2) &
                     + dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)

               gamma = alpha*a(i,j) &
                    +   dhx*( bX(i,j) + bX(i+1,j) ) &
                    +   dhy*( bY(i,j) + bY(i,j+1) )

               rho_x = dhx*(bX(i,j)*phi(i-1,j,n) + bX(i+1,j)*phi(i+1,j,n))

               a_ls(j-lo(2)) = -dhy*bY(i,j)
               b_ls(j-lo(2)) = gamma - delta
               c_ls(j-lo(2)) = -dhy*bY(i,j+1)
               r_ls(j-lo(2)) = rhs(i,j,n) + rho_x - phi(i,j,n)*delta

               if (j .eq. lo(2)) &
                  r_ls(j-lo(2)) = r_ls(j-lo(2)) + dhy*bY(i,j)*phi(i,j-1,n)

               if (j .eq. hi(2)) &
                  r_ls(j-lo(2)) = r_ls(j-lo(2)) + dhy*bY(i,j+1)*phi(i,j+1,n)

             end do

             call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,jlen)

             do j = lo(2), hi(2)
               phi(i,j,n) = u_ls(j-lo(2))
             end do
         end do

       else if (do_line .eq. 1) then

           joff = MOD(lo(2) + redblack, 2)
           do j = lo(2) + joff,hi(2),2
             do i = lo(1), hi(1)

               cf0 = merge(f0(blo(1),j), 0.0D0, &
                    (i .eq. blo(1)) .and. (m0(blo(1)-1,j).gt.0))
               cf1 = merge(f1(i,blo(2)), 0.0D0, &
                    (j .eq. blo(2)) .and. (m1(i,blo(2)-1).gt.0))
               cf2 = merge(f2(bhi(1),j), 0.0D0, &
                    (i .eq. bhi(1)) .and. (m2(bhi(1)+1,j).gt.0))
               cf3 = merge(f3(i,bhi(2)), 0.0D0, &
                    (j .eq. bhi(2)) .and. (m3(i,bhi(2)+1).gt.0))

               delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2) &
                     + dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)

               gamma = alpha*a(i,j) &
                    +   dhx*( bX(i,j) + bX(i+1,j) ) &
                    +   dhy*( bY(i,j) + bY(i,j+1) )

               rho_y = dhy*(bY(i,j)*phi(i,j-1,n) + bY(i,j+1)*phi(i,j+1,n))

               a_ls(i-lo(1)) = -dhx*bX(i,j)
               b_ls(i-lo(1)) = gamma - delta
               c_ls(i-lo(1)) = -dhx*bX(i+1,j)
               r_ls(i-lo(1)) = rhs(i,j,n) + rho_y - phi(i,j,n)*delta

               if (i .eq. lo(1)) &
                  r_ls(i-lo(1)) = r_ls(i-lo(1)) + dhx*bX(i,j)*phi(i-1,j,n)

               if (i .eq. hi(1)) &
                  r_ls(i-lo(1)) = r_ls(i-lo(1)) + dhx*bX(i+1,j)*phi(i+1,j,n)
             end do

             call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,ilen)

             do i = lo(1), hi(1)
               phi(i,j,n) = u_ls(i-lo(1))
             end do
         end do

       else
#ifdef AMREX_DEBUG
         print *,'BOGUS DO_LINE '
#endif
         call bl_error("stop")
       end if
      end do

    end subroutine amrex_abec_gsrb

!-----------------------------------------------------------------------
!>
!>     JACOBI:
!>     Apply the JACOBI relaxation to the state phi for the equation
!>     ``L(phi) = alpha*a(x)*phi(x) - beta*Div(b(x)Grad(phi(x))) = rhs(x)``
!>     central differenced, according to the arrays of boundary
!>     masks (m#) and auxiliary data (f#).
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
    subroutine amrex_abec_jacobi ( &
           phi,phi_l1,phi_l2,phi_h1,phi_h2, &
           rhs,rhs_l1,rhs_l2,rhs_h1,rhs_h2, &
           alpha, beta, &
           a,  a_l1,a_l2,a_h1,a_h2, &
           bX, bX_l1,bX_l2,bX_h1,bX_h2, &
           bY, bY_l1,bY_l2,bY_h1,bY_h2, &
           f0, f0_l1,f0_l2,f0_h1,f0_h2, &
           m0, m0_l1,m0_l2,m0_h1,m0_h2, &
           f1, f1_l1,f1_l2,f1_h1,f1_h2, &
           m1, m1_l1,m1_l2,m1_h1,m1_h2, &
           f2, f2_l1,f2_l2,f2_h1,f2_h2, &
           m2, m2_l1,m2_l2,m2_h1,m2_h2, &
           f3, f3_l1,f3_l2,f3_h1,f3_h2, &
           m3, m3_l1,m3_l2,m3_h1,m3_h2, &
           lo,hi,nc, &
           h &
           ) bind(c,name='amrex_abec_jacobi')

      implicit none

      real(amrex_real) alpha, beta
      integer phi_l1,phi_l2,phi_h1,phi_h2
      integer rhs_l1,rhs_l2,rhs_h1,rhs_h2
      integer a_l1,a_l2,a_h1,a_h2
      integer bX_l1,bX_l2,bX_h1,bX_h2
      integer bY_l1,bY_l2,bY_h1,bY_h2
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer nc
      integer f0_l1,f0_l2,f0_h1,f0_h2
      real(amrex_real) f0(f0_l1:f0_h1,f0_l2:f0_h2)
      integer f1_l1,f1_l2,f1_h1,f1_h2
      real(amrex_real) f1(f1_l1:f1_h1,f1_l2:f1_h2)
      integer f2_l1,f2_l2,f2_h1,f2_h2
      real(amrex_real) f2(f2_l1:f2_h1,f2_l2:f2_h2)
      integer f3_l1,f3_l2,f3_h1,f3_h2
      real(amrex_real) f3(f3_l1:f3_h1,f3_l2:f3_h2)
      integer m0_l1,m0_l2,m0_h1,m0_h2
      integer m0(m0_l1:m0_h1,m0_l2:m0_h2)
      integer m1_l1,m1_l2,m1_h1,m1_h2
      integer m1(m1_l1:m1_h1,m1_l2:m1_h2)
      integer m2_l1,m2_l2,m2_h1,m2_h2
      integer m2(m2_l1:m2_h1,m2_l2:m2_h2)
      integer m3_l1,m3_l2,m3_h1,m3_h2
      integer m3(m3_l1:m3_h1,m3_l2:m3_h2)
      real(amrex_real)  h(BL_SPACEDIM)
      real(amrex_real)   phi(phi_l1:phi_h1,phi_l2:phi_h2,nc)
      real(amrex_real)   rhs(rhs_l1:rhs_h1,rhs_l2:rhs_h2,nc)
      real(amrex_real)     a(a_l1:a_h1,a_l2:a_h2)
      real(amrex_real)    bX(bX_l1:bX_h1,bX_l2:bX_h2)
      real(amrex_real)    bY(bY_l1:bY_h1,bY_l2:bY_h2)

      integer  i, j, n

      real(amrex_real) dhx, dhy, cf0, cf1, cf2, cf3
      real(amrex_real) delta, gamma, rho

      real(amrex_real), allocatable :: phinew(:,:)

      allocate(phinew(lo(1):hi(1),lo(2):hi(2)))

      dhx = beta/h(1)**2
      dhy = beta/h(2)**2

      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1),hi(1)

               cf0 = merge(f0(lo(1),j), 0.0D0, &
                    (i .eq. lo(1)) .and. (m0(lo(1)-1,j).gt.0))
               cf1 = merge(f1(i,lo(2)), 0.0D0, &
                    (j .eq. lo(2)) .and. (m1(i,lo(2)-1).gt.0))
               cf2 = merge(f2(hi(1),j), 0.0D0, &
                    (i .eq. hi(1)) .and. (m2(hi(1)+1,j).gt.0))
               cf3 = merge(f3(i,hi(2)), 0.0D0, &
                    (j .eq. hi(2)) .and. (m3(i,hi(2)+1).gt.0))

               delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2) &
                    +  dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)

               gamma = alpha*a(i,j) &
                    +   dhx*( bX(i,j) + bX(i+1,j) ) &
                    +   dhy*( bY(i,j) + bY(i,j+1) )

               rho = dhx*(bX(i,j)*phi(i-1,j,n) + bX(i+1,j)*phi(i+1,j,n)) &
                    +dhy*(bY(i,j)*phi(i,j-1,n) + bY(i,j+1)*phi(i,j+1,n))

               phinew(i,j) = (rhs(i,j,n) + rho - phi(i,j,n)*delta) &
                    /                (gamma - delta)

            end do
         end do

         phi(lo(1):hi(1),lo(2):hi(2),n) = phinew(lo(1):hi(1),lo(2):hi(2))

      end do

      deallocate(phinew)

    end subroutine amrex_abec_jacobi

!-----------------------------------------------------------------------
!>
!>     Fill in a matrix x vector operator here
!>
    subroutine amrex_abec_adotx( &
           y,y_l1,y_l2,y_h1,y_h2, &
           x,x_l1,x_l2,x_h1,x_h2, &
           alpha, beta, &
           a, a_l1,a_l2,a_h1,a_h2, &
           bX,bX_l1,bX_l2,bX_h1,bX_h2, &
           bY,bY_l1,bY_l2,bY_h1,bY_h2, &
           lo,hi,nc, &
           h &
           ) bind(c,name='amrex_abec_adotx')

      implicit none

      real(amrex_real) alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer y_l1,y_l2,y_h1,y_h2
      integer x_l1,x_l2,x_h1,x_h2
      integer a_l1,a_l2,a_h1,a_h2
      integer bX_l1,bX_l2,bX_h1,bX_h2
      integer bY_l1,bY_l2,bY_h1,bY_h2
      real(amrex_real)  y(y_l1:y_h1,y_l2:y_h2,nc)
      real(amrex_real)  x(x_l1:x_h1,x_l2:x_h2,nc)
      real(amrex_real)  a(a_l1:a_h1,a_l2:a_h2)
      real(amrex_real) bX(bX_l1:bX_h1,bX_l2:bX_h2)
      real(amrex_real) bY(bY_l1:bY_h1,bY_l2:bY_h2)
      real(amrex_real) h(BL_SPACEDIM)

      integer i,j,n
      real(amrex_real) dhx,dhy

      dhx = beta/h(1)**2
      dhy = beta/h(2)**2

      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               y(i,j,n) = alpha*a(i,j)*x(i,j,n) &
                    - dhx* &
                    (   bX(i+1,j)*( x(i+1,j,n) - x(i  ,j,n) ) &
                    -   bX(i  ,j)*( x(i  ,j,n) - x(i-1,j,n) ) ) &
                    - dhy* &
                    (   bY(i,j+1)*( x(i,j+1,n) - x(i,j  ,n) ) &
                    -   bY(i,j  )*( x(i,j  ,n) - x(i,j-1,n) ) )
            end do
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
           a, a_l1,a_l2,a_h1,a_h2, &
           bX,bX_l1,bX_l2,bX_h1,bX_h2, &
           bY,bY_l1,bY_l2,bY_h1,bY_h2, &
           lo,hi,nc, &
           h &
           ) bind(c,name='amrex_abec_norma')

      implicit none

      real(amrex_real) res
      real(amrex_real) alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer a_l1,a_l2,a_h1,a_h2
      integer bX_l1,bX_l2,bX_h1,bX_h2
      integer bY_l1,bY_l2,bY_h1,bY_h2
      real(amrex_real)  a(a_l1:a_h1,a_l2:a_h2)
      real(amrex_real) bX(bX_l1:bX_h1,bX_l2:bX_h2)
      real(amrex_real) bY(bY_l1:bY_h1,bY_l2:bY_h2)
      real(amrex_real) h(BL_SPACEDIM)

      integer i,j,n
      real(amrex_real) dhx,dhy

      dhx = beta/h(1)**2
      dhy = beta/h(2)**2

      res = 0.0D0
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               res = max(res, &
                    + abs( alpha*a(i,j) &
                         + dhx * (bX(i+1,j) + bX(i,j)) &
                         + dhy * (bY(i,j+1) + bY(i,j))) &
                    + abs(-dhx*bX(i+1,j)) + abs(-dhx*bX(i,j)) &
                    + abs(-dhy*bY(i,j+1)) + abs(-dhy*bY(i,j)))
            end do
         end do
      end do

    end subroutine amrex_abec_norma

!-----------------------------------------------------------------------
!>
!>     Fill in fluxes
!>
    subroutine amrex_abec_flux( &
           x,x_l1,x_l2,x_h1,x_h2, &
           alpha, beta, &
           a, a_l1,a_l2,a_h1,a_h2, &
           bX,bX_l1,bX_l2,bX_h1,bX_h2, &
           bY,bY_l1,bY_l2,bY_h1,bY_h2, &
           xlo,xhi, &
           ylo,yhi, &
           nc, &
           h, &
           xflux,xflux_l1,xflux_l2,xflux_h1,xflux_h2, &
           yflux,yflux_l1,yflux_l2,yflux_h1,yflux_h2 &
           ) bind(c,name='amrex_abec_flux')

      implicit none

      real(amrex_real) alpha, beta
      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM)
      integer ylo(BL_SPACEDIM), yhi(BL_SPACEDIM)
      integer nc
      integer x_l1,x_l2,x_h1,x_h2
      integer a_l1,a_l2,a_h1,a_h2
      integer bX_l1,bX_l2,bX_h1,bX_h2
      integer bY_l1,bY_l2,bY_h1,bY_h2
      integer xflux_l1,xflux_l2,xflux_h1,xflux_h2
      integer yflux_l1,yflux_l2,yflux_h1,yflux_h2
      real(amrex_real)  x(x_l1:x_h1,x_l2:x_h2,nc)
      real(amrex_real)  a(a_l1:a_h1,a_l2:a_h2)
      real(amrex_real) bX(bX_l1:bX_h1,bX_l2:bX_h2)
      real(amrex_real) bY(bY_l1:bY_h1,bY_l2:bY_h2)
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
               xflux(i,j,n) = - dhx*bX(i,j)*( x(i,j,n) - x(i-1,j,n) )
            end do
         end do
      end do
      do n = 1, nc
         do j = ylo(2), yhi(2)
            do i = ylo(1), yhi(1)
               yflux(i,j,n) = - dhy*bY(i,j)*( x(i,j,n) - x(i,j-1,n) )
            end do
         end do
      end do

    end subroutine amrex_abec_flux

end module amrex_abec_module
