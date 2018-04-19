
module amrex_abec2_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

  include 'AMReX_lo_bctypes.fi'

contains

      subroutine amrex_ab2_gsrb (lo, hi, blo, bhi, &
     &     phi, phi_l1, phi_l2, phi_h1, phi_h2,&
     &     res, res_l1, res_l2, res_h1, res_h2,&
     &     alpha, beta,                        &
     &     a,   a_l1,   a_l2,   a_h1,   a_h2,  &
     &     bX,  bX_l1,  bX_l2,  bX_h1,  bX_h2, &
     &     bY,  bY_l1,  bY_l2,  bY_h1,  bY_h2, &
     &     f0,  f0_l1,  f0_l2,  f0_h1,  f0_h2, &
     &     m0,  m0_l1,  m0_l2,  m0_h1,  m0_h2, &
     &     f1,  f1_l1,  f1_l2,  f1_h1,  f1_h2, &
     &     m1,  m1_l1,  m1_l2,  m1_h1,  m1_h2, &
     &     f2,  f2_l1,  f2_l2,  f2_h1,  f2_h2, &
     &     m2,  m2_l1,  m2_l2,  m2_h1,  m2_h2, &
     &     f3,  f3_l1,  f3_l2,  f3_h1,  f3_h2, &
     &     m3,  m3_l1,  m3_l2,  m3_h1,  m3_h2, &
     &     nc, h, redBlack) bind(c,name='amrex_ab2_gsrb')

      use amrex_abec_util_module, only : tridiag

      implicit none
      integer lo(2), hi(2), blo(2), bhi(2), nc, redblack
      integer phi_l1, phi_l2, phi_h1, phi_h2
      integer res_l1, res_l2, res_h1, res_h2      
      integer a_l1,   a_h1,   a_l2,    a_h2
      integer bX_l1,  bX_h1,  bX_l2,   bX_h2
      integer bY_l1,  bY_h1,  bY_l2,   bY_h2
      integer f0_l1,  f0_h1,  f0_l2,   f0_h2
      integer f1_l1,  f1_h1,  f1_l2,   f1_h2
      integer f2_l1,  f2_h1,  f2_l2,   f2_h2
      integer f3_l1,  f3_h1,  f3_l2,   f3_h2
      integer m0_l1,  m0_h1,  m0_l2,   m0_h2
      integer m1_l1,  m1_h1,  m1_l2,   m1_h2
      integer m2_l1,  m2_h1,  m2_l2,   m2_h2
      integer m3_l1,  m3_h1,  m3_l2,   m3_h2

      real(amrex_real)  phi(phi_l1:phi_h1, phi_l2:phi_h2,1:nc)
      real(amrex_real)  res(res_l1:res_h1, res_l2:res_h2,1:nc)
      real(amrex_real)  alpha, beta

      real(amrex_real)   a( a_l1:a_h1,  a_l2:a_h2)
      real(amrex_real)  bX(bX_l1:bX_h1, bX_l2:bX_h2)
      real(amrex_real)  bY(bY_l1:bY_h1, bY_l2:bY_h2)
      real(amrex_real)  f0(f0_l1:f0_h1, f0_l2:f0_h2)
      real(amrex_real)  f1(f1_l1:f1_h1, f1_l2:f1_h2)
      real(amrex_real)  f2(f2_l1:f2_h1, f2_l2:f2_h2)
      real(amrex_real)  f3(f3_l1:f3_h1, f3_l2:f3_h2)
      integer  m0(m0_l1:m0_h1, m0_l2:m0_h2)
      integer  m1(m1_l1:m1_h1, m1_l2:m1_h2)
      integer  m2(m2_l1:m2_h1, m2_l2:m2_h2)
      integer  m3(m3_l1:m3_h1, m3_l2:m3_h2)
      real(amrex_real)  h(2)

      integer  i, j, ioff, n

      real(amrex_real) dhx, dhy, cf0, cf1, cf2, cf3
      real(amrex_real) delta, gamma, rho_x, rho_y

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
          print *,'TOO BIG FOR LINE SOLVE IN GSRB: ilen = ',ilen
          call bl_error("stop")
        end if
      else if (h(1) .gt. 1.5D0*h(2)) then
        do_line = 2
        jlen = hi(2)-lo(2)+1
        if (jlen .gt. LSDIM) then
          print *,'TOO BIG FOR LINE SOLVE IN GSRB: jlen = ',jlen
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
            ioff = MOD(j + redblack,2)
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
     
     
               phi(i,j,n) = phi(i,j,n) + res(i,j,n)/(gamma - delta)
     
            end do
         end do
       else if (do_line .eq. 2) then

          print *,'Need to fix line solver in ABec2_2D.F'
          call bl_abort()

         do i = lo(1) + redblack,hi(1),2
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
               r_ls(j-lo(2)) = res(i,j,n) + rho_x - phi(i,j,n)*delta

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

          print *,'Need to fix line solver in ABec2_2D.F'
          call bl_abort()

           do j = lo(2) + redblack,hi(2),2
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
               r_ls(i-lo(1)) = res(i,j,n) + rho_y - phi(i,j,n)*delta

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
         print *,'BOGUS DO_LINE '
         call bl_error("stop")
       end if
      end do

      end


      subroutine amrex_ab2_jacobi (lo, hi, &
     &     phi, phi_l1, phi_l2, phi_h1, phi_h2,&
     &     res, res_l1, res_l2, res_h1, res_h2,&
     &     alpha, beta,                        &
     &     a,   a_l1,   a_l2,   a_h1,   a_h2,  &
     &     bX,  bX_l1,  bX_l2,  bX_h1,  bX_h2, &
     &     bY,  bY_l1,  bY_l2,  bY_h1,  bY_h2, &
     &     f0,  f0_l1,  f0_l2,  f0_h1,  f0_h2, &
     &     m0,  m0_l1,  m0_l2,  m0_h1,  m0_h2, &
     &     f1,  f1_l1,  f1_l2,  f1_h1,  f1_h2, &
     &     m1,  m1_l1,  m1_l2,  m1_h1,  m1_h2, &
     &     f2,  f2_l1,  f2_l2,  f2_h1,  f2_h2, &
     &     m2,  m2_l1,  m2_l2,  m2_h1,  m2_h2, &
     &     f3,  f3_l1,  f3_l2,  f3_h1,  f3_h2, &
     &     m3,  m3_l1,  m3_l2,  m3_h1,  m3_h2, &
     &     nc, h) bind(c,name='amrex_ab2_jacobi')

      implicit none
      integer lo(2), hi(2), nc
      integer phi_l1, phi_l2, phi_h1, phi_h2
      integer res_l1, res_l2, res_h1, res_h2      
      integer a_l1,   a_h1,   a_l2,    a_h2
      integer bX_l1,  bX_h1,  bX_l2,   bX_h2
      integer bY_l1,  bY_h1,  bY_l2,   bY_h2
      integer f0_l1,  f0_h1,  f0_l2,   f0_h2
      integer f1_l1,  f1_h1,  f1_l2,   f1_h2
      integer f2_l1,  f2_h1,  f2_l2,   f2_h2
      integer f3_l1,  f3_h1,  f3_l2,   f3_h2
      integer m0_l1,  m0_h1,  m0_l2,   m0_h2
      integer m1_l1,  m1_h1,  m1_l2,   m1_h2
      integer m2_l1,  m2_h1,  m2_l2,   m2_h2
      integer m3_l1,  m3_h1,  m3_l2,   m3_h2

      real(amrex_real)  phi(phi_l1:phi_h1, phi_l2:phi_h2,1:nc)
      real(amrex_real)  res(res_l1:res_h1, res_l2:res_h2,1:nc)
      real(amrex_real)  alpha, beta

      real(amrex_real)   a( a_l1:a_h1,  a_l2:a_h2)
      real(amrex_real)  bX(bX_l1:bX_h1, bX_l2:bX_h2)
      real(amrex_real)  bY(bY_l1:bY_h1, bY_l2:bY_h2)
      real(amrex_real)  f0(f0_l1:f0_h1, f0_l2:f0_h2)
      real(amrex_real)  f1(f1_l1:f1_h1, f1_l2:f1_h2)
      real(amrex_real)  f2(f2_l1:f2_h1, f2_l2:f2_h2)
      real(amrex_real)  f3(f3_l1:f3_h1, f3_l2:f3_h2)
      integer  m0(m0_l1:m0_h1, m0_l2:m0_h2)
      integer  m1(m1_l1:m1_h1, m1_l2:m1_h2)
      integer  m2(m2_l1:m2_h1, m2_l2:m2_h2)
      integer  m3(m3_l1:m3_h1, m3_l2:m3_h2)
      real(amrex_real)  h(2)

      integer  i, j, n

      real(amrex_real) dhx, dhy, cf0, cf1, cf2, cf3
      real(amrex_real) delta, gamma
      real(amrex_real), allocatable :: phit(:,:)
      allocate(phit(lo(1):hi(1),lo(2):hi(2)))

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
     
               phit = phi(i,j,n) + res(i,j,n) / (gamma - delta)
     
            end do
         end do
         phi(lo(1):hi(1),lo(2):hi(2),n) = phit         
      end do
      deallocate(phit)
      end


      subroutine amrex_ab2_bndryrlx(lo, hi, &
     &     den,  f_l1, f_l2, f_h1, f_h2, &
     &     mask, m_l1, m_l2, m_h1, m_h2, &
     &     cdir, bct, bcl, maxorder, h) bind(c,name='amrex_ab2_bndryrlx')

      use amrex_lo_util_module, only : polyInterpCoeff
      implicit none
      integer lo(2), hi(2)
      integer f_l1, f_l2, f_h1, f_h2
      integer m_l1, m_l2, m_h1, m_h2
      
      real(amrex_real)   den(f_l1:f_h1, f_l2:f_h2)
      integer mask(m_l1:m_h1, m_l2:m_h2)

      integer cdir, bct, maxorder
      real(amrex_real) bcl
      real(amrex_real) h(2)

      integer i, j
      logical is_dirichlet
      logical is_neumann

      integer m, lenx, leny

      integer Lmaxorder
      integer maxmaxorder
      parameter(maxmaxorder=4)
      real(amrex_real) x(-1:maxmaxorder-2)
      real(amrex_real) coef(-1:maxmaxorder-2)
      real(amrex_real) xInt

      is_dirichlet(i) = ( i .eq. LO_DIRICHLET )
      is_neumann(i)   = ( i .eq. LO_NEUMANN )

      if ( maxorder .eq. -1 ) then
         Lmaxorder = maxmaxorder
      else
         Lmaxorder = MIN(maxorder,maxmaxorder)
      end if
      lenx = MIN(hi(1)-lo(1), Lmaxorder-2)
      leny = MIN(hi(2)-lo(2), Lmaxorder-2)
!     
!     The Left face of the grid
!
      if(cdir .eq. 0) then
         if (is_neumann(bct)) then
            do j = lo(2), hi(2)
               den(lo(1),j) = 1.0D0
            end do
         else if (is_dirichlet(bct)) then
            do m=0,lenx
               x(m) = m + 0.5D0
            end do
            x(-1) = - bcl/h(1)
            xInt = - 0.5D0
            call polyInterpCoeff(xInt, x, lenx+2, coef)
            do j = lo(2), hi(2)
               den(lo(1),j) = merge(coef(0), 0.0D0, &
                   mask(lo(1)-1,j) .gt. 0)
            end do
            
         else if ( bct .eq. LO_REFLECT_ODD ) then

            do j = lo(2), hi(2)
               den(lo(1),j) = merge(-1.0D0, 0.0D0, &
                   mask(lo(1)-1,j) .gt. 0)
            end do

         else 
            print *,'UNKNOWN BC ON LEFT FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
!     
!     The Right face of the grid
!
      if(cdir .eq. 2) then
         if(is_neumann(bct)) then
            do j = lo(2), hi(2)
               den(hi(1),j) = 1.0D0
            end do
         else if (is_dirichlet(bct)) then
            do m=0,lenx
               x(m) = m + 0.5D0
            end do
            x(-1) = - bcl/h(1)
            xInt = - 0.5D0
            call polyInterpCoeff(xInt, x, lenx+2, coef)
            do j = lo(2), hi(2)
               den(hi(1),j)   = merge(coef(0), 0.0D0, &
                   mask(hi(1)+1,j) .gt. 0)
            end do

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do j = lo(2), hi(2)
               den(hi(1),j) = merge(-1.0D0, 0.0D0, &
                   mask(hi(1)+1,j) .gt. 0)
            end do

         else
            print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
!
!     The Bottom of the Grid
!
      if(cdir .eq. 1) then
         if(is_neumann(bct)) then
            do i = lo(1),hi(1)
               den(i,lo(2))   = 1.0D0
            end do
         else if (is_dirichlet(bct)) then
            do m=0,leny
               x(m) = m + 0.5D0
            end do
            x(-1) = - bcl/h(2)
            xInt = - 0.5D0
            call polyInterpCoeff(xInt, x, leny+2, coef)
            do i = lo(1), hi(1)
               den(i, lo(2))   = merge(coef(0), 0.0D0, &
                   mask(i, lo(2)-1) .gt. 0)
            end do

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do i = lo(1), hi(1)
               den(i,lo(2)) = merge(-1.0D0, 0.0D0, &
                   mask(i,lo(2)-1) .gt. 0)
            end do

         else
            print *,'UNKNOWN BC ON BOTTOM FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
!
!     The top of the grid
!
      if (cdir .eq. 3) then
         if(is_neumann(bct)) then
            do i = lo(1), hi(1)
               den(i,hi(2))   = 1.0D0
            end do
         else if (is_dirichlet(bct)) then
            if ( bct .eq. LO_REFLECT_ODD ) leny = 0
            do m=0,leny
               x(m) = m + 0.5D0
            end do
            x(-1) = - bcl/h(2)
            xInt = - 0.5D0
            call polyInterpCoeff(xInt, x, leny+2, coef)
            do i = lo(1), hi(1)
               den(i,hi(2))   = merge(coef(0), 0.0D0, &
                   mask(i,hi(2)+1) .gt. 0)
            end do

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do i = lo(1), hi(1)
               den(i,hi(2)) = merge(-1.0D0, 0.0D0, &
                   mask(i,hi(2)+1) .gt. 0)
            end do

         else
            print *,'UNKNOWN BC ON TOP FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if

      end

end module amrex_abec2_module

