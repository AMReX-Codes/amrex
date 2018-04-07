#include <AMReX_REAL.H>
#include "AMReX_ArrayLim.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_LO_BCTYPES.H"

      subroutine ab2_gsrb (lo, hi, blo, bhi,
     &     phi, phi_l1, phi_l2, phi_h1, phi_h2,
     &     res, res_l1, res_l2, res_h1, res_h2,
     &     alpha, beta,
     &     a,   a_l1,   a_l2,   a_h1,   a_h2,
     &     bX,  bX_l1,  bX_l2,  bX_h1,  bX_h2,
     &     bY,  bY_l1,  bY_l2,  bY_h1,  bY_h2,
     &     f0,  f0_l1,  f0_l2,  f0_h1,  f0_h2,
     &     m0,  m0_l1,  m0_l2,  m0_h1,  m0_h2,
     &     f1,  f1_l1,  f1_l2,  f1_h1,  f1_h2,
     &     m1,  m1_l1,  m1_l2,  m1_h1,  m1_h2,
     &     f2,  f2_l1,  f2_l2,  f2_h1,  f2_h2,
     &     m2,  m2_l1,  m2_l2,  m2_h1,  m2_h2,
     &     f3,  f3_l1,  f3_l2,  f3_h1,  f3_h2,
     &     m3,  m3_l1,  m3_l2,  m3_h1,  m3_h2,
     &     nc, h, redBlack)

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

      REAL_T  phi(phi_l1:phi_h1, phi_l2:phi_h2,1:nc)
      REAL_T  res(res_l1:res_h1, res_l2:res_h2,1:nc)
      REAL_T  alpha, beta

      REAL_T   a( a_l1:a_h1,  a_l2:a_h2)
      REAL_T  bX(bX_l1:bX_h1, bX_l2:bX_h2)
      REAL_T  bY(bY_l1:bY_h1, bY_l2:bY_h2)
      REAL_T  f0(f0_l1:f0_h1, f0_l2:f0_h2)
      REAL_T  f1(f1_l1:f1_h1, f1_l2:f1_h2)
      REAL_T  f2(f2_l1:f2_h1, f2_l2:f2_h2)
      REAL_T  f3(f3_l1:f3_h1, f3_l2:f3_h2)
      integer  m0(m0_l1:m0_h1, m0_l2:m0_h2)
      integer  m1(m1_l1:m1_h1, m1_l2:m1_h2)
      integer  m2(m2_l1:m2_h1, m2_l2:m2_h2)
      integer  m3(m3_l1:m3_h1, m3_l2:m3_h2)
      REAL_T  h(2)
c
      integer  i, j, ioff, n
c
      REAL_T dhx, dhy, cf0, cf1, cf2, cf3
      REAL_T delta, gamma, rho_x, rho_y

      integer LSDIM
      parameter(LSDIM=127)
      REAL_T a_ls(0:LSDIM)
      REAL_T b_ls(0:LSDIM)
      REAL_T c_ls(0:LSDIM)
      REAL_T r_ls(0:LSDIM)
      REAL_T u_ls(0:LSDIM)

      integer do_line
      integer ilen,jlen
      
      if (h(2). gt. 1.5D0*h(1)) then 
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

c
      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      do n = 1, nc
       if (do_line .eq. 0) then
         do j = lo(2), hi(2)
            ioff = MOD(j + redblack,2)
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
               delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2)
     $              +  dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)

               gamma = alpha*a(i,j)
     $              +   dhx*( bX(i,j) + bX(i+1,j) )
     $              +   dhy*( bY(i,j) + bY(i,j+1) )
c     
c     
               phi(i,j,n) = phi(i,j,n) + res(i,j,n)/(gamma - delta)
c     
            end do
         end do
       else if (do_line .eq. 2) then

          print *,'Need to fix line solver in ABec2_2D.F'
          call bl_abort()

         do i = lo(1) + redblack,hi(1),2
             do j = lo(2), hi(2)
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
               delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2)
     $               + dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)
c     
               gamma = alpha*a(i,j)
     $              +   dhx*( bX(i,j) + bX(i+1,j) )
     $              +   dhy*( bY(i,j) + bY(i,j+1) )
c     
               rho_x = dhx*(bX(i,j)*phi(i-1,j,n) + bX(i+1,j)*phi(i+1,j,n))

               a_ls(j-lo(2)) = -dhy*bY(i,j)
               b_ls(j-lo(2)) = gamma - delta
               c_ls(j-lo(2)) = -dhy*bY(i,j+1)
               r_ls(j-lo(2)) = res(i,j,n) + rho_x - phi(i,j,n)*delta

               if (j .eq. lo(2)) 
     $            r_ls(j-lo(2)) = r_ls(j-lo(2)) + dhy*bY(i,j)*phi(i,j-1,n)

               if (j .eq. hi(2)) 
     $            r_ls(j-lo(2)) = r_ls(j-lo(2)) + dhy*bY(i,j+1)*phi(i,j+1,n)

             end do

             call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,jlen)
c     
             do j = lo(2), hi(2)
               phi(i,j,n) = u_ls(j-lo(2))
             end do
         end do

       else if (do_line .eq. 1) then

          print *,'Need to fix line solver in ABec2_2D.F'
          call bl_abort()

           do j = lo(2) + redblack,hi(2),2
             do i = lo(1), hi(1)
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
               delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2)
     $               + dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)
c     
               gamma = alpha*a(i,j)
     $              +   dhx*( bX(i,j) + bX(i+1,j) )
     $              +   dhy*( bY(i,j) + bY(i,j+1) )
c     
               rho_y = dhy*(bY(i,j)*phi(i,j-1,n) + bY(i,j+1)*phi(i,j+1,n))

               a_ls(i-lo(1)) = -dhx*bX(i,j)
               b_ls(i-lo(1)) = gamma - delta
               c_ls(i-lo(1)) = -dhx*bX(i+1,j)
               r_ls(i-lo(1)) = res(i,j,n) + rho_y - phi(i,j,n)*delta

               if (i .eq. lo(1)) 
     $            r_ls(i-lo(1)) = r_ls(i-lo(1)) + dhx*bX(i,j)*phi(i-1,j,n)

               if (i .eq. hi(1)) 
     $            r_ls(i-lo(1)) = r_ls(i-lo(1)) + dhx*bX(i+1,j)*phi(i+1,j,n)
             end do

             call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,ilen)
c     
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


      subroutine ab2_jacobi (lo, hi,
     &     phi, phi_l1, phi_l2, phi_h1, phi_h2,
     &     res, res_l1, res_l2, res_h1, res_h2,
     &     alpha, beta,
     &     a,   a_l1,   a_l2,   a_h1,   a_h2,
     &     bX,  bX_l1,  bX_l2,  bX_h1,  bX_h2,
     &     bY,  bY_l1,  bY_l2,  bY_h1,  bY_h2,
     &     f0,  f0_l1,  f0_l2,  f0_h1,  f0_h2,
     &     m0,  m0_l1,  m0_l2,  m0_h1,  m0_h2,
     &     f1,  f1_l1,  f1_l2,  f1_h1,  f1_h2,
     &     m1,  m1_l1,  m1_l2,  m1_h1,  m1_h2,
     &     f2,  f2_l1,  f2_l2,  f2_h1,  f2_h2,
     &     m2,  m2_l1,  m2_l2,  m2_h1,  m2_h2,
     &     f3,  f3_l1,  f3_l2,  f3_h1,  f3_h2,
     &     m3,  m3_l1,  m3_l2,  m3_h1,  m3_h2,
     &     nc, h)

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

      REAL_T  phi(phi_l1:phi_h1, phi_l2:phi_h2,1:nc)
      REAL_T  res(res_l1:res_h1, res_l2:res_h2,1:nc)
      REAL_T  alpha, beta

      REAL_T   a( a_l1:a_h1,  a_l2:a_h2)
      REAL_T  bX(bX_l1:bX_h1, bX_l2:bX_h2)
      REAL_T  bY(bY_l1:bY_h1, bY_l2:bY_h2)
      REAL_T  f0(f0_l1:f0_h1, f0_l2:f0_h2)
      REAL_T  f1(f1_l1:f1_h1, f1_l2:f1_h2)
      REAL_T  f2(f2_l1:f2_h1, f2_l2:f2_h2)
      REAL_T  f3(f3_l1:f3_h1, f3_l2:f3_h2)
      integer  m0(m0_l1:m0_h1, m0_l2:m0_h2)
      integer  m1(m1_l1:m1_h1, m1_l2:m1_h2)
      integer  m2(m2_l1:m2_h1, m2_l2:m2_h2)
      integer  m3(m3_l1:m3_h1, m3_l2:m3_h2)
      REAL_T  h(2)
c
      integer  i, j, n
c
      REAL_T dhx, dhy, cf0, cf1, cf2, cf3
      REAL_T delta, gamma
      REAL_T, allocatable :: phit(:,:)
      allocate(phit(lo(1):hi(1),lo(2):hi(2)))
c
      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
c
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1),hi(1)
c     
               cf0 = merge(f0(lo(1),j), 0.0D0,
     $              (i .eq. lo(1)) .and. (m0(lo(1)-1,j).gt.0))
               cf1 = merge(f1(i,lo(2)), 0.0D0,
     $              (j .eq. lo(2)) .and. (m1(i,lo(2)-1).gt.0))
               cf2 = merge(f2(hi(1),j), 0.0D0,
     $              (i .eq. hi(1)) .and. (m2(hi(1)+1,j).gt.0))
               cf3 = merge(f3(i,hi(2)), 0.0D0,
     $              (j .eq. hi(2)) .and. (m3(i,hi(2)+1).gt.0))
c     
               delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2)
     $              +  dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)
c     
               gamma = alpha*a(i,j)
     $              +   dhx*( bX(i,j) + bX(i+1,j) )
     $              +   dhy*( bY(i,j) + bY(i,j+1) )
c     
               phit = phi(i,j,n) + res(i,j,n) / (gamma - delta)
c     
            end do
         end do
         phi(lo(1):hi(1),lo(2):hi(2),n) = phit         
      end do
      deallocate(phit)
      end


      subroutine ab2_bndryrlx(lo, hi,
     &     den,  f_l1, f_l2, f_h1, f_h2,
     &     mask, m_l1, m_l2, m_h1, m_h2,
     &     cdir, bct, bcl, maxorder, h)

      use amrex_lo_util_module, only : polyInterpCoeff
      implicit none
      integer lo(2), hi(2)
      integer f_l1, f_l2, f_h1, f_h2
      integer m_l1, m_l2, m_h1, m_h2
      
      REAL_T   den(f_l1:f_h1, f_l2:f_h2)
      integer mask(m_l1:m_h1, m_l2:m_h2)

      integer cdir, bct, maxorder
      REAL_T bcl
      REAL_T h(2)
c
      integer i, j
      logical is_dirichlet
      logical is_neumann
c
      integer m, lenx, leny
c
      integer Lmaxorder
      integer maxmaxorder
      parameter(maxmaxorder=4)
      REAL_T x(-1:maxmaxorder-2)
      REAL_T coef(-1:maxmaxorder-2)
      REAL_T xInt
c
      is_dirichlet(i) = ( i .eq. LO_DIRICHLET )
      is_neumann(i)   = ( i .eq. LO_NEUMANN )
c
      if ( maxorder .eq. -1 ) then
         Lmaxorder = maxmaxorder
      else
         Lmaxorder = MIN(maxorder,maxmaxorder)
      end if
      lenx = MIN(hi(1)-lo(1), Lmaxorder-2)
      leny = MIN(hi(2)-lo(2), Lmaxorder-2)
c     
c     The Left face of the grid
c
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
               den(lo(1),j) = merge(coef(0), 0.0D0,
     $              mask(lo(1)-1,j) .gt. 0)
            end do
            
         else if ( bct .eq. LO_REFLECT_ODD ) then

            do j = lo(2), hi(2)
               den(lo(1),j) = merge(-1.0D0, 0.0D0,
     $              mask(lo(1)-1,j) .gt. 0)
            end do

         else 
            print *,'UNKNOWN BC ON LEFT FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
c     
c     The Right face of the grid
c
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
               den(hi(1),j)   = merge(coef(0), 0.0D0,
     $              mask(hi(1)+1,j) .gt. 0)
            end do

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do j = lo(2), hi(2)
               den(hi(1),j) = merge(-1.0D0, 0.0D0,
     $              mask(hi(1)+1,j) .gt. 0)
            end do

         else
            print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
c
c     The Bottom of the Grid
c
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
               den(i, lo(2))   = merge(coef(0), 0.0D0,
     $              mask(i, lo(2)-1) .gt. 0)
            end do

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do i = lo(1), hi(1)
               den(i,lo(2)) = merge(-1.0D0, 0.0D0,
     $              mask(i,lo(2)-1) .gt. 0)
            end do

         else
            print *,'UNKNOWN BC ON BOTTOM FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
c
c     The top of the grid
c
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
               den(i,hi(2))   = merge(coef(0), 0.0D0,
     $              mask(i,hi(2)+1) .gt. 0)
            end do

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do i = lo(1), hi(1)
               den(i,hi(2)) = merge(-1.0D0, 0.0D0,
     $              mask(i,hi(2)+1) .gt. 0)
            end do

         else
            print *,'UNKNOWN BC ON TOP FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
c
      end
