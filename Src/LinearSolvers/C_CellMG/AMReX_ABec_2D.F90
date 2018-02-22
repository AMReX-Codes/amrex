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
!     Gauss-Seidel Red-Black (GSRB):
!     Apply the GSRB relaxation to the state phi for the equation
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
    subroutine FORT_GSRB ( &
           phi,DIMS(phi), &
           rhs,DIMS(rhs), &
           alpha, beta, &
           a,  DIMS(a), &
           bX, DIMS(bX), &
           bY, DIMS(bY), &
           f0, DIMS(f0), &
           m0, DIMS(m0), &
           f1, DIMS(f1), &
           m1, DIMS(m1), &
           f2, DIMS(f2), &
           m2, DIMS(m2), &
           f3, DIMS(f3), &
           m3, DIMS(m3), &
           lo,hi,blo,bhi, &
           nc,h,redblack &
           )

      implicit none

      REAL_T alpha, beta
      integer DIMDEC(phi)
      integer DIMDEC(rhs)
      integer DIMDEC(a)
      integer DIMDEC(bX)
      integer DIMDEC(bY)
      integer  lo(BL_SPACEDIM),  hi(BL_SPACEDIM)
      integer blo(BL_SPACEDIM), bhi(BL_SPACEDIM)
      integer nc
      integer redblack
      integer DIMDEC(f0)
      REAL_T f0(DIMV(f0))
      integer DIMDEC(f1)
      REAL_T f1(DIMV(f1))
      integer DIMDEC(f2)
      REAL_T f2(DIMV(f2))
      integer DIMDEC(f3)
      REAL_T f3(DIMV(f3))
      integer DIMDEC(m0)
      integer m0(DIMV(m0))
      integer DIMDEC(m1)
      integer m1(DIMV(m1))
      integer DIMDEC(m2)
      integer m2(DIMV(m2))
      integer DIMDEC(m3)
      integer m3(DIMV(m3))
      REAL_T  h(BL_SPACEDIM)
      REAL_T   phi(DIMV(phi),nc)
      REAL_T   rhs(DIMV(rhs),nc)
      REAL_T     a(DIMV(a))
      REAL_T    bX(DIMV(bX))
      REAL_T    bY(DIMV(bY))

      integer  i, j, ioff, joff, n

      REAL_T dhx, dhy, cf0, cf1, cf2, cf3
      REAL_T delta, gamma, rho, rho_x, rho_y

      integer LSDIM
      parameter(LSDIM=127)
      REAL_T a_ls(0:LSDIM)
      REAL_T b_ls(0:LSDIM)
      REAL_T c_ls(0:LSDIM)
      REAL_T r_ls(0:LSDIM)
      REAL_T u_ls(0:LSDIM)

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
         print *,'BOGUS DO_LINE '
         call bl_error("stop")
       end if
      end do

    end subroutine FORT_GSRB

!-----------------------------------------------------------------------
!      
!     JACOBI:
!     Apply the JACOBI relaxation to the state phi for the equation
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
    subroutine FORT_JACOBI ( &
           phi,DIMS(phi), &
           rhs,DIMS(rhs), &
           alpha, beta, &
           a,  DIMS(a), &
           bX, DIMS(bX), &
           bY, DIMS(bY), &
           f0, DIMS(f0), &
           m0, DIMS(m0), &
           f1, DIMS(f1), &
           m1, DIMS(m1), &
           f2, DIMS(f2), &
           m2, DIMS(m2), &
           f3, DIMS(f3), &
           m3, DIMS(m3), &
           lo,hi,nc, &
           h &
           )

      implicit none

      REAL_T alpha, beta
      integer DIMDEC(phi)
      integer DIMDEC(rhs)
      integer DIMDEC(a)
      integer DIMDEC(bX)
      integer DIMDEC(bY)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer nc
      integer DIMDEC(f0)
      REAL_T f0(DIMV(f0))
      integer DIMDEC(f1)
      REAL_T f1(DIMV(f1))
      integer DIMDEC(f2)
      REAL_T f2(DIMV(f2))
      integer DIMDEC(f3)
      REAL_T f3(DIMV(f3))
      integer DIMDEC(m0)
      integer m0(DIMV(m0))
      integer DIMDEC(m1)
      integer m1(DIMV(m1))
      integer DIMDEC(m2)
      integer m2(DIMV(m2))
      integer DIMDEC(m3)
      integer m3(DIMV(m3))
      REAL_T  h(BL_SPACEDIM)
      REAL_T   phi(DIMV(phi),nc)
      REAL_T   rhs(DIMV(rhs),nc)
      REAL_T     a(DIMV(a))
      REAL_T    bX(DIMV(bX))
      REAL_T    bY(DIMV(bY))

      integer  i, j, n

      REAL_T dhx, dhy, cf0, cf1, cf2, cf3
      REAL_T delta, gamma, rho

      REAL_T, allocatable :: phinew(:,:)

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

    end subroutine FORT_JACOBI

!-----------------------------------------------------------------------
!
!     Fill in a matrix x vector operator here
!
    subroutine FORT_ADOTX( &
           y,DIMS(y), &
           x,DIMS(x), &
           alpha, beta, &
           a, DIMS(a), &
           bX,DIMS(bX), &
           bY,DIMS(bY), &
           lo,hi,nc, &
           h &
           )

      implicit none

      REAL_T alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer DIMDEC(y)
      integer DIMDEC(x)
      integer DIMDEC(a)
      integer DIMDEC(bX)
      integer DIMDEC(bY)
      REAL_T  y(DIMV(y),nc)
      REAL_T  x(DIMV(x),nc)
      REAL_T  a(DIMV(a))
      REAL_T bX(DIMV(bX))
      REAL_T bY(DIMV(bY))
      REAL_T h(BL_SPACEDIM)

      integer i,j,n
      REAL_T dhx,dhy

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

    end subroutine FORT_ADOTX

!-----------------------------------------------------------------------
!
!     Fill in a matrix x vector operator here
!
    subroutine FORT_NORMA( &
           res, &
           alpha, beta, &
           a, DIMS(a), &
           bX,DIMS(bX), &
           bY,DIMS(bY), &
           lo,hi,nc, &
           h &
           )

      implicit none

      REAL_T res
      REAL_T alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer DIMDEC(a)
      integer DIMDEC(bX)
      integer DIMDEC(bY)
      REAL_T  a(DIMV(a))
      REAL_T bX(DIMV(bX))
      REAL_T bY(DIMV(bY))
      REAL_T h(BL_SPACEDIM)

      integer i,j,n
      REAL_T dhx,dhy

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

    end subroutine FORT_NORMA

!-----------------------------------------------------------------------
!
!     Fill in fluxes
!
    subroutine FORT_FLUX( &
           x,DIMS(x), &
           alpha, beta, &
           a, DIMS(a), &
           bX,DIMS(bX), &
           bY,DIMS(bY), &
           xlo,xhi, &
           ylo,yhi, &
           nc, &
           h, &
           xflux,DIMS(xflux), &
           yflux,DIMS(yflux) &
           )

      implicit none

      REAL_T alpha, beta
      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM)
      integer ylo(BL_SPACEDIM), yhi(BL_SPACEDIM)
      integer nc
      integer DIMDEC(x)
      integer DIMDEC(a)
      integer DIMDEC(bX)
      integer DIMDEC(bY)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      REAL_T  x(DIMV(x),nc)
      REAL_T  a(DIMV(a))
      REAL_T bX(DIMV(bX))
      REAL_T bY(DIMV(bY))
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T h(BL_SPACEDIM)

      REAL_T dhx, dhy
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

    end subroutine FORT_FLUX
