c
c $Id: ABec_2D.F,v 1.1 1998-03-24 07:06:46 almgren Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>

#include "ABec_F.H"
#include "ArrayLim.H"

c-----------------------------------------------------------------------
c      
c     Gauss-Seidel Red-Black (GSRB):
c     Apply the GSRB relaxation to the state phi for the equation
c     L(phi) = alpha*a(x)*phi(x) - beta*Div(b(x)Grad(phi(x))) = rhs(x)
c     central differenced, according to the arrays of boundary
c     masks (m#) and auxiliary data (f#).
c     
c     In general, if the linear operator L=gamma*y-rho, the GS relaxation
c     is y = (R - rho)/gamma.  Near a boundary, the ghost data is filled
c     using a polynomial interpolant based on the "old" phi values, so
c     L=(gamma-delta)*y - rho + delta*yOld.   The resulting iteration is
c     
c     y = (R - delta*yOld + rho)/(gamma - delta)
c     
c     This expression is valid additionally in the interior provided
c     delta->0 there.  delta is constructed by summing all the
c     contributions to the central stencil element coming from boundary 
c     interpolants.  The f#s contain the corresponding coefficient of 
c     the interpolating polynomial.  The masks are set > 0 if the boundary 
c     value was filled with an interpolant involving the central stencil 
c     element.
c     
c-----------------------------------------------------------------------
      subroutine FORT_GSRB (
     $     phi,DIMS(phi),
     $     rhs,DIMS(rhs),
     $     alpha, beta,
     $     a,  DIMS(a),
     $     bX, DIMS(bX),
     $     bY, DIMS(bY),
     $     f0, DIMS(f0),
     $     m0, DIMS(m0),
     $     f1, DIMS(f1),
     $     m1, DIMS(m1),
     $     f2, DIMS(f2),
     $     m2, DIMS(m2),
     $     f3, DIMS(f3),
     $     m3, DIMS(m3),
     $     lo,hi,nc,
     $     h,redblack
     $     )
      REAL_T alpha, beta
      integer DIMDEC(phi)
      integer DIMDEC(rhs)
      integer DIMDEC(a)
      integer DIMDEC(bX)
      integer DIMDEC(bY)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
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
c
      integer  i, j, ioff, n
c
      REAL_T dhx, dhy, cf0, cf1, cf2, cf3
      REAL_T delta, gamma, rho, rho_x, rho_y

      REAL_T a_ls(0:4096)
      REAL_T b_ls(0:4096)
      REAL_T c_ls(0:4096)
      REAL_T r_ls(0:4096)
      REAL_T u_ls(0:4096)

      integer do_line
      integer ilen,jlen
      
      if (h(2). gt. 1.5*h(1)) then 
        do_line = 1
        ilen = hi(1)-lo(1)+1
        if (ilen .gt. 4096) then
          print *,'TOO BIG FOR LINE SOLVE IN GSRB: ilen = ',ilen
          stop
        end if
      else if (h(1) .gt. 1.5*h(2)) then
        do_line = 2
        jlen = hi(2)-lo(2)+1
        if (jlen .gt. 4096) then
          print *,'TOO BIG FOR LINE SOLVE IN GSRB: jlen = ',jlen
          stop
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
               cf0 = cvmgt(f0(lo(1),j), 0.0,
     $              (i .eq. lo(1)) .and. (m0(lo(1)-1,j).gt.0))
               cf1 = cvmgt(f1(i,lo(2)), 0.0,
     $              (j .eq. lo(2)) .and. (m1(i,lo(2)-1).gt.0))
               cf2 = cvmgt(f2(hi(1),j), 0.0,
     $              (i .eq. hi(1)) .and. (m2(hi(1)+1,j).gt.0))
               cf3 = cvmgt(f3(i,hi(2)), 0.0,
     $              (j .eq. hi(2)) .and. (m3(i,hi(2)+1).gt.0))
c     
               delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2)
     $              +  dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)
c     
               gamma = alpha*a(i,j)
     $              +   dhx*( bX(i,j) + bX(i+1,j) )
     $              +   dhy*( bY(i,j) + bY(i,j+1) )
c     
               rho = dhx*(bX(i,j)*phi(i-1,j,n) + bX(i+1,j)*phi(i+1,j,n))
     $              +dhy*(bY(i,j)*phi(i,j-1,n) + bY(i,j+1)*phi(i,j+1,n))
c     
               phi(i,j,n) = (rhs(i,j,n) + rho - phi(i,j,n)*delta)
     $              /                (gamma - delta)
c     
            end do
         end do
       else if (do_line .eq. 2) then
         do i = lo(1) + redblack,hi(1),2
             do j = lo(2), hi(2)
c     
               cf0 = cvmgt(f0(lo(1),j), 0.0,
     $              (i .eq. lo(1)) .and. (m0(lo(1)-1,j).gt.0))
               cf1 = cvmgt(f1(i,lo(2)), 0.0,
     $              (j .eq. lo(2)) .and. (m1(i,lo(2)-1).gt.0))
               cf2 = cvmgt(f2(hi(1),j), 0.0,
     $              (i .eq. hi(1)) .and. (m2(hi(1)+1,j).gt.0))
               cf3 = cvmgt(f3(i,hi(2)), 0.0,
     $              (j .eq. hi(2)) .and. (m3(i,hi(2)+1).gt.0))
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
               r_ls(j-lo(2)) = rhs(i,j,n) + rho_x - phi(i,j,n)*delta

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

           do j = lo(2) + redblack,hi(2),2
             do i = lo(1), hi(1)
c     
               cf0 = cvmgt(f0(lo(1),j), 0.0,
     $              (i .eq. lo(1)) .and. (m0(lo(1)-1,j).gt.0))
               cf1 = cvmgt(f1(i,lo(2)), 0.0,
     $              (j .eq. lo(2)) .and. (m1(i,lo(2)-1).gt.0))
               cf2 = cvmgt(f2(hi(1),j), 0.0,
     $              (i .eq. hi(1)) .and. (m2(hi(1)+1,j).gt.0))
               cf3 = cvmgt(f3(i,hi(2)), 0.0,
     $              (j .eq. hi(2)) .and. (m3(i,hi(2)+1).gt.0))
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
               r_ls(i-lo(1)) = rhs(i,j,n) + rho_y - phi(i,j,n)*delta

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
         stop
       end if
      end do

      end
      
c-----------------------------------------------------------------------
c
c     Solve Preconditioned system here
c
      subroutine FORT_CGPRECND(
     $     zz,DIMS(zz),
     $     rho,
     $     rr,DIMS(rr),
     $     alpha, beta,
     $     a, DIMS(a),
     $     bX,DIMS(bX),
     $     bY,DIMS(bY),
     $     lo,hi,nc,
     $     h
     $     )
      REAL_T alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer DIMDEC(zz)
      integer DIMDEC(rr)
      integer DIMDEC(a)
      integer DIMDEC(bX)
      integer DIMDEC(bY)
      REAL_T zz(DIMV(zz),nc)
      REAL_T rr(DIMV(rr),nc)
      REAL_T  a(DIMV(a))
      REAL_T bX(DIMV(bX))
      REAL_T bY(DIMV(bY))
      REAL_T h(BL_SPACEDIM), rho
c
      integer i,j, n
      REAL_T de, dhx, dhy
c
      rho = 0.0d0
      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               de = 1.0/(alpha*a(i,j) + dhx*(bX(i+1,j)+bX(i,j))
     $              + dhy*(bY(i,j+1)+bY(i,j))
     $              )
               zz(i,j,n) = rr(i,j,n)*de
               rho = rho + zz(i,j,n)*rr(i,j,n)
            end do
         end do
      end do
c
      end
c-----------------------------------------------------------------------
c
c     Fill in a matrix x vector operator here
c
      subroutine FORT_ADOTX(
     $     y,DIMS(y),
     $     x,DIMS(x),
     $     alpha, beta,
     $     a, DIMS(a),
     $     bX,DIMS(bX),
     $     bY,DIMS(bY),
     $     lo,hi,nc,
     $     h
     $     )
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
c
      integer i,j,n
      REAL_T uxx,uyy,dhx,dhy
c
      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
c
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               y(i,j,n) = alpha*a(i,j)*x(i,j,n)
     $              - dhx*
     $              (   bX(i+1,j)*( x(i+1,j,n) - x(i  ,j,n) )
     $              -   bX(i  ,j)*( x(i  ,j,n) - x(i-1,j,n) ) )
     $              - dhy*
     $              (   bY(i,j+1)*( x(i,j+1,n) - x(i,j  ,n) )
     $              -   bY(i,j  )*( x(i,j  ,n) - x(i,j-1,n) ) )
            end do
         end do
      end do
      end

c-----------------------------------------------------------------------
c
c     Tridiagonal solve
c
      subroutine tridiag(a,b,c,r,u,n)
 
      integer n
      integer nmax
  
      REAL_T a(n)
      REAL_T b(n)
      REAL_T c(n)
      REAL_T r(n)
      REAL_T u(n)

      parameter (nmax = 4098)

      integer j
      REAL_T bet
      REAL_T gam(nmax)
      if (b(1) .eq. 0) print *,'CANT HAVE B(1) = ZERO'

      bet = b(1)
      u(1) = r(1)/bet

      do j = 2,n
        gam(j) = c(j-1)/bet
        bet = b(j) - a(j)*gam(j)
        if (bet .eq. 0) then
          print *,'TRIDIAG FAILED '
          stop
        end if
        u(j) = (r(j)-a(j)*u(j-1))/bet
      end do

      do j = n-1,1,-1
        u(j) = u(j) - gam(j+1)*u(j+1)
      end do

      return
      end

