c
c $Id: ABec_3D.F,v 1.6 2000-08-26 19:41:00 car Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>

#include "ABec_F.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"

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
c     L=(gamma-delta)*y - rho + delta*yOld.  The resulting iteration is
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
     $     bZ, DIMS(bZ),
     $     f0, DIMS(f0),
     $     m0, DIMS(m0),
     $     f1, DIMS(f1),
     $     m1, DIMS(m1),
     $     f2, DIMS(f2),
     $     m2, DIMS(m2),
     $     f3, DIMS(f3),
     $     m3, DIMS(m3),
     $     f4, DIMS(f4),
     $     m4, DIMS(m4),
     $     f5, DIMS(f5),
     $     m5, DIMS(m5),
     $     lo,hi,nc,
     $     h,redblack
     $     )
      REAL_T alpha, beta
      integer DIMDEC(phi)
      integer DIMDEC(rhs)
      integer DIMDEC(a)
      integer DIMDEC(bX)
      integer DIMDEC(bY)
      integer DIMDEC(bZ)
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
      integer DIMDEC(f4)
      REAL_T f4(DIMV(f4))
      integer DIMDEC(f5)
      REAL_T f5(DIMV(f5))
      integer DIMDEC(m0)
      integer m0(DIMV(m0))
      integer DIMDEC(m1)
      integer m1(DIMV(m1))
      integer DIMDEC(m2)
      integer m2(DIMV(m2))
      integer DIMDEC(m3)
      integer m3(DIMV(m3))
      integer DIMDEC(m4)
      integer m4(DIMV(m4))
      integer DIMDEC(m5)
      integer m5(DIMV(m5))
      REAL_T  h(BL_SPACEDIM)
      REAL_T   phi(DIMV(phi),nc)
      REAL_T   rhs(DIMV(rhs),nc)
      REAL_T     a(DIMV(a))
      REAL_T    bX(DIMV(bX))
      REAL_T    bY(DIMV(bY))
      REAL_T    bZ(DIMV(bZ))
c
      integer  i, j, k, ioff, n
c
      REAL_T dhx, dhy, dhz, cf0, cf1, cf2, cf3, cf4, cf5
      REAL_T delta, gamma, rho, rho_xy, rho_yz, rho_xz
c
      integer do_line
      integer LSDIM
      parameter (LSDIM=127)
      REAL_T a_ls(0:LSDIM)
      REAL_T b_ls(0:LSDIM)
      REAL_T c_ls(0:LSDIM)
      REAL_T r_ls(0:LSDIM)
      REAL_T u_ls(0:LSDIM)
      integer ilen,jlen,klen

      ilen = hi(1)-lo(1)+1
      jlen = hi(2)-lo(2)+1
      klen = hi(3)-lo(3)+1

#if 0
      if (klen .gt. ilen .and. klen .gt. jlen) then
        do_line = 3
        if (klen .gt. LSDIM) then
          print *,'TOO BIG FOR LINE SOLVE IN GSRB: klen = ',klen
          stop
        endif
      else if (jlen .gt. ilen .and. jlen .gt. klen) then
        do_line = 2
        if (jlen .gt. LSDIM) then
          print *,'JLEN TOO BIG FOR LINE SOLVE IN GSRB: jlen = ',jlen
          stop
        endif
      else if (ilen .gt. jlen .and. ilen .gt. klen) then
        do_line = 1
        if (ilen .gt. LSDIM) then
          print *,'ILEN TOO BIG FOR LINE SOLVE IN GSRB: ilen = ',ilen
          stop
        endif
      else 
        do_line = 0
      endif
#endif

      do_line = 0

      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      dhz = beta/h(3)**2

      do n = 1, nc
         if (do_line .eq. 0) then
          do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               ioff = MOD(j + k + redblack,2)
               do i = lo(1) + ioff,hi(1),2
c
                  cf0 = cvmgt(f0(lo(1),j,k), 0.0,
     $                 (i .eq. lo(1)) .and. (m0(lo(1)-1,j,k).gt.0))
                  cf1 = cvmgt(f1(i,lo(2),k), 0.0,
     $                 (j .eq. lo(2)) .and. (m1(i,lo(2)-1,k).gt.0))
                  cf2 = cvmgt(f2(i,j,lo(3)), 0.0,
     $                 (k .eq. lo(3)) .and. (m2(i,j,lo(3)-1).gt.0))
                  cf3 = cvmgt(f3(hi(1),j,k), 0.0,
     $                 (i .eq. hi(1)) .and. (m3(hi(1)+1,j,k).gt.0))
                  cf4 = cvmgt(f4(i,hi(2),k), 0.0,
     $                 (j .eq. hi(2)) .and. (m4(i,hi(2)+1,k).gt.0))
                  cf5 = cvmgt(f5(i,j,hi(3)), 0.0,
     $                 (k .eq. hi(3)) .and. (m5(i,j,hi(3)+1).gt.0))
c
                  delta = dhx*(bX(i,j,k)*cf0 + bX(i+1,j,k)*cf3)
     $                 +  dhy*(bY(i,j,k)*cf1 + bY(i,j+1,k)*cf4)
     $                 +  dhz*(bZ(i,j,k)*cf2 + bZ(i,j,k+1)*cf5)
c                  
                  gamma = alpha*a(i,j,k)
     $                 +   dhx*(bX(i,j,k)+bX(i+1,j,k))
     $                 +   dhy*(bY(i,j,k)+bY(i,j+1,k))
     $                 +   dhz*(bZ(i,j,k)+bZ(i,j,k+1))
c
                  rho =  dhx*( bX(i  ,j,k)*phi(i-1,j,k,n)
     $                 +       bX(i+1,j,k)*phi(i+1,j,k,n) )
     $                 + dhy*( bY(i,j  ,k)*phi(i,j-1,k,n)
     $                 +       bY(i,j+1,k)*phi(i,j+1,k,n) )
     $                 + dhz*( bZ(i,j,k  )*phi(i,j,k-1,n)
     $                 +       bZ(i,j,k+1)*phi(i,j,k+1,n) )
c
                  phi(i,j,k,n) = (rhs(i,j,k,n)+rho-phi(i,j,k,n)*delta)
     $                 /                   (gamma - delta)
c                  
               end do
            end do
          end do

         else if (do_line .eq. 3) then

           do j = lo(2), hi(2)
             ioff = MOD(j + redblack,2)
             do i = lo(1) + ioff,hi(1),2
c
               do k = lo(3),hi(3)
                  cf0 = cvmgt(f0(lo(1),j,k), 0.0,
     $                 (i .eq. lo(1)) .and. (m0(lo(1)-1,j,k).gt.0))
                  cf1 = cvmgt(f1(i,lo(2),k), 0.0,
     $                 (j .eq. lo(2)) .and. (m1(i,lo(2)-1,k).gt.0))
                  cf2 = cvmgt(f2(i,j,lo(3)), 0.0,
     $                 (k .eq. lo(3)) .and. (m2(i,j,lo(3)-1).gt.0))
                  cf3 = cvmgt(f3(hi(1),j,k), 0.0,
     $                 (i .eq. hi(1)) .and. (m3(hi(1)+1,j,k).gt.0))
                  cf4 = cvmgt(f4(i,hi(2),k), 0.0,
     $                 (j .eq. hi(2)) .and. (m4(i,hi(2)+1,k).gt.0))
                  cf5 = cvmgt(f5(i,j,hi(3)), 0.0,
     $                 (k .eq. hi(3)) .and. (m5(i,j,hi(3)+1).gt.0))
c
                  delta = dhx*(bX(i,j,k)*cf0 + bX(i+1,j,k)*cf3)
     $                 +  dhy*(bY(i,j,k)*cf1 + bY(i,j+1,k)*cf4)
     $                 +  dhz*(bZ(i,j,k)*cf2 + bZ(i,j,k+1)*cf5)
c                  
                  gamma = alpha*a(i,j,k)
     $                 +   dhx*(bX(i,j,k)+bX(i+1,j,k))
     $                 +   dhy*(bY(i,j,k)+bY(i,j+1,k))
     $                 +   dhz*(bZ(i,j,k)+bZ(i,j,k+1))
c
                  rho_xy =  dhx*( bX(i  ,j,k)*phi(i-1,j,k,n)
     $                    +       bX(i+1,j,k)*phi(i+1,j,k,n) )
     $                    + dhy*( bY(i,j  ,k)*phi(i,j-1,k,n)
     $                    +       bY(i,j+1,k)*phi(i,j+1,k,n) )

                  a_ls(k-lo(3)) = -dhz*bZ(i,j,k)
                  b_ls(k-lo(3)) = gamma - delta
                  c_ls(k-lo(3)) = -dhz*bZ(i,j,k+1)
                  r_ls(k-lo(3)) = rhs(i,j,k,n) + rho_xy - phi(i,j,k,n)*delta

                  if (k .eq. lo(3))
     $               r_ls(k-lo(3)) = r_ls(k-lo(3)) + dhz*bZ(i,j,k)*phi(i,j,k-1,n)

                  if (k .eq. hi(3))
     $               r_ls(k-lo(3)) = r_ls(k-lo(3)) + dhz*bZ(i,j,k+1)*phi(i,j,k+1,n)
                end do
c
                call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,klen)

                do k = lo(3),hi(3)
                  phi(i,j,k,n) = u_ls(k-lo(3))
                end do
c                  
               end do
            end do

         else if (do_line .eq. 2) then

           do k = lo(3), hi(3)
             ioff = MOD(k + redblack,2)
             do i = lo(1) + ioff,hi(1),2
c
               do j = lo(2),hi(2)
                  cf0 = cvmgt(f0(lo(1),j,k), 0.0,
     $                 (i .eq. lo(1)) .and. (m0(lo(1)-1,j,k).gt.0))
                  cf1 = cvmgt(f1(i,lo(2),k), 0.0,
     $                 (j .eq. lo(2)) .and. (m1(i,lo(2)-1,k).gt.0))
                  cf2 = cvmgt(f2(i,j,lo(3)), 0.0,
     $                 (k .eq. lo(3)) .and. (m2(i,j,lo(3)-1).gt.0))
                  cf3 = cvmgt(f3(hi(1),j,k), 0.0,
     $                 (i .eq. hi(1)) .and. (m3(hi(1)+1,j,k).gt.0))
                  cf4 = cvmgt(f4(i,hi(2),k), 0.0,
     $                 (j .eq. hi(2)) .and. (m4(i,hi(2)+1,k).gt.0))
                  cf5 = cvmgt(f5(i,j,hi(3)), 0.0,
     $                 (k .eq. hi(3)) .and. (m5(i,j,hi(3)+1).gt.0))
c
                  delta = dhx*(bX(i,j,k)*cf0 + bX(i+1,j,k)*cf3)
     $                 +  dhy*(bY(i,j,k)*cf1 + bY(i,j+1,k)*cf4)
     $                 +  dhz*(bZ(i,j,k)*cf2 + bZ(i,j,k+1)*cf5)
c                  
                  gamma = alpha*a(i,j,k)
     $                 +   dhx*(bX(i,j,k)+bX(i+1,j,k))
     $                 +   dhy*(bY(i,j,k)+bY(i,j+1,k))
     $                 +   dhz*(bZ(i,j,k)+bZ(i,j,k+1))
c
                  rho_xz =  dhx*( bX(i  ,j,k)*phi(i-1,j,k,n)
     $                    +       bX(i+1,j,k)*phi(i+1,j,k,n) )
     $                    + dhz*( bZ(i,j,k  )*phi(i,j,k-1,n)
     $                    +       bZ(i,j,k+1)*phi(i,j,k+1,n) )

                  a_ls(j-lo(2)) = -dhy*bY(i,j,k)
                  b_ls(j-lo(2)) = gamma - delta
                  c_ls(j-lo(2)) = -dhy*bY(i,j+1,k)
                  r_ls(j-lo(2)) = rhs(i,j,k,n) + rho_xz - phi(i,j,k,n)*delta

                  if (j .eq. lo(2))
     $               r_ls(j-lo(2)) = r_ls(j-lo(2)) + dhy*bY(i,j,k)*phi(i,j-1,k,n)

                  if (j .eq. hi(2))
     $               r_ls(j-lo(2)) = r_ls(j-lo(2)) + dhy*bY(i,j+1,k)*phi(i,j+1,k,n)
                end do
c
                call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,jlen)

                do j = lo(2),hi(2)
                  phi(i,j,k,n) = u_ls(j-lo(2))
                end do
c                  
               end do
            end do

         else if (do_line .eq. 1) then

           do k = lo(3), hi(3)
             ioff = MOD(k + redblack,2)
             do j = lo(2) + ioff,hi(2),2
c
               do i = lo(1),hi(1)
                  cf0 = cvmgt(f0(lo(1),j,k), 0.0,
     $                 (i .eq. lo(1)) .and. (m0(lo(1)-1,j,k).gt.0))
                  cf1 = cvmgt(f1(i,lo(2),k), 0.0,
     $                 (j .eq. lo(2)) .and. (m1(i,lo(2)-1,k).gt.0))
                  cf2 = cvmgt(f2(i,j,lo(3)), 0.0,
     $                 (k .eq. lo(3)) .and. (m2(i,j,lo(3)-1).gt.0))
                  cf3 = cvmgt(f3(hi(1),j,k), 0.0,
     $                 (i .eq. hi(1)) .and. (m3(hi(1)+1,j,k).gt.0))
                  cf4 = cvmgt(f4(i,hi(2),k), 0.0,
     $                 (j .eq. hi(2)) .and. (m4(i,hi(2)+1,k).gt.0))
                  cf5 = cvmgt(f5(i,j,hi(3)), 0.0,
     $                 (k .eq. hi(3)) .and. (m5(i,j,hi(3)+1).gt.0))
c
                  delta = dhx*(bX(i,j,k)*cf0 + bX(i+1,j,k)*cf3)
     $                 +  dhy*(bY(i,j,k)*cf1 + bY(i,j+1,k)*cf4)
     $                 +  dhz*(bZ(i,j,k)*cf2 + bZ(i,j,k+1)*cf5)
c                  
                  gamma = alpha*a(i,j,k)
     $                 +   dhx*(bX(i,j,k)+bX(i+1,j,k))
     $                 +   dhy*(bY(i,j,k)+bY(i,j+1,k))
     $                 +   dhz*(bZ(i,j,k)+bZ(i,j,k+1))
c
                  rho_yz =  dhy*( bY(i,j  ,k)*phi(i,j-1,k,n)
     $                    +       bY(i,j+1,k)*phi(i,j+1,k,n) )
     $                    + dhz*( bZ(i,j,k  )*phi(i,j,k-1,n)
     $                    +       bZ(i,j,k+1)*phi(i,j,k+1,n) )

                  a_ls(i-lo(1)) = -dhx*bX(i,j,k)
                  b_ls(i-lo(1)) = gamma - delta
                  c_ls(i-lo(1)) = -dhx*bX(i+1,j,k)
                  r_ls(i-lo(1)) = rhs(i,j,k,n) + rho_yz - phi(i,j,k,n)*delta

                  if (i .eq. lo(1))
     $               r_ls(i-lo(1)) = r_ls(i-lo(1)) + dhx*bX(i,j,k)*phi(i-1,j,k,n)

                  if (i .eq. hi(1))
     $               r_ls(i-lo(1)) = r_ls(i-lo(1)) + dhx*bX(i+1,j,k)*phi(i+1,j,k,n)
                end do
c
                call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,ilen)

                do i = lo(1),hi(1)
                  phi(i,j,k,n) = u_ls(i-lo(1))
                end do
c                  
               end do
            end do
         endif
      end do
c     
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
     $     bZ,DIMS(bZ),
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
      integer DIMDEC(bZ)
      REAL_T zz(DIMV(zz),nc)
      REAL_T rr(DIMV(rr),nc)
      REAL_T  a(DIMV(a))
      REAL_T bX(DIMV(bX))
      REAL_T bY(DIMV(bY))
      REAL_T bZ(DIMV(bZ))
      REAL_T h(BL_SPACEDIM), rho
c
      integer i,j, k, n
      REAL_T de, dhx, dhy, dhz
c
      rho = 0.0
      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      dhz = beta/h(3)**2
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  de = 1.0/(alpha*a(i,j,k)
     $                 + dhx*(bX(i+1,j,k)+bX(i,j,k))
     $                 + dhy*(bY(i,j+1,k)+bY(i,j,k))
     $                 + dhz*(bZ(i,j,k+1)+bZ(i,j,k))
     $                 )
                  zz(i,j,k,n) = rr(i,j,k,n)*de
                  rho = rho + zz(i,j,k,n)*rr(i,j,k,n)
               end do
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
     $     bZ,DIMS(bZ),
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
      integer DIMDEC(bZ)
      REAL_T  y(DIMV(y),nc)
      REAL_T  x(DIMV(x),nc)
      REAL_T  a(DIMV(a))
      REAL_T bX(DIMV(bX))
      REAL_T bY(DIMV(bY))
      REAL_T bZ(DIMV(bZ))
      REAL_T h(BL_SPACEDIM)
c
      integer i,j,k,n
      REAL_T dhx,dhy,dhz
c
      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      dhz = beta/h(3)**2
c
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n)
     $                 - dhx*
     $                 (   bX(i+1,j,k)*( x(i+1,j,k,n) - x(i  ,j,k,n) )
     $                 -   bX(i  ,j,k)*( x(i  ,j,k,n) - x(i-1,j,k,n) ) )
     $                 - dhy*
     $                 (   bY(i,j+1,k)*( x(i,j+1,k,n) - x(i,j  ,k,n) )
     $                 -   bY(i,j  ,k)*( x(i,j  ,k,n) - x(i,j-1,k,n) ) )
     $                 - dhz*
     $                 (   bZ(i,j,k+1)*( x(i,j,k+1,n) - x(i,j,k  ,n) )
     $                 -   bZ(i,j,k  )*( x(i,j,k  ,n) - x(i,j,k-1,n) ) )
               end do
            end do
         end do
      end do
      end

c-----------------------------------------------------------------------
c
c     Fill in a matrix x vector operator here
c
      subroutine FORT_NORMA(
     &     res,
     $     alpha, beta,
     $     a, DIMS(a),
     $     bX,DIMS(bX),
     $     bY,DIMS(bY),
     $     bZ,DIMS(bZ),
     $     lo,hi,nc,
     $     h
     $     )
      REAL_T alpha, beta, res
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer DIMDEC(a)
      integer DIMDEC(bX)
      integer DIMDEC(bY)
      integer DIMDEC(bZ)
      REAL_T  a(DIMV(a))
      REAL_T bX(DIMV(bX))
      REAL_T bY(DIMV(bY))
      REAL_T bZ(DIMV(bZ))
      REAL_T h(BL_SPACEDIM)
c
      integer i,j,k,n
      REAL_T dhx,dhy,dhz
c
      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      dhz = beta/h(3)**2
c
      res = 0.0
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  res = max(res, 
     &                 abs(alpha*a(i,j,k)
     &                 + dhx*(bX(i+1,j,k) + bX(i,j,k))
     &                 + dhy*(bY(i,j+1,k) + bY(i,j,k))
     $                 - dhz*(bZ(i,j,k+1) + bZ(i,j,k)))
     &                 + abs( -dhx*bX(i+1,j,k)) + abs( -dhx*bX(i,j,k))
     &                 + abs( -dhy*bY(i,j+1,k)) + abs( -dhy*bY(i,j,k))
     &                 + abs( -dhz*bZ(i,j,k+1)) + abs( -dhz*bZ(i,j,k)))
               end do
            end do
         end do
      end do
      end
c-----------------------------------------------------------------------
c
c     Fill in fluxes
c
      subroutine FORT_FLUX(
     $     x,DIMS(x),
     $     alpha, beta,
     $     a, DIMS(a),
     $     bX,DIMS(bX),
     $     bY,DIMS(bY),
     $     bZ,DIMS(bZ),
     $     lo,hi,nc,
     $     h,
     $     xflux,DIMS(xflux),
     $     yflux,DIMS(yflux),
     $     zflux,DIMS(zflux)
     $     )
      implicit none
      REAL_T alpha, beta
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM), nc
      integer DIMDEC(x)
      integer DIMDEC(a)
      integer DIMDEC(bX)
      integer DIMDEC(bY)
      integer DIMDEC(bZ)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
      REAL_T  x(DIMV(x),nc)
      REAL_T  a(DIMV(a))
      REAL_T bX(DIMV(bX))
      REAL_T bY(DIMV(bY))
      REAL_T bZ(DIMV(bZ))
      REAL_T xflux(DIMV(xflux),nc)
      REAL_T yflux(DIMV(yflux),nc)
      REAL_T zflux(DIMV(zflux),nc)
      REAL_T h(BL_SPACEDIM)
c
      REAL_T dhx, dhy, dhz
      integer i,j,k,n
c
      dhx = one/h(1)
      dhy = one/h(2)
      dhz = one/h(3)
c
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)+1
                  xflux(i,j,k,n) = - dhx*bX(i,j,k)*( x(i,j,k,n) - x(i-1,j,k,n) )
               end do
            end do
         end do
      end do
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)+1
               do i = lo(1), hi(1)
                  yflux(i,j,k,n) = - dhy*bY(i,j,k)*( x(i,j,k,n) - x(i,j-1,k,n) )
               end do
            end do
         end do
      end do
      do n = 1, nc
         do k = lo(3), hi(3)+1
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  zflux(i,j,k,n) = - dhz*bZ(i,j,k)*( x(i,j,k,n) - x(i,j,k-1,n) )
               end do
            end do
         end do
      end do
      end

      
