c
c $Id: ABec_3D.F,v 1.2 1998-12-14 20:23:47 lijewski Exp $
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
      REAL_T delta, gamma, rho
c
      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      dhz = beta/h(3)**2
      do n = 1, nc
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
