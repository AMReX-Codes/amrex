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
!     L=(gamma-delta)*y - rho + delta*yOld.  The resulting iteration is
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
           phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
           rhs,rhs_l1,rhs_l2,rhs_l3,rhs_h1,rhs_h2,rhs_h3, &
           alpha, beta, &
           a,  a_l1,a_l2,a_l3,a_h1,a_h2,a_h3, &
           bX, bX_l1,bX_l2,bX_l3,bX_h1,bX_h2,bX_h3, &
           bY, bY_l1,bY_l2,bY_l3,bY_h1,bY_h2,bY_h3, &
           bZ, bZ_l1,bZ_l2,bZ_l3,bZ_h1,bZ_h2,bZ_h3, &
           f0, f0_l1,f0_l2,f0_l3,f0_h1,f0_h2,f0_h3, &
           m0, m0_l1,m0_l2,m0_l3,m0_h1,m0_h2,m0_h3, &
           f1, f1_l1,f1_l2,f1_l3,f1_h1,f1_h2,f1_h3, &
           m1, m1_l1,m1_l2,m1_l3,m1_h1,m1_h2,m1_h3, &
           f2, f2_l1,f2_l2,f2_l3,f2_h1,f2_h2,f2_h3, &
           m2, m2_l1,m2_l2,m2_l3,m2_h1,m2_h2,m2_h3, &
           f3, f3_l1,f3_l2,f3_l3,f3_h1,f3_h2,f3_h3, &
           m3, m3_l1,m3_l2,m3_l3,m3_h1,m3_h2,m3_h3, &
           f4, f4_l1,f4_l2,f4_l3,f4_h1,f4_h2,f4_h3, &
           m4, m4_l1,m4_l2,m4_l3,m4_h1,m4_h2,m4_h3, &
           f5, f5_l1,f5_l2,f5_l3,f5_h1,f5_h2,f5_h3, &
           m5, m5_l1,m5_l2,m5_l3,m5_h1,m5_h2,m5_h3, &
           lo,hi,blo,bhi, &
           nc, h,redblack &
           )

      implicit none

      REAL_T alpha, beta
      integer DIMDEC(phi)
      integer DIMDEC(rhs)
      integer DIMDEC(a)
      integer DIMDEC(bX)
      integer DIMDEC(bY)
      integer DIMDEC(bZ)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
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

      integer  i, j, k, ioff, n

      REAL_T dhx, dhy, dhz, cf0, cf1, cf2, cf3, cf4, cf5
      REAL_T g_m_d, gamma, rho, res

!     This factor of 1.15 in 3D does over-relaxation but seems to consistently reduce the number of V-cycles needed.
      REAL_T omega
      omega = 1.15d0

      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      dhz = beta/h(3)**2

      do n = 1, nc
          do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               ioff = MOD(lo(1) + j + k + redblack,2)
               do i = lo(1) + ioff,hi(1),2

                  cf0 = merge(f0(blo(1),j,k), 0.0D0, &
                       (i .eq. blo(1)) .and. (m0(blo(1)-1,j,k).gt.0))
                  cf1 = merge(f1(i,blo(2),k), 0.D00, &
                       (j .eq. blo(2)) .and. (m1(i,blo(2)-1,k).gt.0))
                  cf2 = merge(f2(i,j,blo(3)), 0.0D0, &
                       (k .eq. blo(3)) .and. (m2(i,j,blo(3)-1).gt.0))
                  cf3 = merge(f3(bhi(1),j,k), 0.0D0, &
                       (i .eq. bhi(1)) .and. (m3(bhi(1)+1,j,k).gt.0))
                  cf4 = merge(f4(i,bhi(2),k), 0.0D0, &
                       (j .eq. bhi(2)) .and. (m4(i,bhi(2)+1,k).gt.0))
                  cf5 = merge(f5(i,j,bhi(3)), 0.0D0, &
                       (k .eq. bhi(3)) .and. (m5(i,j,bhi(3)+1).gt.0))

                  gamma = alpha*a(i,j,k) &
                       +   dhx*(bX(i,j,k)+bX(i+1,j,k)) &
                       +   dhy*(bY(i,j,k)+bY(i,j+1,k)) &
                       +   dhz*(bZ(i,j,k)+bZ(i,j,k+1))

                  g_m_d = gamma &
                       - (dhx*(bX(i,j,k)*cf0 + bX(i+1,j,k)*cf3) &
                       +  dhy*(bY(i,j,k)*cf1 + bY(i,j+1,k)*cf4) &
                       +  dhz*(bZ(i,j,k)*cf2 + bZ(i,j,k+1)*cf5))

                  rho =  dhx*( bX(i  ,j,k)*phi(i-1,j,k,n) &
                       +       bX(i+1,j,k)*phi(i+1,j,k,n) ) &
                       + dhy*( bY(i,j  ,k)*phi(i,j-1,k,n) &
                       +       bY(i,j+1,k)*phi(i,j+1,k,n) ) &
                       + dhz*( bZ(i,j,k  )*phi(i,j,k-1,n) &
                       +       bZ(i,j,k+1)*phi(i,j,k+1,n) )

                  res =  rhs(i,j,k,n) - (gamma*phi(i,j,k,n) - rho)
                  phi(i,j,k,n) = phi(i,j,k,n) + omega/g_m_d * res

               end do
            end do
          end do
      end do

    end subroutine FORT_GSRB

!-----------------------------------------------------------------------
!      
!     Jacobi:
!     Apply the Jacobi relaxation to the state phi for the equation
!     L(phi) = alpha*a(x)*phi(x) - beta*Div(b(x)Grad(phi(x))) = rhs(x)
!     central differenced, according to the arrays of boundary
!     masks (m#) and auxiliary data (f#).
!     
!     In general, if the linear operator L=gamma*y-rho, the GS relaxation
!     is y = (R - rho)/gamma.  Near a boundary, the ghost data is filled
!     using a polynomial interpolant based on the "old" phi values, so
!     L=(gamma-delta)*y - rho + delta*yOld.  The resulting iteration is
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
           phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
           rhs,rhs_l1,rhs_l2,rhs_l3,rhs_h1,rhs_h2,rhs_h3, &
           alpha, beta, &
           a,  a_l1,a_l2,a_l3,a_h1,a_h2,a_h3, &
           bX, bX_l1,bX_l2,bX_l3,bX_h1,bX_h2,bX_h3, &
           bY, bY_l1,bY_l2,bY_l3,bY_h1,bY_h2,bY_h3, &
           bZ, bZ_l1,bZ_l2,bZ_l3,bZ_h1,bZ_h2,bZ_h3, &
           f0, f0_l1,f0_l2,f0_l3,f0_h1,f0_h2,f0_h3, &
           m0, m0_l1,m0_l2,m0_l3,m0_h1,m0_h2,m0_h3, &
           f1, f1_l1,f1_l2,f1_l3,f1_h1,f1_h2,f1_h3, &
           m1, m1_l1,m1_l2,m1_l3,m1_h1,m1_h2,m1_h3, &
           f2, f2_l1,f2_l2,f2_l3,f2_h1,f2_h2,f2_h3, &
           m2, m2_l1,m2_l2,m2_l3,m2_h1,m2_h2,m2_h3, &
           f3, f3_l1,f3_l2,f3_l3,f3_h1,f3_h2,f3_h3, &
           m3, m3_l1,m3_l2,m3_l3,m3_h1,m3_h2,m3_h3, &
           f4, f4_l1,f4_l2,f4_l3,f4_h1,f4_h2,f4_h3, &
           m4, m4_l1,m4_l2,m4_l3,m4_h1,m4_h2,m4_h3, &
           f5, f5_l1,f5_l2,f5_l3,f5_h1,f5_h2,f5_h3, &
           m5, m5_l1,m5_l2,m5_l3,m5_h1,m5_h2,m5_h3, &
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
      integer DIMDEC(bZ)
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

      integer  i, j, k, n

      REAL_T dhx, dhy, dhz, cf0, cf1, cf2, cf3, cf4, cf5
      REAL_T delta, gamma, rho

      REAL_T, allocatable :: phinew(:,:,:)

      allocate(phinew(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      dhz = beta/h(3)**2

      do n = 1, nc
          do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1),hi(1)

                  cf0 = merge(f0(lo(1),j,k), 0.0D0, &
                       (i .eq. lo(1)) .and. (m0(lo(1)-1,j,k).gt.0))
                  cf1 = merge(f1(i,lo(2),k), 0.D00, &
                       (j .eq. lo(2)) .and. (m1(i,lo(2)-1,k).gt.0))
                  cf2 = merge(f2(i,j,lo(3)), 0.0D0, &
                       (k .eq. lo(3)) .and. (m2(i,j,lo(3)-1).gt.0))
                  cf3 = merge(f3(hi(1),j,k), 0.0D0, &
                       (i .eq. hi(1)) .and. (m3(hi(1)+1,j,k).gt.0))
                  cf4 = merge(f4(i,hi(2),k), 0.0D0, &
                       (j .eq. hi(2)) .and. (m4(i,hi(2)+1,k).gt.0))
                  cf5 = merge(f5(i,j,hi(3)), 0.0D0, &
                       (k .eq. hi(3)) .and. (m5(i,j,hi(3)+1).gt.0))

                  delta = dhx*(bX(i,j,k)*cf0 + bX(i+1,j,k)*cf3) &
                       +  dhy*(bY(i,j,k)*cf1 + bY(i,j+1,k)*cf4) &
                       +  dhz*(bZ(i,j,k)*cf2 + bZ(i,j,k+1)*cf5)

                  gamma = alpha*a(i,j,k) &
                       +   dhx*(bX(i,j,k)+bX(i+1,j,k)) &
                       +   dhy*(bY(i,j,k)+bY(i,j+1,k)) &
                       +   dhz*(bZ(i,j,k)+bZ(i,j,k+1))

                  rho =  dhx*( bX(i  ,j,k)*phi(i-1,j,k,n) &
                       +       bX(i+1,j,k)*phi(i+1,j,k,n) ) &
                       + dhy*( bY(i,j  ,k)*phi(i,j-1,k,n) &
                       +       bY(i,j+1,k)*phi(i,j+1,k,n) ) &
                       + dhz*( bZ(i,j,k  )*phi(i,j,k-1,n) &
                       +       bZ(i,j,k+1)*phi(i,j,k+1,n) )

                  phinew(i,j,k) = (rhs(i,j,k,n)+rho-phi(i,j,k,n)*delta) &
                       /                   (gamma - delta)

               end do
            end do
          end do

         phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),n) = &
               phinew(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

      end do

      deallocate(phinew)

    end subroutine FORT_JACOBI

!-----------------------------------------------------------------------
!
!     Fill in a matrix x vector operator here
!
    subroutine FORT_ADOTX( &
           y,y_l1,y_l2,y_l3,y_h1,y_h2,y_h3, &
           x,x_l1,x_l2,x_l3,x_h1,x_h2,x_h3, &
           alpha, beta, &
           a, a_l1,a_l2,a_l3,a_h1,a_h2,a_h3, &
           bX,bX_l1,bX_l2,bX_l3,bX_h1,bX_h2,bX_h3, &
           bY,bY_l1,bY_l2,bY_l3,bY_h1,bY_h2,bY_h3, &
           bZ,bZ_l1,bZ_l2,bZ_l3,bZ_h1,bZ_h2,bZ_h3, &
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
      integer DIMDEC(bZ)
      REAL_T  y(DIMV(y),nc)
      REAL_T  x(DIMV(x),nc)
      REAL_T  a(DIMV(a))
      REAL_T bX(DIMV(bX))
      REAL_T bY(DIMV(bY))
      REAL_T bZ(DIMV(bZ))
      REAL_T h(BL_SPACEDIM)

      integer i,j,k,n
      REAL_T dhx,dhy,dhz

      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      dhz = beta/h(3)**2

      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n) &
                       - dhx* &
                       (   bX(i+1,j,k)*( x(i+1,j,k,n) - x(i  ,j,k,n) ) &
                       -   bX(i  ,j,k)*( x(i  ,j,k,n) - x(i-1,j,k,n) ) ) &
                       - dhy* &
                       (   bY(i,j+1,k)*( x(i,j+1,k,n) - x(i,j  ,k,n) ) &
                       -   bY(i,j  ,k)*( x(i,j  ,k,n) - x(i,j-1,k,n) ) ) &
                       - dhz* &
                       (   bZ(i,j,k+1)*( x(i,j,k+1,n) - x(i,j,k  ,n) ) &
                       -   bZ(i,j,k  )*( x(i,j,k  ,n) - x(i,j,k-1,n) ) )
               end do
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
           a, a_l1,a_l2,a_l3,a_h1,a_h2,a_h3, &
           bX,bX_l1,bX_l2,bX_l3,bX_h1,bX_h2,bX_h3, &
           bY,bY_l1,bY_l2,bY_l3,bY_h1,bY_h2,bY_h3, &
           bZ,bZ_l1,bZ_l2,bZ_l3,bZ_h1,bZ_h2,bZ_h3, &
           lo,hi,nc, &
           h &
           )

      implicit none

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

      integer i,j,k,n
      REAL_T dhx,dhy,dhz

      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      dhz = beta/h(3)**2

      res = 0.0D0

      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  res = max(res, abs(alpha*a(i,j,k) &
                       + dhx*(bX(i+1,j,k) + bX(i,j,k)) &
                       + dhy*(bY(i,j+1,k) + bY(i,j,k)) &
                       + dhz*(bZ(i,j,k+1) + bZ(i,j,k))) &
                       + abs( -dhx*bX(i+1,j,k)) + abs( -dhx*bX(i,j,k)) &
                       + abs( -dhy*bY(i,j+1,k)) + abs( -dhy*bY(i,j,k)) &
                       + abs( -dhz*bZ(i,j,k+1)) + abs( -dhz*bZ(i,j,k)))
               end do
            end do
         end do
      end do

    end subroutine FORT_NORMA

!-----------------------------------------------------------------------
!
!     Fill in fluxes
!
    subroutine FORT_FLUX( &
           x,x_l1,x_l2,x_l3,x_h1,x_h2,x_h3, &
           alpha, beta, &
           a, a_l1,a_l2,a_l3,a_h1,a_h2,a_h3, &
           bX,bX_l1,bX_l2,bX_l3,bX_h1,bX_h2,bX_h3, &
           bY,bY_l1,bY_l2,bY_l3,bY_h1,bY_h2,bY_h3, &
           bZ,bZ_l1,bZ_l2,bZ_l3,bZ_h1,bZ_h2,bZ_h3, &
           xlo,xhi, &
           ylo,yhi, &
           zlo,zhi, &
           nc, &
           h, &
           xflux,xflux_l1,xflux_l2,xflux_l3,xflux_h1,xflux_h2,xflux_h3, &
           yflux,yflux_l1,yflux_l2,yflux_l3,yflux_h1,yflux_h2,yflux_h3, &
           zflux,zflux_l1,zflux_l2,zflux_l3,zflux_h1,zflux_h2,zflux_h3 &
           )

      implicit none

      REAL_T alpha, beta
      integer xlo(BL_SPACEDIM), xhi(BL_SPACEDIM)
      integer ylo(BL_SPACEDIM), yhi(BL_SPACEDIM)
      integer zlo(BL_SPACEDIM), zhi(BL_SPACEDIM)
      integer nc
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

      REAL_T dhx, dhy, dhz
      integer i,j,k,n

      dhx = one/h(1)
      dhy = one/h(2)
      dhz = one/h(3)

      do n = 1, nc
         do k = xlo(3), xhi(3)
            do j = xlo(2), xhi(2)
               do i = xlo(1), xhi(1)
                  xflux(i,j,k,n) = - dhx*bX(i,j,k)*( x(i,j,k,n) - x(i-1,j,k,n) )
               end do
            end do
         end do

         do k = ylo(3), yhi(3)
            do j = ylo(2), yhi(2)
               do i = ylo(1), yhi(1)
                  yflux(i,j,k,n) = - dhy*bY(i,j,k)*( x(i,j,k,n) - x(i,j-1,k,n) )
               end do
            end do
         end do

         do k = zlo(3), zhi(3)
            do j = zlo(2), zhi(2)
               do i = zlo(1), zhi(1)
                  zflux(i,j,k,n) = - dhz*bZ(i,j,k)*( x(i,j,k,n) - x(i,j,k-1,n) )
               end do
            end do
         end do
      end do

    end subroutine FORT_FLUX
  
