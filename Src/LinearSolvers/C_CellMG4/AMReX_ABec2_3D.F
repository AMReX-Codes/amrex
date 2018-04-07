#include <AMReX_REAL.H>
#include "AMReX_ArrayLim.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_LO_BCTYPES.H"

      subroutine ab2_gsrb (lo, hi, blo, bhi,
     &     phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3,
     &     res,res_l1,res_l2,res_l3,res_h1,res_h2,res_h3,
     &     alpha, beta,
     &     a, a_l1, a_l2, a_l3, a_h1, a_h2, a_h3, 
     &     bX,bX_l1,bX_l2,bX_l3,bX_h1,bX_h2,bX_h3,
     &     bY,bY_l1,bY_l2,bY_l3,bY_h1,bY_h2,bY_h3,
     &     bZ,bZ_l1,bZ_l2,bZ_l3,bZ_h1,bZ_h2,bZ_h3,
     &     f0,f0_l1,f0_l2,f0_l3,f0_h1,f0_h2,f0_h3,
     &     m0,m0_l1,m0_l2,m0_l3,m0_h1,m0_h2,m0_h3,
     &     f1,f1_l1,f1_l2,f1_l3,f1_h1,f1_h2,f1_h3,
     &     m1,m1_l1,m1_l2,m1_l3,m1_h1,m1_h2,m1_h3,
     &     f2,f2_l1,f2_l2,f2_l3,f2_h1,f2_h2,f2_h3,
     &     m2,m2_l1,m2_l2,m2_l3,m2_h1,m2_h2,m2_h3,
     &     f3,f3_l1,f3_l2,f3_l3,f3_h1,f3_h2,f3_h3,
     &     m3,m3_l1,m3_l2,m3_l3,m3_h1,m3_h2,m3_h3,
     &     f4,f4_l1,f4_l2,f4_l3,f4_h1,f4_h2,f4_h3,
     &     m4,m4_l1,m4_l2,m4_l3,m4_h1,m4_h2,m4_h3,
     &     f5,f5_l1,f5_l2,f5_l3,f5_h1,f5_h2,f5_h3,
     &     m5,m5_l1,m5_l2,m5_l3,m5_h1,m5_h2,m5_h3,
     &     nc, h, redBlack)

      implicit none
      integer lo(3), hi(3), blo(3), bhi(3), nc, redblack

      integer phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      integer res_l1,res_l2,res_l3,res_h1,res_h2,res_h3
      REAL_T alpha, beta
      integer a_l1,  a_l2, a_l3, a_h1, a_h2, a_h3
      integer bX_l1,bX_l2,bX_l3,bX_h1,bX_h2,bX_h3
      integer bY_l1,bY_l2,bY_l3,bY_h1,bY_h2,bY_h3
      integer bZ_l1,bZ_l2,bZ_l3,bZ_h1,bZ_h2,bZ_h3
      integer f0_l1,f0_l2,f0_l3,f0_h1,f0_h2,f0_h3
      integer m0_l1,m0_l2,m0_l3,m0_h1,m0_h2,m0_h3
      integer f1_l1,f1_l2,f1_l3,f1_h1,f1_h2,f1_h3
      integer m1_l1,m1_l2,m1_l3,m1_h1,m1_h2,m1_h3
      integer f2_l1,f2_l2,f2_l3,f2_h1,f2_h2,f2_h3
      integer m2_l1,m2_l2,m2_l3,m2_h1,m2_h2,m2_h3
      integer f3_l1,f3_l2,f3_l3,f3_h1,f3_h2,f3_h3
      integer m3_l1,m3_l2,m3_l3,m3_h1,m3_h2,m3_h3
      integer f4_l1,f4_l2,f4_l3,f4_h1,f4_h2,f4_h3
      integer m4_l1,m4_l2,m4_l3,m4_h1,m4_h2,m4_h3
      integer f5_l1,f5_l2,f5_l3,f5_h1,f5_h2,f5_h3
      integer m5_l1,m5_l2,m5_l3,m5_h1,m5_h2,m5_h3

      REAL_T  phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3,1:nc)
      REAL_T  res(res_l1:res_h1,res_l2:res_h2,res_l3:res_h3,1:nc)
      REAL_T   a( a_l1: a_h1, a_l2: a_h2, a_l3: a_h3)
      REAL_T  bX(bX_l1:bX_h1,bX_l2:bX_h2,bX_l3:bX_h3)
      REAL_T  bY(bY_l1:bY_h1,bY_l2:bY_h2,bY_l3:bY_h3)
      REAL_T  bZ(bZ_l1:bZ_h1,bZ_l2:bZ_h2,bZ_l3:bZ_h3)
      REAL_T  f0(f0_l1:f0_h1,f0_l2:f0_h2,f0_l3:f0_h3)
      REAL_T  f1(f1_l1:f1_h1,f1_l2:f1_h2,f1_l3:f1_h3)
      REAL_T  f2(f2_l1:f2_h1,f2_l2:f2_h2,f2_l3:f2_h3)
      REAL_T  f3(f3_l1:f3_h1,f3_l2:f3_h2,f3_l3:f3_h3)
      REAL_T  f4(f4_l1:f4_h1,f4_l2:f4_h2,f4_l3:f4_h3)
      REAL_T  f5(f5_l1:f5_h1,f5_l2:f5_h2,f5_l3:f5_h3)
      integer  m0(m0_l1:m0_h1,m0_l2:m0_h2,m0_l3:m0_h3)
      integer  m1(m1_l1:m1_h1,m1_l2:m1_h2,m1_l3:m1_h3)
      integer  m2(m2_l1:m2_h1,m2_l2:m2_h2,m2_l3:m2_h3)
      integer  m3(m3_l1:m3_h1,m3_l2:m3_h2,m3_l3:m3_h3)
      integer  m4(m4_l1:m4_h1,m4_l2:m4_h2,m4_l3:m4_h3)
      integer  m5(m5_l1:m5_h1,m5_l2:m5_h2,m5_l3:m5_h3)
      REAL_T  h(3)
c
      integer  i, j, k, ioff, n
c
      REAL_T dhx, dhy, dhz, cf0, cf1, cf2, cf3, cf4, cf5
      REAL_T delta, gamma
c
      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      dhz = beta/h(3)**2

      do n = 1, nc
          do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               ioff = MOD(j + k + redblack,2)
               do i = lo(1) + ioff,hi(1),2

                  cf0 = merge(f0(blo(1),j,k), 0.0D0,
     $                 (i .eq. blo(1)) .and. (m0(blo(1)-1,j,k).gt.0))
                  cf1 = merge(f1(i,blo(2),k), 0.D00,
     $                 (j .eq. blo(2)) .and. (m1(i,blo(2)-1,k).gt.0))
                  cf2 = merge(f2(i,j,blo(3)), 0.0D0,
     $                 (k .eq. blo(3)) .and. (m2(i,j,blo(3)-1).gt.0))
                  cf3 = merge(f3(bhi(1),j,k), 0.0D0,
     $                 (i .eq. bhi(1)) .and. (m3(bhi(1)+1,j,k).gt.0))
                  cf4 = merge(f4(i,bhi(2),k), 0.0D0,
     $                 (j .eq. bhi(2)) .and. (m4(i,bhi(2)+1,k).gt.0))
                  cf5 = merge(f5(i,j,bhi(3)), 0.0D0,
     $                 (k .eq. bhi(3)) .and. (m5(i,j,bhi(3)+1).gt.0))

                  delta = dhx*(bX(i,j,k)*cf0 + bX(i+1,j,k)*cf3)
     $                 +  dhy*(bY(i,j,k)*cf1 + bY(i,j+1,k)*cf4)
     $                 +  dhz*(bZ(i,j,k)*cf2 + bZ(i,j,k+1)*cf5)

                  gamma = alpha*a(i,j,k)
     $                 +   dhx*(bX(i,j,k)+bX(i+1,j,k))
     $                 +   dhy*(bY(i,j,k)+bY(i,j+1,k))
     $                 +   dhz*(bZ(i,j,k)+bZ(i,j,k+1))

                  phi(i,j,k,n) = phi(i,j,k,n) + res(i,j,k,n)/(gamma - delta)

               end do
            end do
          end do
      end do
      end


      subroutine ab2_jacobi (lo, hi,
     &     phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3,
     &     res,res_l1,res_l2,res_l3,res_h1,res_h2,res_h3,
     &     alpha, beta,
     &     a, a_l1, a_l2, a_l3, a_h1, a_h2, a_h3, 
     &     bX,bX_l1,bX_l2,bX_l3,bX_h1,bX_h2,bX_h3,
     &     bY,bY_l1,bY_l2,bY_l3,bY_h1,bY_h2,bY_h3,
     &     bZ,bZ_l1,bZ_l2,bZ_l3,bZ_h1,bZ_h2,bZ_h3,
     &     f0,f0_l1,f0_l2,f0_l3,f0_h1,f0_h2,f0_h3,
     &     m0,m0_l1,m0_l2,m0_l3,m0_h1,m0_h2,m0_h3,
     &     f1,f1_l1,f1_l2,f1_l3,f1_h1,f1_h2,f1_h3,
     &     m1,m1_l1,m1_l2,m1_l3,m1_h1,m1_h2,m1_h3,
     &     f2,f2_l1,f2_l2,f2_l3,f2_h1,f2_h2,f2_h3,
     &     m2,m2_l1,m2_l2,m2_l3,m2_h1,m2_h2,m2_h3,
     &     f3,f3_l1,f3_l2,f3_l3,f3_h1,f3_h2,f3_h3,
     &     m3,m3_l1,m3_l2,m3_l3,m3_h1,m3_h2,m3_h3,
     &     f4,f4_l1,f4_l2,f4_l3,f4_h1,f4_h2,f4_h3,
     &     m4,m4_l1,m4_l2,m4_l3,m4_h1,m4_h2,m4_h3,
     &     f5,f5_l1,f5_l2,f5_l3,f5_h1,f5_h2,f5_h3,
     &     m5,m5_l1,m5_l2,m5_l3,m5_h1,m5_h2,m5_h3,
     &     nc, h)

      implicit none
      integer lo(3), hi(3), nc

      integer phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
      integer res_l1,res_l2,res_l3,res_h1,res_h2,res_h3
      REAL_T alpha, beta
      integer a_l1,  a_l2, a_l3, a_h1, a_h2, a_h3
      integer bX_l1,bX_l2,bX_l3,bX_h1,bX_h2,bX_h3
      integer bY_l1,bY_l2,bY_l3,bY_h1,bY_h2,bY_h3
      integer bZ_l1,bZ_l2,bZ_l3,bZ_h1,bZ_h2,bZ_h3
      integer f0_l1,f0_l2,f0_l3,f0_h1,f0_h2,f0_h3
      integer m0_l1,m0_l2,m0_l3,m0_h1,m0_h2,m0_h3
      integer f1_l1,f1_l2,f1_l3,f1_h1,f1_h2,f1_h3
      integer m1_l1,m1_l2,m1_l3,m1_h1,m1_h2,m1_h3
      integer f2_l1,f2_l2,f2_l3,f2_h1,f2_h2,f2_h3
      integer m2_l1,m2_l2,m2_l3,m2_h1,m2_h2,m2_h3
      integer f3_l1,f3_l2,f3_l3,f3_h1,f3_h2,f3_h3
      integer m3_l1,m3_l2,m3_l3,m3_h1,m3_h2,m3_h3
      integer f4_l1,f4_l2,f4_l3,f4_h1,f4_h2,f4_h3
      integer m4_l1,m4_l2,m4_l3,m4_h1,m4_h2,m4_h3
      integer f5_l1,f5_l2,f5_l3,f5_h1,f5_h2,f5_h3
      integer m5_l1,m5_l2,m5_l3,m5_h1,m5_h2,m5_h3

      REAL_T  phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3,1:nc)
      REAL_T  res(res_l1:res_h1,res_l2:res_h2,res_l3:res_h3,1:nc)
      REAL_T   a( a_l1: a_h1, a_l2: a_h2, a_l3: a_h3)
      REAL_T  bX(bX_l1:bX_h1,bX_l2:bX_h2,bX_l3:bX_h3)
      REAL_T  bY(bY_l1:bY_h1,bY_l2:bY_h2,bY_l3:bY_h3)
      REAL_T  bZ(bZ_l1:bZ_h1,bZ_l2:bZ_h2,bZ_l3:bZ_h3)
      REAL_T  f0(f0_l1:f0_h1,f0_l2:f0_h2,f0_l3:f0_h3)
      REAL_T  f1(f1_l1:f1_h1,f1_l2:f1_h2,f1_l3:f1_h3)
      REAL_T  f2(f2_l1:f2_h1,f2_l2:f2_h2,f2_l3:f2_h3)
      REAL_T  f3(f3_l1:f3_h1,f3_l2:f3_h2,f3_l3:f3_h3)
      REAL_T  f4(f4_l1:f4_h1,f4_l2:f4_h2,f4_l3:f4_h3)
      REAL_T  f5(f5_l1:f5_h1,f5_l2:f5_h2,f5_l3:f5_h3)
      integer  m0(m0_l1:m0_h1,m0_l2:m0_h2,m0_l3:m0_h3)
      integer  m1(m1_l1:m1_h1,m1_l2:m1_h2,m1_l3:m1_h3)
      integer  m2(m2_l1:m2_h1,m2_l2:m2_h2,m2_l3:m2_h3)
      integer  m3(m3_l1:m3_h1,m3_l2:m3_h2,m3_l3:m3_h3)
      integer  m4(m4_l1:m4_h1,m4_l2:m4_h2,m4_l3:m4_h3)
      integer  m5(m5_l1:m5_h1,m5_l2:m5_h2,m5_l3:m5_h3)
      REAL_T  h(3)
c
      integer  i, j, k, n
c
      REAL_T dhx, dhy, dhz, cf0, cf1, cf2, cf3, cf4, cf5
      REAL_T delta, gamma
      REAL_T, allocatable :: phit(:,:,:)
      allocate(phit(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
c
      dhx = beta/h(1)**2
      dhy = beta/h(2)**2
      dhz = beta/h(3)**2

      do n = 1, nc
          do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1),hi(1)

                  cf0 = merge(f0(lo(1),j,k), 0.0D0,
     $                 (i .eq. lo(1)) .and. (m0(lo(1)-1,j,k).gt.0))
                  cf1 = merge(f1(i,lo(2),k), 0.D00,
     $                 (j .eq. lo(2)) .and. (m1(i,lo(2)-1,k).gt.0))
                  cf2 = merge(f2(i,j,lo(3)), 0.0D0,
     $                 (k .eq. lo(3)) .and. (m2(i,j,lo(3)-1).gt.0))
                  cf3 = merge(f3(hi(1),j,k), 0.0D0,
     $                 (i .eq. hi(1)) .and. (m3(hi(1)+1,j,k).gt.0))
                  cf4 = merge(f4(i,hi(2),k), 0.0D0,
     $                 (j .eq. hi(2)) .and. (m4(i,hi(2)+1,k).gt.0))
                  cf5 = merge(f5(i,j,hi(3)), 0.0D0,
     $                 (k .eq. hi(3)) .and. (m5(i,j,hi(3)+1).gt.0))

                  delta = dhx*(bX(i,j,k)*cf0 + bX(i+1,j,k)*cf3)
     $                 +  dhy*(bY(i,j,k)*cf1 + bY(i,j+1,k)*cf4)
     $                 +  dhz*(bZ(i,j,k)*cf2 + bZ(i,j,k+1)*cf5)

                  gamma = alpha*a(i,j,k)
     $                 +   dhx*(bX(i,j,k)+bX(i+1,j,k))
     $                 +   dhy*(bY(i,j,k)+bY(i,j+1,k))
     $                 +   dhz*(bZ(i,j,k)+bZ(i,j,k+1))

                  phit = phi(i,j,k,n) + res(i,j,k,n) / (gamma - delta)

               end do
            end do
          end do
          phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),n) = phit
      end do

      deallocate(phit)
      end

      subroutine ab2_bndryrlx(lo, hi,
     &     den, f_l1,f_l2,f_l3,f_h1,f_h2,f_h3,
     &     mask,m_l1,m_l2,m_l3,m_h1,m_h2,m_h3,
     &     cdir, bct, bcl, maxorder, h)

      use amrex_lo_util_module, only : polyInterpCoeff
      implicit none
      integer lo(3), hi(3)
      integer f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
      integer m_l1,m_l2,m_l3,m_h1,m_h2,m_h3
      
      REAL_T   den(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3)
      integer mask(m_l1:m_h1,m_l2:m_h2,m_l3:m_h3)

      integer cdir, bct, maxorder
      REAL_T bcl
      REAL_T h(3)
c
      integer i, j, k
      logical is_dirichlet
      logical is_neumann
c
      integer m, lenx, leny, lenz
c
      integer Lmaxorder
      integer maxmaxorder
      parameter(maxmaxorder=4)
      REAL_T x(-1:maxmaxorder-2)
      REAL_T coef(-1:maxmaxorder-2)
      REAL_T xInt

      is_dirichlet(i) = ( i .eq. LO_DIRICHLET )
      is_neumann(i) = (i .eq. LO_NEUMANN)

      if ( maxorder .eq. -1 ) then
         Lmaxorder = maxmaxorder
      else
         Lmaxorder = MIN(maxorder,maxmaxorder)
      end if
      lenx = MIN(hi(1)-lo(1), Lmaxorder-2)
      leny = MIN(hi(2)-lo(2), Lmaxorder-2)
      lenz = MIN(hi(3)-lo(3), Lmaxorder-2)

      do m=0,maxmaxorder-2
         x(m) = m + 0.5D0
      end do

      select case (cdir)

      case (0)
         !     
         ! The Left face of the grid
         !
         if (is_neumann(bct)) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  den(lo(1),j,k) = 1.0D0
               end do
            end do
         else if (is_dirichlet(bct)) then
            x(-1) = - bcl/h(1)
            xInt = -0.5D0
            call polyInterpCoeff(xInt, x, lenx+2, coef)
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  den(lo(1),j,k) = merge(coef(0), 0.0D0,
     $                 mask(lo(1)-1,j,k) .gt. 0)
               end do
            end do
         else if ( bct .eq. LO_REFLECT_ODD ) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  den(lo(1),j,k) = merge(-1.0D0, 0.0D0,
     $                 mask(lo(1)-1,j,k) .gt. 0)
               end do
            end do
         else
            print *,'UNKNOWN BC ON LEFT FACE IN APPLYBC'
            call bl_error("stop")
         end if

      case (3)
         !
         ! The Right face of the grid
         ! 
         if(is_neumann(bct)) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  den(hi(1),j,k) = 1.0D0
               end do
            end do
         else if (is_dirichlet(bct)) then
            x(-1) = - bcl/h(1)
            xInt = -0.5D0
            call polyInterpCoeff(xInt, x, lenx+2, coef)
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  den(hi(1),j,k)   = merge(coef(0), 0.0D0,
     $                 mask(hi(1)+1,j,k) .gt. 0)
               end do
            end do
         else if ( bct .eq. LO_REFLECT_ODD ) then
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  den(hi(1),j,k) = merge(-1.0D0, 0.0D0,
     $                 mask(hi(1)+1,j,k) .gt. 0)
               end do
            end do
         else
            print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
            call bl_error("stop")
         end if

         case (1)
            !
            ! The Bottom of the Grid
            !
            if(is_neumann(bct)) then
               do k = lo(3), hi(3)
                  do i = lo(1),hi(1)
                     den(i,lo(2),k)   = 1.0D0
                  end do
               end do
            else if (is_dirichlet(bct)) then
               x(-1) = - bcl/h(2)
               xInt = -0.5D0
               call polyInterpCoeff(xInt, x, leny+2, coef)
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     den(i, lo(2),k)   = merge(coef(0), 0.0D0,
     $                    mask(i, lo(2)-1,k) .gt. 0)
                  end do
               end do
            else if ( bct .eq. LO_REFLECT_ODD ) then
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     den(i,lo(2),k) = merge(-1.0D0, 0.0D0,
     $                    mask(i,lo(2)-1,k) .gt. 0)
                  end do
               end do
            else
               print *,'UNKNOWN BC ON BOTTOM FACE IN APPLYBC'
               call bl_error("stop")
            end if

         case (4)
            !
            ! The top of the grid
            !
            if(is_neumann(bct)) then
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     den(i,hi(2),k)   = 1.0D0
                  end do
               end do
            else if (is_dirichlet(bct)) then
               x(-1) = - bcl/h(2)
               xInt = -0.5D0
               call polyInterpCoeff(xInt, x, leny+2, coef)
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     den(i,hi(2),k)   = merge(coef(0), 0.0D0,
     $                    mask(i,hi(2)+1,k) .gt. 0)
                  end do
               end do
            else if ( bct .eq. LO_REFLECT_ODD ) then
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     den(i,hi(2),k) = merge(-1.0D0, 0.0D0,
     $                    mask(i,hi(2)+1,k) .gt. 0)
                  end do
               end do
            else
               print *,'UNKNOWN BC ON TOP FACE IN APPLYBC'
               call bl_error("stop")
            end if

         case (2)
            !
            ! The Front of the Grid
            !
            if(is_neumann(bct)) then
               do j = lo(2), hi(2)
                  do i = lo(1),hi(1)
                     den(i,j,lo(3))   = 1.0D0
                  end do
               end do
            else if (is_dirichlet(bct)) then
               x(-1) = - bcl/h(3)
               xInt = -0.5D0
               call polyInterpCoeff(xInt, x, lenz+2, coef)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     den(i, j, lo(3))   = merge(coef(0), 0.0D0,
     $                    mask(i, j, lo(3)-1) .gt. 0)
                  end do
               end do
            else if ( bct .eq. LO_REFLECT_ODD ) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     den(i,j,lo(3)) = merge(-1.0D0, 0.0D0,
     $                    mask(i,j,lo(3)-1) .gt. 0)
                  end do
               end do
            else
               print *,'UNKNOWN BC ON FRONT FACE IN APPLYBC'
               call bl_error("stop")
            end if

         case (5)
            !
            ! The back of the grid
            !
            if(is_neumann(bct)) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     den(i,j, hi(3))   = 1.0D0
                  end do
               end do
            else if (is_dirichlet(bct)) then
               x(-1) = - bcl/h(3)
               xInt = -0.5D0
               call polyInterpCoeff(xInt, x, lenz+2, coef)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     den(i,j, hi(3))   = merge(coef(0), 0.0D0,
     $                    mask(i,j, hi(3)+1) .gt. 0)
                  end do
               end do
            else if ( bct .eq. LO_REFLECT_ODD ) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     den(i,j,hi(3)) = merge(-1.0D0, 0.0D0,
     $                    mask(i,j,hi(3)+1) .gt. 0)
                  end do
               end do
            else
               print *,'UNKNOWN BC ON BACK FACE IN APPLYBC'
               call bl_error("stop")
            end if
         end select

      end

