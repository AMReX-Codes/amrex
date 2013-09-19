
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include <CONSTANTS.H>
#include "MG_F.H"
#include "ArrayLim.H"

      subroutine FORT_AVERAGE (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc)
      implicit none
      integer nc
      integer DIMDEC(c)
      integer DIMDEC(f)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)

      integer i, i2, i2p1, j, j2, j2p1, k, k2, k2p1, n

!$omp parallel private(i,j,k,k2,k2p1,j2,j2p1,i2,i2p1,n)
      do n = 1, nc
!$omp do
         do k = lo(3), hi(3)
            k2 = 2*k
            k2p1 = k2 + 1
	    do j = lo(2), hi(2)
               j2 = 2*j
               j2p1 = j2 + 1
               do i = lo(1), hi(1)
                  i2 = 2*i
                  i2p1 = i2 + 1
                  c(i,j,k,n) =  (
     $                 + f(i2p1,j2p1,k2  ,n) + f(i2,j2p1,k2  ,n)
     $                 + f(i2p1,j2  ,k2  ,n) + f(i2,j2  ,k2  ,n)
     $                 + f(i2p1,j2p1,k2p1,n) + f(i2,j2p1,k2p1,n)
     $                 + f(i2p1,j2  ,k2p1,n) + f(i2,j2  ,k2p1,n)
     $                 )*eighth
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel

      end
!
! This can't be OpenMP'd.
! Instead we thread over FABs in MultiGrid::interpolate()
!
      subroutine FORT_INTERP (
     $     f, DIMS(f),
     $     c, DIMS(c),
     $     lo, hi, nc)
      implicit none
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)

      integer :: i, j, k, n, twoi, twoj, twoip1, twojp1, twok, twokp1
      logical :: interior_i, interior_j, interior_k

      REAL_T, parameter ::   ONE64TH  = 1.0d0 / 64.0d0
      REAL_T, parameter :: THREE64THS = 3.0d0 / 64.0d0
      !
      ! Trilinear interpolotion.
      ! Don't assume we have any grow cells.
      ! First do all face points using piecewise-constant interpolation.
      !
      do n = 1, nc
         do k = lo(3),hi(3)
            twok   = 2*k
            twokp1 = twok+1
            interior_k = ( k > lo(3) .and. k < hi(3) )

            do j = lo(2),hi(2)
               twoj   = 2*j
               twojp1 = twoj+1
               interior_j = ( j > lo(2) .and. j < hi(2) )

               do i = lo(1),hi(1)
                  interior_i = ( i > lo(1) .and. i < hi(1) )

                  if ( interior_i .and. interior_j .and. interior_k ) cycle

                  twoi   = 2*i
                  twoip1 = twoi+1

                  f(twoip1, twojp1, twokp1,n) = f(twoip1, twojp1, twokp1,n) + c(i,j,k,n)
                  f(twoi,   twojp1, twokp1,n) = f(twoi,   twojp1, twokp1,n) + c(i,j,k,n)
                  f(twoip1, twoj,   twokp1,n) = f(twoip1, twoj,   twokp1,n) + c(i,j,k,n)
                  f(twoi,   twoj,   twokp1,n) = f(twoi,   twoj,   twokp1,n) + c(i,j,k,n)
                  f(twoip1, twojp1, twok  ,n) = f(twoip1, twojp1, twok  ,n) + c(i,j,k,n)
                  f(twoi,   twojp1, twok  ,n) = f(twoi,   twojp1, twok  ,n) + c(i,j,k,n)
                  f(twoip1, twoj,   twok  ,n) = f(twoip1, twoj,   twok  ,n) + c(i,j,k,n)
                  f(twoi,   twoj,   twok  ,n) = f(twoi,   twoj,   twok  ,n) + c(i,j,k,n)
               end do
            end do
         end do
      end do
      !
      ! Now linearly interp the interior.
      !
      do n = 1, nc
         do k = lo(3)+1,hi(3)-1
            twok   = 2*k
            twokp1 = twok+1

            do j = lo(2)+1,hi(2)-1
               twoj   = 2*j
               twojp1 = twoj+1

               do i = lo(1)+1,hi(1)-1
                  twoi   = 2*i
                  twoip1 = twoi+1

                  f(twoip1, twojp1, twokp1,n) = f(twoip1, twojp1, twokp1,n) +
     &                 THREE64THS * ( 9*c(i,j,k  ,n) + 3*c(i+1,j,k  ,n) + 3*c(i,j+1,k  ,n) + c(i+1,j+1,k  ,n) ) +
     &                 ONE64TH    * ( 9*c(i,j,k+1,n) + 3*c(i+1,j,k+1,n) + 3*c(i,j+1,k+1,n) + c(i+1,j+1,k+1,n) )
                  f(twoi,   twojp1, twokp1,n) = f(twoi,   twojp1, twokp1,n) +
     &                 THREE64THS * ( 9*c(i,j,k  ,n) + 3*c(i-1,j,k  ,n) + 3*c(i,j+1,k  ,n) + c(i-1,j+1,k  ,n) ) +
     &                 ONE64TH    * ( 9*c(i,j,k+1,n) + 3*c(i-1,j,k+1,n) + 3*c(i,j+1,k+1,n) + c(i-1,j+1,k+1,n) )
                  f(twoip1, twoj,   twokp1,n) = f(twoip1, twoj,   twokp1,n) +
     &                 THREE64THS * ( 9*c(i,j,k  ,n) + 3*c(i+1,j,k  ,n) + 3*c(i,j-1,k  ,n) + c(i+1,j-1,k  ,n) ) +
     &                 ONE64TH    * ( 9*c(i,j,k+1,n) + 3*c(i+1,j,k+1,n) + 3*c(i,j-1,k+1,n) + c(i+1,j-1,k+1,n) )
                  f(twoi,   twoj,   twokp1,n) = f(twoi,   twoj,   twokp1,n) +
     &                 THREE64THS * ( 9*c(i,j,k  ,n) + 3*c(i-1,j,k  ,n) + 3*c(i,j-1,k  ,n) + c(i-1,j-1,k  ,n) ) + 
     &                 ONE64TH    * ( 9*c(i,j,k+1,n) + 3*c(i-1,j,k+1,n) + 3*c(i,j-1,k+1,n) + c(i-1,j-1,k+1,n) )
                  f(twoip1, twojp1, twok,n) = f(twoip1, twojp1, twok,n) +
     &                 THREE64THS * ( 9*c(i,j,k  ,n) + 3*c(i+1,j,k  ,n) + 3*c(i,j+1,k  ,n) + c(i+1,j+1,k  ,n) ) +
     &                 ONE64TH    * ( 9*c(i,j,k-1,n) + 3*c(i+1,j,k-1,n) + 3*c(i,j+1,k-1,n) + c(i+1,j+1,k-1,n) )
                  f(twoi,   twojp1, twok,n) = f(twoi,   twojp1, twok,n) +
     &                 THREE64THS * ( 9*c(i,j,k  ,n) + 3*c(i-1,j,k  ,n) + 3*c(i,j+1,k  ,n) + c(i-1,j+1,k  ,n) ) +
     &                 ONE64TH    * ( 9*c(i,j,k-1,n) + 3*c(i-1,j,k-1,n) + 3*c(i,j+1,k-1,n) + c(i-1,j+1,k-1,n) )
                  f(twoip1, twoj,   twok,n) = f(twoip1, twoj,   twok,n) +
     &                 THREE64THS * ( 9*c(i,j,k  ,n) + 3*c(i+1,j,k  ,n) + 3*c(i,j-1,k  ,n) + c(i+1,j-1,k  ,n) ) +
     &                 ONE64TH    * ( 9*c(i,j,k-1,n) + 3*c(i+1,j,k-1,n) + 3*c(i,j-1,k-1,n) + c(i+1,j-1,k-1,n) )
                  f(twoi,   twoj,   twok,n) = f(twoi,   twoj,   twok,n) +
     &                 THREE64THS * ( 9*c(i,j,k  ,n) + 3*c(i-1,j,k  ,n) + 3*c(i,j-1,k  ,n) + c(i-1,j-1,k  ,n) ) +
     &                 ONE64TH    * ( 9*c(i,j,k-1,n) + 3*c(i-1,j,k-1,n) + 3*c(i,j-1,k-1,n) + c(i-1,j-1,k-1,n) )
               end do
            end do
         end do
      end do

#if 0
      integer i, i2, i2p1, j, j2, j2p1, k, k2, k2p1, n

      do n = 1, nc
         do k = lo(3), hi(3)
            k2 = 2*k
            k2p1 = k2 + 1
	    do j = lo(2), hi(2)
               j2 = 2*j
               j2p1 = j2 + 1

               do i = lo(1), hi(1)
                  i2 = 2*i
                  i2p1 = i2 + 1

                  f(i2p1,j2p1,k2  ,n) = c(i,j,k,n) + f(i2p1,j2p1,k2  ,n)
                  f(i2  ,j2p1,k2  ,n) = c(i,j,k,n) + f(i2  ,j2p1,k2  ,n)
                  f(i2p1,j2  ,k2  ,n) = c(i,j,k,n) + f(i2p1,j2  ,k2  ,n)
                  f(i2  ,j2  ,k2  ,n) = c(i,j,k,n) + f(i2  ,j2  ,k2  ,n)
                  f(i2p1,j2p1,k2p1,n) = c(i,j,k,n) + f(i2p1,j2p1,k2p1,n)
                  f(i2  ,j2p1,k2p1,n) = c(i,j,k,n) + f(i2  ,j2p1,k2p1,n)
                  f(i2p1,j2  ,k2p1,n) = c(i,j,k,n) + f(i2p1,j2  ,k2p1,n)
                  f(i2  ,j2  ,k2p1,n) = c(i,j,k,n) + f(i2  ,j2  ,k2p1,n)

               end do
            end do
         end do
      end do
#endif

      end
