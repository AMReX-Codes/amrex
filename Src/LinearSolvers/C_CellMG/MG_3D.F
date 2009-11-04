c
c $Id: MG_3D.F,v 1.7 2009-11-04 16:55:37 lijewski Exp $
c
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

      end
