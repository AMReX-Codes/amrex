c
c $Id: MG_3D.F,v 1.3 2003-09-18 19:21:00 lijewski Exp $
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
      integer nc
      integer DIMDEC(c)
      integer DIMDEC(f)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)

      integer i, i2, i2p1
      integer j, j2, j2p1
      integer k, k2, k2p1
      integer n
      REAL_T denom
      parameter(denom=eighth)

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

                  c(i,j,k,n) =  (
     $                 + f(i2p1,j2p1,k2  ,n) + f(i2,j2p1,k2  ,n)
     $                 + f(i2p1,j2  ,k2  ,n) + f(i2,j2  ,k2  ,n)
     $                 + f(i2p1,j2p1,k2p1,n) + f(i2,j2p1,k2p1,n)
     $                 + f(i2p1,j2  ,k2p1,n) + f(i2,j2  ,k2p1,n)
     $                 )*denom

               end do
            end do
         end do
      end do

      end

      subroutine FORT_INTERP (
     $     f, DIMS(f),
     $     c, DIMS(c),
     $     lo, hi, nc)
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)

      integer i, i2, i2p1
      integer j, j2, j2p1
      integer k, k2, k2p1
      integer n

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
