c
c $Id: MG_2D.F,v 1.1 1998-03-24 07:07:35 almgren Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include "MG_F.H"
#include "ArrayLim.H"

c-----------------------------------------------------------------------
      subroutine FORT_AVERAGE (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc
     $     )
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)
c
      integer i
      integer j
      integer n
      REAL_T denom
      parameter(denom=0.25d0)
c
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               c(i,j,n) =  (
     $              f(2*i+1,2*j+1,n) + f(2*i  ,2*j+1,n)
     $              + f(2*i+1,2*j,n ) + f(2*i  ,2*j ,n))*denom
            end do
         end do
      end do
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_INTERP (
     $     f, DIMS(f),
     $     c, DIMS(c),
     $     lo, hi, nc
     $     )
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)
c
      integer i
      integer j
      integer n
c
      do n = 1, nc
         do j = lo(2), hi(2)
cdir$   ivdep
            do i = lo(1), hi(1)
               f(2*i+1,2*j+1,n) = c(i,j,n) + f(2*i+1,2*j+1,n)
               f(2*i  ,2*j+1,n) = c(i,j,n) + f(2*i  ,2*j+1,n)
               f(2*i+1,2*j  ,n) = c(i,j,n) + f(2*i+1,2*j  ,n)
               f(2*i  ,2*j  ,n) = c(i,j,n) + f(2*i  ,2*j  ,n)
            end do
         end do
      end do
c
      end
