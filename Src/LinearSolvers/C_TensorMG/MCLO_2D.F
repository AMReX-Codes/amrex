c     $Id: MCLO_2D.F,v 1.1 1998-03-26 19:15:35 marc Exp $

#include <REAL.H>
#include "LO_BCTYPES.H"
#include "MCLO_F.H"
#include "ArrayLim.H"

c-----------------------------------------------------------------------
      subroutine FORT_RESIDL (
     $     res, DIMS(res), 
     $     rhs, DIMS(rhs),
     $     phi, DIMS(phi),
     $     lo, hi, nc
     $     )
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(phi)
      REAL_T phi(DIMV(phi),nc)
      integer DIMDEC(rhs)
      REAL_T rhs(DIMV(rhs),nc)
      integer DIMDEC(res)
      REAL_T res(DIMV(res),nc)
c
      integer i
      integer j
      integer n
c
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               res(i,j,n) = rhs(i,j,n) - phi(i,j,n)
            enddo
         enddo
      enddo
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_HARMONIC_AVERAGEEC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc,
     $     cdir
     $     )
c
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer DIMDEC(f)
      REAL_T f(DIMV(f),nc)
      integer DIMDEC(c)
      REAL_T c(DIMV(c),nc)
c
      REAL_T factor, den
      parameter(factor=2.00)
      integer n
      integer i
      integer j
c
      if ( cdir .eq. 0 ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)+1
                  den = f(2*i,2*j,n) + f(2*i,2*j+1,n)
                  if (den .ne. 0.0) then
                    c(i,j,n) =  factor*f(2*i,2*j,n)*f(2*i,2*j+1,n)/den
                  else
                    c(i,j,n) =  0.0
                  endif
               enddo
            enddo
         enddo
      else if (cdir .eq. 1 ) then
         do n = 1, nc
            do j = lo(2), hi(2)+1
               do i = lo(1), hi(1)
                  den = f(2*i,2*j,n) + f(2*i+1,2*j,n)
                  if (den .ne. 0.0) then
                    c(i,j,n) =  factor*f(2*i,2*j,n)*f(2*i+1,2*j,n)/den
                  else
                    c(i,j,n) =  0.0
                  endif
               enddo
            enddo
         enddo
      endif
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_AVERAGEEC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc,
     $     cdir
     $     )
c
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer DIMDEC(f)
      REAL_T f(DIMV(f),nc)
      integer DIMDEC(c)
      REAL_T c(DIMV(c),nc)
c
      integer n
      integer i
      integer j
      REAL_T denom
      parameter(denom=0.50d0)
c
      if (cdir .eq. 0 ) then
         do n = 1, nc
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)+1
                  c(i,j,n) = (f(2*i,2*j,n) + f(2*i,2*j+1,n))*denom
               enddo
            enddo
         enddo
      else if (cdir .eq. 1) then
         do n = 1, nc
            do j = lo(2), hi(2)+1
               do i = lo(1), hi(1)
                  c(i,j,n) = (f(2*i,2*j,n) + f(2*i+1,2*j,n))*denom
               enddo
            enddo
         enddo
      endif
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_AVERAGECC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc
     $     )
c
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
     $              + f(2*i+1,2*j  ,n) + f(2*i  ,2*j  ,n))*denom
            enddo
         enddo
      enddo
c
      end
