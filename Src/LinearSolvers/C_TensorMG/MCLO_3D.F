
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
      integer k
      integer n
c
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  res(i,j,k,n) = rhs(i,j,k,n) - phi(i,j,k,n)
               enddo
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
      REAL_T  factor
      parameter(factor=4.0d0)

      integer n
      integer i
      integer j
      integer k
c
      if ( cdir .eq. 0 ) then
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
c                  
                     c(i,j,k,n) = factor/(
     $                    + 1.0/f(2*i,2*j  ,2*k  ,n)
     $                    + 1.0/f(2*i,2*j+1,2*k  ,n)
     $                    + 1.0/f(2*i,2*j  ,2*k+1,n)
     $                    + 1.0/f(2*i,2*j+1,2*k+1,n) )
c                     
                  enddo
               enddo
            enddo
         enddo
      else if (cdir .eq. 1 ) then
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
c
                     c(i,j,k,n) = factor/(
     $                    + 1.0/f(2*i  ,2*j,2*k  ,n)
     $                    + 1.0/f(2*i+1,2*j,2*k  ,n)
     $                    + 1.0/f(2*i  ,2*j,2*k+1,n)
     $                    + 1.0/f(2*i+1,2*j,2*k+1,n) )
c                     
                  enddo
               enddo
            enddo
         enddo
      else if (cdir .eq. 2 ) then
         do n = 1, nc
            do k = lo(3), hi(3)+1
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
c                     
                     c(i,j,k,n) = factor/(
     $                    + 1.0/f(2*i  ,2*j  ,2*k,n)
     $                    + 1.0/f(2*i+1,2*j  ,2*k,n)
     $                    + 1.0/f(2*i  ,2*j+1,2*k,n)
     $                    + 1.0/f(2*i+1,2*j+1,2*k,n) )
c                     
                  enddo
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
      integer k
      REAL_T factor
      parameter(factor = 0.25d0)
c     
      if ( cdir .eq. 0 ) then
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
c     
                     c(i,j,k,n) = factor*(
     $                    + f(2*i,2*j  ,2*k  ,n)
     $                    + f(2*i,2*j+1,2*k  ,n)
     $                    + f(2*i,2*j  ,2*k+1,n)
     $                    + f(2*i,2*j+1,2*k+1,n) )
c                     
                  enddo
               enddo
            enddo
         enddo
      else if (cdir .eq. 1 ) then
         do n = 1, nc
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
c                     
                     c(i,j,k,n) = factor*(
     $                    + f(2*i  ,2*j,2*k  ,n)
     $                    + f(2*i+1,2*j,2*k  ,n)
     $                    + f(2*i  ,2*j,2*k+1,n)
     $                    + f(2*i+1,2*j,2*k+1,n) )
c                     
                  enddo
               enddo
            enddo
         enddo
      else if (cdir .eq. 2 ) then
         do n = 1, nc
            do k = lo(3), hi(3)+1
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
c     
                     c(i,j,k,n) = factor*(
     $                    + f(2*i  ,2*j  ,2*k,n)
     $                    + f(2*i+1,2*j  ,2*k,n)
     $                    + f(2*i  ,2*j+1,2*k,n)
     $                    + f(2*i+1,2*j+1,2*k,n) )
c                     
                  enddo
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
      integer k
      integer n
      REAL_T factor
      parameter(factor=0.125d0)
c
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
c                  
                  c(i,j,k,n) =  factor*(
     $                 + f(2*i+1,2*j+1,2*k  ,n)
     $                 + f(2*i  ,2*j+1,2*k  ,n)
     $                 + f(2*i+1,2*j  ,2*k  ,n)
     $                 + f(2*i  ,2*j  ,2*k  ,n)
     $                 + f(2*i+1,2*j+1,2*k+1,n)
     $                 + f(2*i  ,2*j+1,2*k+1,n)
     $                 + f(2*i+1,2*j  ,2*k+1,n)
     $                 + f(2*i  ,2*j  ,2*k+1,n) )
c                     
               enddo
            enddo
         enddo
      enddo
c     
      end
