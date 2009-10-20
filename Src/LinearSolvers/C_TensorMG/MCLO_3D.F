
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

      integer i, j, k, n

      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  res(i,j,k,n) = rhs(i,j,k,n) - phi(i,j,k,n)
               enddo
            enddo
         enddo
      enddo
     
      end
c-----------------------------------------------------------------------
      subroutine FORT_HARMONIC_AVERAGEEC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc,
     $     cdir
     $     )

      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer DIMDEC(f)
      REAL_T f(DIMV(f),nc)
      integer DIMDEC(c)
      REAL_T c(DIMV(c),nc)

      integer n, i, j, k

!$omp parallel private(i,j,k)

      if ( cdir .eq. 0 ) then
         do n = 1, nc
!$omp do
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
                  
                     c(i,j,k,n) = 4.0d0/(
     $                    + 1.0d0/f(2*i,2*j  ,2*k  ,n)
     $                    + 1.0d0/f(2*i,2*j+1,2*k  ,n)
     $                    + 1.0d0/f(2*i,2*j  ,2*k+1,n)
     $                    + 1.0d0/f(2*i,2*j+1,2*k+1,n) )
                     
                  enddo
               enddo
            enddo
!$omp end do nowait
         enddo
      else if (cdir .eq. 1 ) then
         do n = 1, nc
!$omp do
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)

                     c(i,j,k,n) = 4.0d0/(
     $                    + 1.0d0/f(2*i  ,2*j,2*k  ,n)
     $                    + 1.0d0/f(2*i+1,2*j,2*k  ,n)
     $                    + 1.0d0/f(2*i  ,2*j,2*k+1,n)
     $                    + 1.0d0/f(2*i+1,2*j,2*k+1,n) )
                     
                  enddo
               enddo
            enddo
!$omp end do nowait
         enddo
      else if (cdir .eq. 2 ) then
         do n = 1, nc
!$omp do
            do k = lo(3), hi(3)+1
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     
                     c(i,j,k,n) = 4.0d0/(
     $                    + 1.0d0/f(2*i  ,2*j  ,2*k,n)
     $                    + 1.0d0/f(2*i+1,2*j  ,2*k,n)
     $                    + 1.0d0/f(2*i  ,2*j+1,2*k,n)
     $                    + 1.0d0/f(2*i+1,2*j+1,2*k,n) )
                     
                  enddo
               enddo
            enddo
!$omp end do nowait
         enddo
      endif
!$omp end parallel

      end
c-----------------------------------------------------------------------
      subroutine FORT_AVERAGEEC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc,
     $     cdir
     $     )

      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer DIMDEC(f)
      REAL_T f(DIMV(f),nc)
      integer DIMDEC(c)
      REAL_T c(DIMV(c),nc)
     
      integer n, i, j, k

!$omp parallel private(i,j,k)
     
      if ( cdir .eq. 0 ) then
         do n = 1, nc
!$omp do
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)+1
     
                     c(i,j,k,n) = 0.25d0*(
     $                    + f(2*i,2*j  ,2*k  ,n)
     $                    + f(2*i,2*j+1,2*k  ,n)
     $                    + f(2*i,2*j  ,2*k+1,n)
     $                    + f(2*i,2*j+1,2*k+1,n) )
                     
                  enddo
               enddo
            enddo
!$omp end do nowait
         enddo
      else if (cdir .eq. 1 ) then
         do n = 1, nc
!$omp do
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)+1
                  do i = lo(1), hi(1)
                     
                     c(i,j,k,n) = 0.25d0*(
     $                    + f(2*i  ,2*j,2*k  ,n)
     $                    + f(2*i+1,2*j,2*k  ,n)
     $                    + f(2*i  ,2*j,2*k+1,n)
     $                    + f(2*i+1,2*j,2*k+1,n) )
                     
                  enddo
               enddo
            enddo
!$omp end do nowait
         enddo
      else if (cdir .eq. 2 ) then
         do n = 1, nc
!$omp do
            do k = lo(3), hi(3)+1
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
     
                     c(i,j,k,n) = 0.25d0*(
     $                    + f(2*i  ,2*j  ,2*k,n)
     $                    + f(2*i+1,2*j  ,2*k,n)
     $                    + f(2*i  ,2*j+1,2*k,n)
     $                    + f(2*i+1,2*j+1,2*k,n) )
                     
                  enddo
               enddo
            enddo
!$omp end do nowait
         enddo
      endif
!$omp end parallel
     
      end
c-----------------------------------------------------------------------
      subroutine FORT_AVERAGECC (
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

      integer i, j, k, n

!$omp parallel private(i,j,k)
      do n = 1, nc
!$omp do
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  
                  c(i,j,k,n) =  0.125d0*(
     $                 + f(2*i+1,2*j+1,2*k  ,n)
     $                 + f(2*i  ,2*j+1,2*k  ,n)
     $                 + f(2*i+1,2*j  ,2*k  ,n)
     $                 + f(2*i  ,2*j  ,2*k  ,n)
     $                 + f(2*i+1,2*j+1,2*k+1,n)
     $                 + f(2*i  ,2*j+1,2*k+1,n)
     $                 + f(2*i+1,2*j  ,2*k+1,n)
     $                 + f(2*i  ,2*j  ,2*k+1,n) )
                     
               enddo
            enddo
         enddo
!$omp end do nowait
      enddo
!$omp end parallel
     
      end
