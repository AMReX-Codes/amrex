
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include <CONSTANTS.H>
#include <SPECIALIZE_F.H>
#include <ArrayLim.H>

c
c     This function copies floating point numbers from one array to another.
c
      subroutine FORT_FASTCOPY (dest,
     &                          DIMS(dest),
     &                          imin, jmin, kmin, imax, jmax, kmax,
     &                          src,
     &                          DIMS(src),
     &                          imn,  jmn,  kmn,
     &                          ncomp)

      implicit none
c
c     Bounds to fill in dest 
c
      integer imin, jmin, kmin, imax, jmax, kmax
      integer DIMDEC(dest)
c
c     Bounds to fill from src
c
      integer imn,  jmn,  kmn
      integer DIMDEC(src)
      integer ncomp

      REAL_T  dest(DIMV(dest),ncomp)
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,j,k,l,ioff,joff,koff

      ioff=imn-imin
      joff=jmn-jmin
      koff=kmn-kmin

      do l = 1, ncomp
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dest(i,j,k,l) = src(i+ioff,j+joff,k+koff,l)
               end do
            end do
         end do
      end do

      end
c
c     This function copies from a 3D array to a 1D one.
c
      subroutine FORT_FASTCOPYTOMEM (lo,hi,data,DIMS(data),ncomp,dst)

      implicit none

      integer ncomp
      integer lo(3), hi(3)
      integer DIMDEC(data)
      REAL_T  data(DIMV(data),ncomp), dst(*)
c
c     Local variables.
c
      integer i,j,k,l,cnt
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      cnt = 1

      do l = 1, ncomp
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dst(cnt) = data(i,j,k,l)
                  cnt = cnt + 1
               end do
            end do
         end do
      end do

      end
c
c     This function copies from a 1D array into a 3D array.
c
      subroutine FORT_FASTCOPYFROMMEM (lo,hi,data,DIMS(data),ncomp,src)

      implicit none

      integer ncomp
      integer lo(3), hi(3)
      integer DIMDEC(data)
      REAL_T  data(DIMV(data),ncomp), src(*)
c
c     Local variables.
c
      integer i,j,k,l,cnt
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      cnt = 1

      do l = 1, ncomp
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  data(i,j,k,l) = src(cnt)
                  cnt = cnt + 1
               end do
            end do
         end do
      end do

      end
c
c     This function sets a section of an array to a value.
c
      subroutine FORT_FASTSETVAL (val,lo,hi,dest,DIMS(dest),ncomp)

      implicit none

      integer ncomp
      integer lo(3), hi(3)
      integer DIMDEC(dest)
      REAL_T  val
      REAL_T  dest(DIMV(dest),ncomp)
c
c     Local variables.
c
      integer i,j,k,l
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      do l = 1, ncomp
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dest(i,j,k,l) = val
               end do
            end do
         end do
      end do

      end
c
c     Calculate max(abs) over section of an array.
c
      subroutine FORT_FASTZERONORM (src,DIMS(src),lo,hi,ncomp,nrm)

      implicit none

      integer ncomp
      integer lo(3), hi(3)
      integer DIMDEC(src)
      REAL_T  src(DIMV(src),ncomp), nrm
c
c     Local variables.
c
      integer i,j,k,l
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      nrm = 0.0d0

!$omp parallel private(i,j,k,l) reduction(max : nrm)
      do l = 1, ncomp
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  nrm = max(nrm, abs(src(i,j,k,l)))
               end do
            end do
         end do
!$omp end do
      end do
!$omp end parallel

      end
c
c     Calculate sum(abs) over section of an array.
c
      subroutine FORT_FASTONENORM (src,DIMS(src),lo,hi,ncomp,nrm)

      implicit none

      integer ncomp
      integer lo(3), hi(3)
      integer DIMDEC(src)
      REAL_T  src(DIMV(src),ncomp), nrm
c
c     Local variables.
c
      integer i,j,k,l
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      nrm = 0.0d0

!$omp parallel private(i,j,k,l) reduction(+ : nrm)
      do l = 1, ncomp
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  nrm = nrm + abs(src(i,j,k,l))
               end do
            end do
         end do
!$omp end do
      end do
!$omp end parallel

      end
c
c     Calculate sum() over section of an array.
c
      subroutine FORT_FASTSUM (src,DIMS(src),lo,hi,ncomp,sm)

      implicit none

      integer ncomp
      integer lo(3), hi(3)
      integer DIMDEC(src)
      REAL_T  src(DIMV(src),ncomp), sm
c
c     Local variables.
c
      integer i,j,k,l
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      sm = 0.0d0

      do l = 1, ncomp
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  sm = sm + src(i,j,k,l)
               end do
            end do
         end do
      end do

      end
c
c     This function adds floating point numbers from one array to another.
c
      subroutine FORT_FASTPLUS (dst,
     &                          DIMS(dst),
     &                          imin, jmin, kmin, imax, jmax, kmax,
     &                          src,
     &                          DIMS(src),
     &                          imn,  jmn,  kmn,
     &                          ncomp)

      implicit none
c
c     Bounds to fill in dst 
c
      integer imin, jmin, kmin, imax, jmax, kmax
      integer DIMDEC(dst)
c
c     Bounds to fill from src
c
      integer imn,  jmn,  kmn
      integer DIMDEC(src)
      integer ncomp

      REAL_T  dst(DIMV(dst),ncomp)
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,j,k,l,ioff,joff,koff

      ioff=imn-imin
      joff=jmn-jmin
      koff=kmn-kmin

      do l = 1, ncomp
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dst(i,j,k,l) = dst(i,j,k,l) + src(i+ioff,j+joff,k+koff,l)
               end do
            end do
         end do
      end do

      end
c
c     This function multiplys floating point numbers from one array to another.
c
      subroutine FORT_FASTMULT (dst,
     &                          DIMS(dst),
     &                          imin, jmin, kmin, imax, jmax, kmax,
     &                          src,
     &                          DIMS(src),
     &                          imn,  jmn,  kmn,
     &                          ncomp)

      implicit none
c
c     Bounds to fill in dst 
c
      integer imin, jmin, kmin, imax, jmax, kmax
      integer DIMDEC(dst)
c
c     Bounds to fill from src
c
      integer imn,  jmn,  kmn
      integer DIMDEC(src)
      integer ncomp

      REAL_T  dst(DIMV(dst),ncomp)
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,j,k,l,ioff,joff,koff

      ioff=imn-imin
      joff=jmn-jmin
      koff=kmn-kmin

      do l = 1, ncomp
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dst(i,j,k,l) = dst(i,j,k,l) * src(i+ioff,j+joff,k+koff,l)
               end do
            end do
         end do
      end do

      end
c
c     This function adds scaled floating point numbers from one array to another.
c
      subroutine FORT_FASTSAXPY (dst,
     &                           DIMS(dst),
     &                           imin, jmin, kmin, imax, jmax, kmax,
     &                           a, src,
     &                           DIMS(src),
     &                           imn,  jmn,  kmn,
     &                           ncomp)

      implicit none
c
c     Bounds to fill in dst 
c
      integer imin, jmin, kmin, imax, jmax, kmax
      integer DIMDEC(dst)
c
c     Bounds to fill from src
c
      integer imn,  jmn,  kmn
      integer DIMDEC(src)
      integer ncomp

      REAL_T  dst(DIMV(dst),ncomp)
      REAL_T  a
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,j,k,l,ioff,joff,koff

      ioff=imn-imin
      joff=jmn-jmin
      koff=kmn-kmin

      do l = 1, ncomp
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dst(i,j,k,l) = dst(i,j,k,l) + a * src(i+ioff,j+joff,k+koff,l)
               end do
            end do
         end do
      end do

      end
c
c     This function subtracts floating point numbers from one array to another.
c
      subroutine FORT_FASTMINUS (dst,
     &                           DIMS(dst),
     &                           imin, jmin, kmin, imax, jmax, kmax,
     &                           src,
     &                           DIMS(src),
     &                           imn,  jmn,  kmn,
     &                           ncomp)

      implicit none
c
c     Bounds to fill in dst 
c
      integer imin, jmin, kmin, imax, jmax, kmax
      integer DIMDEC(dst)
c
c     Bounds to fill from src
c
      integer imn,  jmn,  kmn
      integer DIMDEC(src)
      integer ncomp

      REAL_T  dst(DIMV(dst),ncomp)
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,j,k,l,ioff,joff,koff

      ioff=imn-imin
      joff=jmn-jmin
      koff=kmn-kmin

      do l = 1, ncomp
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dst(i,j,k,l) = dst(i,j,k,l) - src(i+ioff,j+joff,k+koff,l)
               end do
            end do
         end do
      end do

      end
c
c     This function divides floating point numbers from one array by another.
c
      subroutine FORT_FASTDIVIDE (dst,
     &                            DIMS(dst),
     &                            imin, jmin, kmin, imax, jmax, kmax,
     &                            src,
     &                            DIMS(src),
     &                            imn,  jmn,  kmn,
     &                            ncomp)

      implicit none
c
c     Bounds to fill in dst 
c
      integer imin, jmin, kmin, imax, jmax, kmax
      integer DIMDEC(dst)
c
c     Bounds to fill from src
c
      integer imn,  jmn,  kmn
      integer DIMDEC(src)
      integer ncomp

      REAL_T  dst(DIMV(dst),ncomp)
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,j,k,l,ioff,joff,koff

      ioff=imn-imin
      joff=jmn-jmin
      koff=kmn-kmin

      do l = 1, ncomp
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  dst(i,j,k,l) = dst(i,j,k,l) / src(i+ioff,j+joff,k+koff,l)
               end do
            end do
         end do
      end do

      end
c
c     This function divides floating point numbers from one array by another
c     wherever the denominator is non-zero.
c
      subroutine FORT_FASTPROTDIVIDE (dst,
     &                                DIMS(dst),
     &                                imin, jmin, kmin, imax, jmax, kmax,
     &                                src,
     &                                DIMS(src),
     &                                imn,  jmn,  kmn,
     &                                ncomp)

      implicit none
c
c     Bounds to fill in dst 
c
      integer imin, jmin, kmin, imax, jmax, kmax
      integer DIMDEC(dst)
c
c     Bounds to fill from src
c
      integer imn,  jmn,  kmn
      integer DIMDEC(src)
      integer ncomp

      REAL_T  dst(DIMV(dst),ncomp)
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,j,k,l,ioff,joff,koff

      ioff=imn-imin
      joff=jmn-jmin
      koff=kmn-kmin

!$omp parallel private(i,j,k,l)
      do l = 1, ncomp
!$omp do
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  if (src(i+ioff,j+joff,k+koff,l) .ne. zero) then
                     dst(i,j,k,l) = dst(i,j,k,l) / src(i+ioff,j+joff,k+koff,l)
                  end if
               end do
            end do
         end do
!$omp end do nowait
      end do
!$omp end parallel

      end

      subroutine FORT_FASTINVERT (src,DIMS(src),lo,hi,val,ncomp)

      implicit none

      integer ncomp
      integer lo(3), hi(3)
      integer DIMDEC(src)
      REAL_T  src(DIMV(src),ncomp), val
c
c     Local variables.
c
      integer i,j,k,l
      integer imin,jmin,kmin,imax,jmax,kmax

      imin = lo(1)
      jmin = lo(2)
      kmin = lo(3)
      imax = hi(1)
      jmax = hi(2)
      kmax = hi(3)

      do l = 1, ncomp
         do k = kmin,kmax
            do j = jmin,jmax
               do i = imin,imax
                  src(i,j,k,l) = val / src(i,j,k,l)
               end do
            end do
         end do
      end do

      end

