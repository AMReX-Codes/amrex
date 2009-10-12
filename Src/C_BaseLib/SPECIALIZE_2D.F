
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
     &                          imin, jmin, imax, jmax,
     &                          src,
     &                          DIMS(src),
     &                          imn,  jmn,  imx,  jmx,
     &                          ncomp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, jmin, imax, jmax
      integer DIMDEC(dest)
c
c     Bounds to fill from src
c
      integer imn,  jmn,  imx,  jmx
      integer DIMDEC(src)
      integer ncomp
      
      REAL_T  dest(DIMV(dest),ncomp)
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,j,k,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               dest(i,j,k) = src(i+ioff,j+joff,k)
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
      integer lo(2), hi(2)
      integer DIMDEC(dest)
      REAL_T  val
      REAL_T  dest(DIMV(dest),ncomp)
c
c     Local variables
c 
      integer i,j,k
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               dest(i,j,k) = val
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
      integer i,j,k
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      nrm = 0.0d0

      do k = 1,ncomp
         do j = jmin,jmax
            do i = imin,imax
               nrm = max(nrm, abs(src(i,j,k)))
            end do
         end do
      end do

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
      integer i,j,k
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      nrm = 0.0d0

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               nrm = nrm + abs(src(i,j,k))
            end do
         end do
      end do

      end
c
c     This function adds floating point numbers from one array to another.
c
      subroutine FORT_FASTPLUS (dest,
     &                          DIMS(dest),
     &                          imin, jmin, imax, jmax,
     &                          src,
     &                          DIMS(src),
     &                          imn,  jmn,  imx,  jmx,
     &                          ncomp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, jmin, imax, jmax
      integer DIMDEC(dest)
c
c     Bounds to fill from src
c
      integer imn,  jmn,  imx,  jmx
      integer DIMDEC(src)
      integer ncomp
      
      REAL_T  dest(DIMV(dest),ncomp)
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,j,k,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               dest(i,j,k) = dest(i,j,k) + src(i+ioff,j+joff,k)
            end do
         end do
      end do

      end
c
c     This function multiplys floating point numbers from one array to another.
c
      subroutine FORT_FASTMULT (dest,
     &                          DIMS(dest),
     &                          imin, jmin, imax, jmax,
     &                          src,
     &                          DIMS(src),
     &                          imn,  jmn,  imx,  jmx,
     &                          ncomp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, jmin, imax, jmax
      integer DIMDEC(dest)
c
c     Bounds to fill from src
c
      integer imn,  jmn,  imx,  jmx
      integer DIMDEC(src)
      integer ncomp
      
      REAL_T  dest(DIMV(dest),ncomp)
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,j,k,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               dest(i,j,k) = dest(i,j,k) * src(i+ioff,j+joff,k)
            end do
         end do
      end do

      end
c
c     This function subtracts floating point numbers from one array to another.
c
      subroutine FORT_FASTMINUS (dest,
     &                           DIMS(dest),
     &                           imin, jmin, imax, jmax,
     &                           src,
     &                           DIMS(src),
     &                           imn,  jmn,  imx,  jmx,
     &                           ncomp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, jmin, imax, jmax
      integer DIMDEC(dest)
c
c     Bounds to fill from src
c
      integer imn,  jmn,  imx,  jmx
      integer DIMDEC(src)
      integer ncomp
      
      REAL_T  dest(DIMV(dest),ncomp)
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,j,k,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               dest(i,j,k) = dest(i,j,k) - src(i+ioff,j+joff,k)
            end do
         end do
      end do

      end
c
c     This function divides floating point numbers from one array by another.
c
      subroutine FORT_FASTDIVIDE (dest,
     &                            DIMS(dest),
     &                            imin, jmin, imax, jmax,
     &                            src,
     &                            DIMS(src),
     &                            imn,  jmn,  imx,  jmx,
     &                            ncomp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, jmin, imax, jmax
      integer DIMDEC(dest)
c
c     Bounds to fill from src
c
      integer imn,  jmn,  imx,  jmx
      integer DIMDEC(src)
      integer ncomp
      
      REAL_T  dest(DIMV(dest),ncomp)
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,j,k,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               dest(i,j,k) = dest(i,j,k) / src(i+ioff,j+joff,k)
            end do
         end do
      end do

      end
