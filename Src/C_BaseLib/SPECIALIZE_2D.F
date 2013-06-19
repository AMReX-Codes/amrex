
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
     &                          imn,  jmn,
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
      integer imn,  jmn
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
c     This function copies from a 2D array to a 1D one.
c
      subroutine FORT_FASTCOPYTOMEM (lo,hi,data,DIMS(data),ncomp,dst)

      implicit none

      integer ncomp
      integer lo(2), hi(2)
      integer DIMDEC(data)
      REAL_T  data(DIMV(data),ncomp), dst(*)
c
c     Local variables.
c
      integer i,j,k,cnt
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      cnt = 1

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               dst(cnt) = data(i,j,k)
               cnt = cnt + 1
            end do
         end do
      end do

      end
c
c     This function copies from a 1D array into a 2D array.
c
      subroutine FORT_FASTCOPYFROMMEM (lo,hi,data,DIMS(data),ncomp,src)

      implicit none

      integer ncomp
      integer lo(2), hi(2)
      integer DIMDEC(data)
      REAL_T  data(DIMV(data),ncomp), src(*)
c
c     Local variables.
c
      integer i,j,k,cnt
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      cnt = 1

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               data(i,j,k) = src(cnt)
               cnt = cnt + 1
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
      integer i,j,k
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      sm = 0.0d0

      do k = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               sm = sm + src(i,j,k)
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
     &                          imn,  jmn,
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
      integer imn,  jmn
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
     &                          imn,  jmn,
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
      integer imn,  jmn
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
c     This function adds scaled floating point numbers from one array to another.
c
      subroutine FORT_FASTSAXPY (dest,
     &                           DIMS(dest),
     &                           imin, jmin, imax, jmax,
     &                           a, src,
     &                           DIMS(src),
     &                           imn,  jmn,
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
      integer imn,  jmn
      integer DIMDEC(src)
      integer ncomp
      
      REAL_T  dest(DIMV(dest),ncomp)
      REAL_T  a
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,j,k,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      do k = 1, ncomp
         do j = jmin, jmax
            do i = imin, imax
               dest(i,j,k) = dest(i,j,k) + a * src(i+ioff,j+joff,k)
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
     &                           imn,  jmn,
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
      integer imn,  jmn
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
c     This function divides floating point numbers from one array by another
c
      subroutine FORT_FASTDIVIDE (dest,
     &                            DIMS(dest),
     &                            imin, jmin, imax, jmax,
     &                            src,
     &                            DIMS(src),
     &                            imn,  jmn,
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
      integer imn,  jmn
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
c
c     This function divides floating point numbers from one array by another
c     wherever the denominator is non-zero.
c
      subroutine FORT_FASTPROTDIVIDE (dest,
     &                                DIMS(dest),
     &                                imin, jmin, imax, jmax,
     &                                src,
     &                                DIMS(src),
     &                                imn,  jmn,
     &                                ncomp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, jmin, imax, jmax
      integer DIMDEC(dest)
c
c     Bounds to fill from src
c
      integer imn,  jmn
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
               if (src(i+ioff,j+joff,k) .ne. zero) then
                  dest(i,j,k) = dest(i,j,k) / src(i+ioff,j+joff,k)
               end if
            end do
         end do
      end do

      end

      subroutine FORT_FASTINVERT (src,DIMS(src),lo,hi,val,ncomp)

      implicit none

      integer ncomp
      integer lo(2), hi(2)
      integer DIMDEC(src)
      REAL_T  src(DIMV(src),ncomp), val
c
c     Local variables.
c
      integer i,j,k
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      do k = 1,ncomp
         do j = jmin,jmax
            do i = imin,imax
               src(i,j,k) = val / src(i,j,k)
            end do
         end do
      end do

      end


c
c     This function does quadratic interpolation 
c
      subroutine FORT_QUADINTERP (dest,
     &                            DIMS(dest),
     &                            imin, jmin, imax, jmax,
     &                            src1,
     &                            DIMS(src1),
     &                            src2,
     &                            DIMS(src2),
     &                            src3,
     &                            DIMS(src3),
     &                            ncomp,
     &                            t,
     &                            t1,
     &                            t2,
     &                            t3)

      implicit none
c
c     Bounds to do interpolation in dest
c
      integer imin, jmin, imax, jmax
      integer DIMDEC(dest)
c
c     Bounds of src
c
      integer DIMDEC(src1)
      integer DIMDEC(src2)
      integer DIMDEC(src3)
      integer ncomp
      
      REAL_T  dest(DIMV(dest),ncomp)
      REAL_T  src1(DIMV(src1),ncomp)
      REAL_T  src2(DIMV(src2),ncomp)
      REAL_T  src3(DIMV(src3),ncomp)

      REAL_T  t, t1, t2, t3
c
c     Local variables
c
      integer i,j,n

      do n = 1, ncomp
         do j = jmin,jmax
            do i = imin,imax
               dest(i,j,n) = 0.d0
            end do
         end do
      end do

      end
