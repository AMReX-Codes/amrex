
#undef  BL_LANG_CC
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
