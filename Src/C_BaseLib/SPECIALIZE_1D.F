
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
     &                          imin, imax,
     &                          src,
     &                          DIMS(src),
     &                          imn, imx,
     &                          ncomp)

      implicit none
c
c     Bounds to fill in dest
c
      integer imin, imax
      integer DIMDEC(dest)
c
c     Bounds to fill from src
c
      integer imn,  imx
      integer DIMDEC(src)
      integer ncomp
      
      REAL_T  dest(DIMV(dest),ncomp)
      REAL_T  src(DIMV(src),ncomp)
c
c     Local variables
c
      integer i,k,ioff

      ioff=imn-imin

      do k = 1, ncomp
            do i = imin,imax
               dest(i,k) = src(i+ioff,k)
            end do
      end do

      end
c
c     This function sets a section of an array to a value.
c
      subroutine FORT_FASTSETVAL (val,lo,hi,dest,DIMS(dest),ncomp)

      implicit none

      integer ncomp
      integer lo(1), hi(1)
      integer DIMDEC(dest)
      REAL_T  val
      REAL_T  dest(DIMV(dest),ncomp)
c
c     Local variables
c 
      integer i,k
      integer imin,imax

      imin = lo(1)
      imax = hi(1)

      do k = 1, ncomp
            do i = imin,imax
               dest(i,k) = val
         end do
      end do

      end
