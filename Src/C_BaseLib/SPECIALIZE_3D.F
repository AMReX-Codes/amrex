
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
     &                          imin, jmin, kmin, imax, jmax, kmax,
     &                          src,
     &                          DIMS(src),
     &                          imn,  jmn,  kmn,  imx,  jmx,  kmx,
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
      integer imn,  jmn,  kmn,  imx,  jmx,  kmx
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
