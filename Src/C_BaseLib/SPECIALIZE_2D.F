c
c $Id: SPECIALIZE_2D.F,v 1.2 1998-02-22 19:07:01 car Exp $
c

#ifdef BL_USE_SPECIALIZE

#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include <CONSTANTS.H>
#include <SPECIALIZE_F.H>
#include <ArrayLim.H>

#define SDIM 2

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
      integer i,j,k,nx,ny,ioff,joff

      ioff=imn-imin
      joff=jmn-jmin

      nx = imax-imin+1
      ny = jmax-jmin+1

      if (nx .ge. ny) then
         do k = 1, ncomp
            do j = jmin,jmax
               do i = imin,imax
                  dest(i,j,k) = src(i-ioff,j-joff,k)
               end do
            end do
         end do
      else
         do k = 1, ncomp
            do i = imin,imax
               do j = jmin,jmax
                  dest(i,j,k) = src(i-ioff,j-joff,k)
               end do
            end do
         end do
      end if

      end
c
c     This function sets a section of an array to a value.
c
      subroutine FORT_FASTSETVAL (val,lo,hi,dest,DIMS(dest),ncomp)

      implicit none

      integer ncomp
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(dest)
      REAL_T  val
      REAL_T  dest(DIMV(dest),ncomp)
c
c     Local variables
c 
      integer i,j,k,nx,ny
      integer imin,jmin,imax,jmax

      imin = lo(1)
      jmin = lo(2)
      imax = hi(1)
      jmax = hi(2)

      nx = imax-imin+1
      ny = jmax-jmin+1

      if (nx .ge. ny) then
         do k = 1, ncomp
            do j = jmin,jmax
               do i = imin,imax
                  dest(i,j,k) = val
               end do
            end do
         end do
      else
         do k = 1, ncomp
            do i = imin,imax
               do j = jmin,jmax
                  dest(i,j,k) = val
               end do
            end do
         end do
      end if

      end

#endif /*BL_USE_SPECIALIZE*/
