c
c $Id: SPECIALIZE_3D.F,v 1.1 1997-12-18 18:44:38 lijewski Exp $
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

#define SDIM 3

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
      integer i,j,k,l,nx,ny,nz,ioff,joff,koff

      integer minnx
      parameter (minnx = 64)

      ioff=imn-imin
      joff=jmn-jmin
      koff=kmn-kmin

      nx = imax-imin+1
      ny = jmax-jmin+1
      nz = kmax-kmin+1

      if ((nx .ge. ny .and. nx .ge. nz) .or. (nx .ge. minnx)) then
         do l = 1, ncomp
            do k = kmin,kmax
               do j = jmin,jmax
                  do i = imin,imax
                     dest(i,j,k,l) = src(i-ioff,j-joff,k-koff,l)
                  enddo
               enddo
            enddo
         enddo
      else if (ny .ge. nx .and. ny .ge. nz) then
         do l = 1, ncomp
            do k = kmin,kmax
               do i = imin,imax
                  do j = jmin,jmax
                     dest(i,j,k,l) = src(i-ioff,j-joff,k-koff,l)
                  enddo
               enddo
            enddo
         enddo
      else 
         do l = 1, ncomp
            do j = jmin,jmax
               do i = imin,imax
                  do k = kmin,kmax
                     dest(i,j,k,l) = src(i-ioff,j-joff,k-koff,l)
                  enddo
               enddo
            enddo
         enddo
      endif

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
c     Local variables.
c
      integer i,j,k,l

      integer nx,ny,nz
      integer minnx
      parameter (minnx = 64)

      nx = hi(1)-lo(1)+1
      ny = hi(2)-lo(2)+1
      nz = hi(3)-lo(3)+1

      if ((nx .ge. ny .and. nx .ge. nz) .or. (nx .ge. minnx)) then
         do l = 1, ncomp
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     dest(i,j,k,l) = val
                  enddo
               enddo
            enddo
         enddo
      else if (ny .ge. nx .and. ny .ge. nz) then
         do l = 1, ncomp
            do k = lo(3),hi(3)
               do i = lo(1),hi(1)
                  do j = lo(2),hi(2)
                     dest(i,j,k,l) = val
                  enddo
               enddo
            enddo
         enddo
      else
         do l = 1, ncomp
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  do k = lo(3),hi(3)
                     dest(i,j,k,l) = val
                  enddo
               enddo
            enddo
         enddo
      endif

      end

#endif /*BL_USE_SPECIALIZE*/
