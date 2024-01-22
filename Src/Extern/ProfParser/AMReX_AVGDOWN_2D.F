#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_AVGDOWN_F.H"
#include "AMReX_ArrayLim.H"

#ifdef BL_USE_FLOAT
#define REAL_T  real
#define zero    0.0e0
#define one     1.0e0
#else
#define REAL_T  double precision
#define zero    0.0d0
#define one     1.0d0
#endif

#define SDIM 2

      subroutine FORT_MAXVAL_AVGDOWN (
     &     crse,crse_l1,crse_l2,crse_h1,crse_h2,nvar,
     &     fine,fine_l1,fine_l2,fine_h1,fine_h2,
     &     lo,hi,ratios)
c     ----------------------------------------------------------
c     Average the fine grid data onto the coarse
c     grid.  Overlap is given in coarse grid coordinates.
c
c     crse      = coarse grid data
c     nvar        = number of components in arrays
c     fine      = fine grid data
c     lo,hi     = index limits of overlap (crse grid)
c     ratios    = IntVect refinement ratio
c     ----------------------------------------------------------
      implicit none
      integer  crse_l1,crse_l2,crse_h1,crse_h2
      integer  fine_l1,fine_l2,fine_h1,fine_h2
      integer  lo(SDIM), hi(SDIM)
      integer  nvar
      integer  ratios(SDIM)
      REAL_T   crse(crse_l1:crse_h1,crse_l2:crse_h2,nvar)
      REAL_T   fine(fine_l1:fine_h1,fine_l2:fine_h2,nvar)

      integer  i, j, n, ic, jc, ioff, joff
      integer  lratx,lraty

      lratx = ratios(1)
      lraty = ratios(2)

      do n = 1, nvar
c
c     set coarse grid to zero on overlap
c
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,n) = zero
            end do
         end do
c
c     sum fine data
c
         do joff = 0, lraty-1
            do jc = lo(2), hi(2)
               j = jc*lraty + joff
               do ioff = 0, lratx-1
                  do ic = lo(1), hi(1)
                     i = ic*lratx + ioff
                     crse(ic,jc,n) = max(crse(ic,jc,n), fine(i,j,n))
                  end do
               end do

            end do
         end do

      end do

      end


      subroutine FORT_CV_AVGDOWN (
     &     crse,crse_l1,crse_l2,crse_h1,crse_h2,nvar,
     &     fine,fine_l1,fine_l2,fine_h1,fine_h2,
     &     lo,hi,ratios)
c     ----------------------------------------------------------
c     Average the fine grid data onto the coarse
c     grid.  Overlap is given in coarse grid coordinates.
c
c     crse      = coarse grid data
c     nvar        = number of components in arrays
c     fine      = fine grid data
c     lo,hi     = index limits of overlap (crse grid)
c     ratios    = IntVect refinement ratio
c     ----------------------------------------------------------
      implicit none
      integer  crse_l1,crse_l2,crse_h1,crse_h2
      integer  fine_l1,fine_l2,fine_h1,fine_h2
      integer  lo(SDIM), hi(SDIM)
      integer  nvar
      integer  ratios(SDIM)
      REAL_T   crse(crse_l1:crse_h1,crse_l2:crse_h2,nvar)
      REAL_T   fine(fine_l1:fine_h1,fine_l2:fine_h2,nvar)

      integer  i, j, n, ic, jc, ioff, joff
      integer  lratx,lraty
      REAL_T   vol_inv

      lratx = ratios(1)
      lraty = ratios(2)
      vol_inv = one / (lratx * lraty)

      do n = 1, nvar
c
c     set coarse grid to zero on overlap
c
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,n) = zero
            end do
         end do
c
c     sum fine data
c
         do joff = 0, lraty-1
            do jc = lo(2), hi(2)
               j = jc*lraty + joff
               do ioff = 0, lratx-1
                  do ic = lo(1), hi(1)
                     i = ic*lratx + ioff
                     crse(ic,jc,n) = crse(ic,jc,n) + fine(i,j,n)
                  end do
               end do

            end do
         end do

         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,n) = crse(ic,jc,n) * vol_inv
            end do
         end do

      end do

      end



      subroutine FORT_CV_AVGDOWN_STAG (
     &     nodal_dir,
     &     crse,crse_l1,crse_l2,crse_h1,crse_h2,nvar,
     &     fine,fine_l1,fine_l2,fine_h1,fine_h2,
     &     lo,hi,ratios)
c     ----------------------------------------------------------
c     Average the fine grid data onto the coarse
c     grid.  Overlap is given in coarse grid coordinates.
c
c     crse      = coarse grid data
c     nvar        = number of components in arrays
c     fine      = fine grid data
c     lo,hi     = index limits of overlap (crse grid)
c     ratios    = IntVect refinement ratio
c     ----------------------------------------------------------
      implicit none
      integer  nodal_dir
      integer  crse_l1,crse_l2,crse_h1,crse_h2
      integer  fine_l1,fine_l2,fine_h1,fine_h2
      integer  lo(SDIM), hi(SDIM)
      integer  nvar
      integer  ratios(SDIM)
      REAL_T   crse(crse_l1:crse_h1,crse_l2:crse_h2,nvar)
      REAL_T   fine(fine_l1:fine_h1,fine_l2:fine_h2,nvar)

      integer  i, j, n, ic, jc, ioff, joff
      integer  lrat
      REAL_T   vol_inv

      if (ratios(1) .ne. ratios(2)) then
         print*,'Error: expecting same refinement ratio in each dir'
         stop
      end if

c     NOTE: switch from C++ 0-based indexing
      lrat = ratios(nodal_dir+1)

      vol_inv = 1.d0 / dble(lrat)

      do n = 1, nvar

c
c     set coarse grid to zero on overlap
c     NOTE: lo and hi already carries the +1 indexing for nodal, so no need to change this
c
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,n) = zero
            end do
         end do

c
c     sum fine data
c

         if (nodal_dir .eq. 0) then

            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  do joff=0,lrat-1
                     crse(i,j,n) = crse(i,j,n) +
     $                    vol_inv*fine(lrat*i,lrat*j+joff,n)
                  end do
               end do
            end do

         else

            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  do ioff=0,lrat-1
                     crse(i,j,n) = crse(i,j,n) +
     $                    vol_inv*fine(lrat*i+ioff,lrat*j,n)
                  end do
               end do
            end do

         end if

      end do

      end

      subroutine FORT_AVGDOWN (
     &     crse,crse_l1,crse_l2,crse_h1,crse_h2,nvar,
     &     fine,fine_l1,fine_l2,fine_h1,fine_h2,
     &     cv,cv_l1,cv_l2,cv_h1,cv_h2,
     &     fv,fv_l1,fv_l2,fv_h1,fv_h2,
     &     lo,hi,ratios)
c     ----------------------------------------------------------
c     Volume-weight average the fine grid data onto the coarse
c     grid.  Overlap is given in coarse grid coordinates.
c
c     crse      =  coarse grid data
c     nvar        = number of components in arrays
c     fine      = fine grid data
c     cv        = coarse grid volume array
c     fv        = fine grid volume array
c     lo,hi     = index limits of overlap (crse grid)
c     ratios    = IntVect refinement ratio
c     ----------------------------------------------------------
      implicit none
      integer  crse_l1,crse_l2,crse_h1,crse_h2
      integer  cv_l1,cv_l2,cv_h1,cv_h2
      integer  fine_l1,fine_l2,fine_h1,fine_h2
      integer  fv_l1,fv_l2,fv_h1,fv_h2
      integer  lo(SDIM), hi(SDIM)
      integer  nvar
      integer  ratios(SDIM)
      REAL_T   crse(crse_l1:crse_h1,crse_l2:crse_h2,nvar)
      REAL_T     cv(cv_l1:cv_h1,cv_l2:cv_h2)
      REAL_T   fine(fine_l1:fine_h1,fine_l2:fine_h2,nvar)
      REAL_T     fv(fv_l1:fv_h1,fv_l2:fv_h2)

      integer  i, j, n, ic, jc, ioff, joff
      integer  lratx,lraty

      lratx = ratios(1)
      lraty = ratios(2)

      do n = 1, nvar
c
c     set coarse grid to zero on overlap
c
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,n) = zero
            end do
         end do
c
c     sum fine data
c
         do joff = 0, lraty-1
            do jc = lo(2), hi(2)
               j = jc*lraty + joff
               do ioff = 0, lratx-1
                  do ic = lo(1), hi(1)
                     i = ic*lratx + ioff
                     crse(ic,jc,n) = crse(ic,jc,n) +
     &                    fv(i,j)*fine(i,j,n)
                  end do
               end do
            end do
         end do
c
c     divide out by volume weight
c
         do ic = lo(1), hi(1)
            do jc = lo(2), hi(2)
               crse(ic,jc,n) = crse(ic,jc,n)/cv(ic,jc)
            end do
         end do
      end do

      end

