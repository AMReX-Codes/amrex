#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "AVGDOWN_F.H"
#include "ArrayLim.H"

#define SDIM 2

      subroutine FORT_CV_AVGDOWN (
     &     crse,DIMS(crse),nvar,
     &     fine,DIMS(fine),
     &     lo,hi,ratios)
c     ----------------------------------------------------------
c     Average the fine grid data onto the coarse
c     grid.  Overlap is given in coarse grid coordinates.
c     
c     crse      = coarse grid data
c     nvar	= number of components in arrays
c     fine      = fine grid data
c     lo,hi     = index limits of overlap (crse grid)
c     ratios    = IntVect refinement ratio
c     ----------------------------------------------------------
      integer  DIMDEC(crse)
      integer  DIMDEC(fine)
      integer  lo(SDIM), hi(SDIM)
      integer  nvar
      integer  ratios(SDIM)
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T   fine(DIMV(fine),nvar)

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

      subroutine FORT_AVGDOWN (
     &     crse,DIMS(crse),nvar,
     &     fine,DIMS(fine),
     &     cv,DIMS(cv),
     &     fv,DIMS(fv),
     &     lo,hi,ratios)
c     ----------------------------------------------------------
c     Volume-weight average the fine grid data onto the coarse
c     grid.  Overlap is given in coarse grid coordinates.
c     
c     crse      =  coarse grid data
c     nvar	= number of components in arrays
c     fine      = fine grid data
c     cv        = coarse grid volume array
c     fv        = fine grid volume array
c     lo,hi     = index limits of overlap (crse grid)
c     ratios    = IntVect refinement ratio
c     ----------------------------------------------------------
      integer  DIMDEC(crse)
      integer  DIMDEC(cv)
      integer  DIMDEC(fine)
      integer  DIMDEC(fv)
      integer  lo(SDIM), hi(SDIM)
      integer  nvar
      integer  ratios(SDIM)
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T     cv(DIMV(cv))
      REAL_T   fine(DIMV(fine),nvar)
      REAL_T     fv(DIMV(fv))

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

