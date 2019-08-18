#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_AVGDOWN_F.H"
#include "AMReX_ArrayLim.H"

#undef AMREX_XSDK

#define SDIM 1

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
      integer  lratx
      REAL_T   vol_inv

      lratx = ratios(1)
      vol_inv = one / (lratx)

      do n = 1, nvar
c
c     set coarse grid to zero on overlap
c
         do ic = lo(1), hi(1)
            crse(ic,n) = zero
         end do
c
c     sum fine data
c
            do ioff = 0, lratx-1
               do ic = lo(1), hi(1)
                  i = ic*lratx + ioff
                  crse(ic,n) = crse(ic,n) + fine(i,n)
               end do
            end do

         do ic = lo(1), hi(1)
            crse(ic,n) = crse(ic,n) * vol_inv
         end do
         
      end do
         
      end

      subroutine FORT_CV_AVGDOWN_STAG (
     &     nodal_dir,
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
      integer  nodal_dir
      integer  DIMDEC(crse)
      integer  DIMDEC(fine)
      integer  lo(SDIM), hi(SDIM)
      integer  nvar
      integer  ratios(SDIM)
      REAL_T   crse(DIMV(crse),nvar)
      REAL_T   fine(DIMV(fine),nvar)

      integer  i, j, n, ic, jc, ioff, joff
      integer  lrat
      REAL_T   vol_inv


c     NOTE: switch from C++ 0-based indexing
      lrat = ratios(nodal_dir+1)

      vol_inv = 1.d0 / dble(lrat)

      do n = 1, nvar

c
c     set coarse grid to zero on overlap
c     NOTE: lo and hi already carries the +1 indexing for nodal, so no need to change this
c
         do ic = lo(1), hi(1)
            crse(ic,n) = zero
         end do

c
c     sum fine data
c

         if (nodal_dir .eq. 0) then

            do i=lo(1),hi(1)
                  crse(i,n) = crse(i,n) + vol_inv*fine(lrat*i,n)
            end do

         else

            do i=lo(1),hi(1)
               do ioff=0,lrat-1
                  crse(i,n) = crse(i,n) + vol_inv*fine(lrat*i+ioff,n)
               end do
            end do

         end if
         
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
      integer  lratx

      lratx = ratios(1)

      do n = 1, nvar
c
c     set coarse grid to zero on overlap
c
         do ic = lo(1), hi(1)
            crse(ic,n) = zero
         end do
c
c     sum fine data
c
            do ioff = 0, lratx-1
               do ic = lo(1), hi(1)
                  i = ic*lratx + ioff
                  crse(ic,n) = crse(ic,n) +
     &                    fv(i)*fine(i,n)
               end do
            end do
c
c     divide out by volume weight
c
         do ic = lo(1), hi(1)
               crse(ic,n) = crse(ic,n)/cv(ic)
         end do
      end do

      end

