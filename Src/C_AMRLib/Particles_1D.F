! ::
! :: ----------------------------------------------------------
! :: Sum the fine grid mass into the coarse grid.
! :: Overlap is given in coarse grid coordinates.
! :: ----------------------------------------------------------
! ::
      subroutine part_summassdown(crse,c_l1,c_h1, &
                                  nvar, &
                                  fine,f_l1,f_h1, &
                                  lo,hi,lrat)
      implicit none

      integer          c_l1,c_h1
      integer          f_l1,f_h1
      integer          lo(1), hi(1)
      integer          nvar, lrat(1)
      double precision crse(c_l1:c_h1,nvar)
      double precision fine(f_l1:f_h1,nvar)

      integer          i, n, ic, ioff
      integer          lratx, lraty, lratz

      lratx = lrat(1)

      do n = 1, nvar

         do ioff = 0, lratx-1
            do ic = lo(1), hi(1)
               i = ic*lratx + ioff
               crse(ic,n) = crse(ic,n) + fine(i,n)
            enddo
         enddo

      enddo

      end subroutine part_summassdown
! ::
! :: ----------------------------------------------------------
! :: Sum the crse grid mass into the fine grid.
! :: Overlap is given in coarse grid coordinates.
! :: ----------------------------------------------------------
! ::
      subroutine part_summassup(crse,c_l1,c_h1, &
                                nvar, &
                                fine,f_l1,f_h1, &
                                lo,hi,lrat)
      implicit none

      integer          c_l1,c_h1
      integer          f_l1,f_h1
      integer          lo(1), hi(1)
      integer          nvar, lrat(1)
      double precision crse(c_l1:c_h1,nvar)
      double precision fine(f_l1:f_h1,nvar)

      integer          i, n, ic, ioff
      integer          lratx, lraty, lratz
      double precision volfrac

      lratx = lrat(1)

      volfrac = 1.d0/float(lratx)

      do n = 1, nvar

         do ioff = 0, lratx-1
            do ic = lo(1), hi(1)
               i = ic*lratx + ioff
               fine(i,j,k,n) = fine(i,j,k,n) + volfrac*crse(ic,jc,kc,n)
            enddo
         enddo

      enddo

      end subroutine part_summassup
