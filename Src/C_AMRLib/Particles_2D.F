! ::
! :: ----------------------------------------------------------
! :: Sum the fine grid mass into the coarse grid.
! :: Overlap is given in coarse grid coordinates.
! :: ----------------------------------------------------------
! ::
      subroutine part_summassdown(crse,c_l1,c_l2,c_h1,c_h2, &
                                  nvar, &
                                  fine,f_l1,f_l2,f_h1,f_h2, &
                                  lo,hi,lrat)
      implicit none

      integer          c_l1,c_l2,c_h1,c_h2
      integer          f_l1,f_l2,f_h1,f_h2
      integer          lo(2), hi(2)
      integer          nvar, lrat(2)
      double precision crse(c_l1:c_h1,c_l2:c_h2,nvar)
      double precision fine(f_l1:f_h1,f_l2:f_h2,nvar)

      integer          i, j, n, ic, jc, ioff, joff
      integer          lratx, lraty, lratz

      lratx = lrat(1)
      lraty = lrat(2)

      do n = 1, nvar

         do joff = 0, lraty-1
            do jc = lo(2), hi(2)
               j = jc*lraty + joff
               do ioff = 0, lratx-1
                  do ic = lo(1), hi(1)
                     i = ic*lratx + ioff
                     crse(ic,jc,n) = crse(ic,jc,n) + fine(i,j,n)
                  enddo
               enddo
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
      subroutine part_summassup(crse,c_l1,c_l2,c_h1,c_h2, &
                                nvar, &
                                fine,f_l1,f_l2,f_h1,f_h2, &
                                lo,hi,lrat)
      implicit none

      integer          c_l1,c_l2,c_h1,c_h2
      integer          f_l1,f_l2,f_h1,f_h2
      integer          lo(2), hi(2)
      integer          nvar, lrat(2)
      double precision crse(c_l1:c_h1,c_l2:c_h2,nvar)
      double precision fine(f_l1:f_h1,f_l2:f_h2,nvar)

      integer          i, j, n, ic, jc, ioff, joff
      integer          lratx, lraty
      double precision volfrac

      lratx = lrat(1)
      lraty = lrat(2)

      volfrac = 1.d0/float(lratx*lraty)

      do n = 1, nvar

         do joff = 0, lraty-1
            do jc = lo(2), hi(2)
               j = jc*lraty + joff
               do ioff = 0, lratx-1
                  do ic = lo(1), hi(1)
                     i = ic*lratx + ioff
                     fine(i,j,n) = fine(i,j,n) + volfrac*crse(ic,jc,n)
                  enddo
               enddo
            enddo
         enddo

      enddo

      end subroutine part_summassup
