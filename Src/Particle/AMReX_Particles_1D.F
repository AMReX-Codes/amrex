#if 0
! ::
! :: ----------------------------------------------------------
! :: Sum the fine grid mass into the coarse grid.
! :: Overlap is given in coarse grid coordinates.
! :: ----------------------------------------------------------
! ::
      subroutine part_summassdown(crse,c_l1,c_h1,
     &                            fine,f_l1,f_h1,
     &                            lo,hi,lrat)
      implicit none

      integer          c_l1,c_h1
      integer          f_l1,f_h1
      integer          lo(1), hi(1), lrat(1)
      double precision crse(c_l1:c_h1)
      double precision fine(f_l1:f_h1)

      integer          i, ic, ioff, lratx

      lratx = lrat(1)

      do ioff = 0, lratx-1
         do ic = lo(1), hi(1)
            i = ic*lratx + ioff
            crse(ic) = crse(ic) + fine(i)
         enddo
      enddo

      end subroutine part_summassdown
! ::
! :: ----------------------------------------------------------
! :: Sum the crse grid mass into the fine grid.
! :: Overlap is given in coarse grid coordinates.
! :: ----------------------------------------------------------
! ::
      subroutine part_summassup(crse,c_l1,c_h1,
     &                          fine,f_l1,f_h1,
     &                          lo,hi,lrat)
      implicit none

      integer          c_l1,c_h1
      integer          f_l1,f_h1
      integer          lo(1), hi(1), lrat(1)
      double precision crse(c_l1:c_h1)
      double precision fine(f_l1:f_h1)

      integer          i, ic, ioff, lratx
      double precision volfrac

      lratx = lrat(1)

      volfrac = 1.d0/float(lratx)

      do ioff = 0, lratx-1
         do ic = lo(1), hi(1)
            i = ic*lratx + ioff
            fine(i) = fine(i) + volfrac*crse(ic)
         enddo
      enddo

      end subroutine part_summassup

! ::
! :: ----------------------------------------------------------
! :: Sum the crse grid density into the fine grid.
! :: Overlap is given in coarse grid coordinates.
! :: This differs from summassup in that the crse contribution
! :: is *not* weighted by the volume fraction, since this is
! :: density, not mass
! :: ----------------------------------------------------------
! ::
      subroutine part_sumdensup(crse,c_l1,c_h1,
     &                          fine,f_l1,f_h1,
     &                          lo,hi,lrat)
      implicit none

      integer          c_l1,c_h1
      integer          f_l1,f_h1
      integer          lo(1), hi(1), lrat(1)
      double precision crse(c_l1:c_h1)
      double precision fine(f_l1:f_h1)

      integer          i, ic, ioff, lratx

      lratx = lrat(1)

      do ioff = 0, lratx-1
         do ic = lo(1), hi(1)
            i = ic*lratx + ioff
            fine(i) = fine(i) + crse(ic)
         enddo
      enddo

      end subroutine part_sumdensup

#endif
