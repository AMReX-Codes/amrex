#if 0
! ::
! :: ----------------------------------------------------------
! :: Sum the fine grid mass into the coarse grid.
! :: Overlap is given in coarse grid coordinates.
! :: ----------------------------------------------------------
! ::
      subroutine part_summassdown(crse,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,
     &                            fine,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3,
     &                            lo,hi,lrat)
      implicit none

      integer          c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
      integer          f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
      integer          lo(3), hi(3), lrat(3)
      double precision crse(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3)
      double precision fine(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3)

      integer          i, j, k, ic, jc, kc, ioff, joff, koff
      integer          lratx, lraty, lratz

      lratx = lrat(1)
      lraty = lrat(2)
      lratz = lrat(3)

      do koff = 0, lratz-1
!$OMP PARALLEL DO PRIVATE(i,j,k,ic,jc,kc,ioff,joff)
         do kc = lo(3),hi(3)
            k = kc*lratz + koff
            do joff = 0, lraty-1
               do jc = lo(2), hi(2)
                  j = jc*lraty + joff
                  do ioff = 0, lratx-1
                     do ic = lo(1), hi(1)
                        i = ic*lratx + ioff
                        crse(ic,jc,kc) = crse(ic,jc,kc) + fine(i,j,k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
      enddo

      end subroutine part_summassdown
! ::
! :: ----------------------------------------------------------
! :: Sum the crse grid mass into the fine grid.
! :: Overlap is given in coarse grid coordinates.
! :: ----------------------------------------------------------
! ::
      subroutine part_summassup(crse,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,
     &                          fine,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3,
     &                          lo,hi,lrat)
      implicit none

      integer          c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
      integer          f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
      integer          lo(3), hi(3), lrat(3)
      double precision crse(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3)
      double precision fine(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3)

      integer          i, j, k, ic, jc, kc, ioff, joff, koff
      integer          lratx, lraty, lratz
      double precision volfrac

      lratx = lrat(1)
      lraty = lrat(2)
      lratz = lrat(3)

      volfrac = 1.d0/float(lratx*lraty*lratz)

      do koff = 0, lratz-1
!$OMP PARALLEL DO PRIVATE(i,j,k,ic,jc,kc,ioff,joff)
         do kc = lo(3),hi(3)
            k = kc*lratz + koff
            do joff = 0, lraty-1
               do jc = lo(2), hi(2)
                  j = jc*lraty + joff
                  do ioff = 0, lratx-1
                     do ic = lo(1), hi(1)
                        i = ic*lratx + ioff
                        fine(i,j,k) = fine(i,j,k) + volfrac*crse(ic,jc,kc)
                     enddo
                  enddo
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
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
      subroutine part_sumdensup(crse,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,
     &                          fine,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3,
     &                          lo,hi,lrat)
      implicit none

      integer          c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
      integer          f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
      integer          lo(3), hi(3), lrat(3)
      double precision crse(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3)
      double precision fine(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3)

      integer          i, j, k, ic, jc, kc, ioff, joff, koff
      integer          lratx, lraty, lratz

      lratx = lrat(1)
      lraty = lrat(2)
      lratz = lrat(3)

      do koff = 0, lratz-1
!$OMP PARALLEL DO PRIVATE(i,j,k,ic,jc,kc,ioff,joff)
         do kc = lo(3),hi(3)
            k = kc*lratz + koff
            do joff = 0, lraty-1
               do jc = lo(2), hi(2)
                  j = jc*lraty + joff
                  do ioff = 0, lratx-1
                     do ic = lo(1), hi(1)
                        i = ic*lratx + ioff
                        fine(i,j,k) = fine(i,j,k) + crse(ic,jc,kc)
                     enddo
                  enddo
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
      enddo

      end subroutine part_sumdensup

#endif
