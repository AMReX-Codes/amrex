! :: ----------------------------------------------------------
! :: Average the fine grid phi onto the coarse
! :: grid.  Overlap is given in coarse grid coordinates.
! :: Note this differs from fort_avgdown in that there is no volume weighting.
! ::
! :: INPUTS / OUTPUTS:
! ::  crse      <=  coarse grid data
! ::  clo,chi    => index limits of crse array interior
! ::  fine       => fine grid data
! ::  flo,fhi    => index limits of fine array interior
! ::  rfine      => (ignore) used in 2-D RZ calc
! ::  lo,hi      => index limits of overlap (crse grid)
! ::  lrat       => refinement ratio
! ::
! :: NOTE:
! ::  Assumes all data cell centered
! :: ----------------------------------------------------------
! ::
      subroutine fort_avgdown_phi (crse,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3, &
                                   fine,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                                   lo,hi,lrat)
      use amrex_fort_module, only : rt => amrex_real
      implicit none
      integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
      integer f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
      integer lo(3), hi(3)
      integer lrat(3)
      real(rt) crse(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3)
      real(rt) fine(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3)

      integer i, j, k, ic, jc, kc, ioff, joff, koff
      integer lratx, lraty, lratz
      real(rt) volfrac

      lratx   = lrat(1)
      lraty   = lrat(2)
      lratz   = lrat(3)
      volfrac = 1.d0/float(lrat(1)*lrat(2)*lrat(3))
      !
      ! ::::: set coarse grid to zero on overlap
      !
      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,kc) = 0.d0
            enddo
         enddo
      enddo
      !
      ! ::::: sum fine data
      !
      do koff = 0, lratz-1
        do kc = lo(3), hi(3)
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
      enddo

      do kc = lo(3), hi(3)
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,kc) = volfrac*crse(ic,jc,kc)
            enddo
         enddo
      enddo

      end subroutine fort_avgdown_phi

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine fort_edge_interp(flo, fhi, nc, ratio, dir, &
           fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)

      use amrex_fort_module, only : rt => amrex_real
      implicit none
      integer flo(0:3-1), fhi(0:3-1), nc, ratio(0:3-1), dir
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      real(rt) &
           fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,n,P,M,L
      real(rt) val, df
      !
      ! Do linear in dir, pc transverse to dir, leave alone the fine values
      ! lining up with coarse edges--assume these have been set to hold the 
      ! values you want to interpolate to the rest.
      !
      if (dir.eq.0) then
         do n=1,nc
            do k=flo(2),fhi(2),ratio(2)
               do j=flo(1),fhi(1),ratio(1)
                  do i=flo(0),fhi(0)-ratio(dir),ratio(0)
                     df = fine(i+ratio(dir),j,k,n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n) &
                             + df*dble(M)/dble(ratio(dir))
                        do P=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                           do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                              fine(i+M,P,L,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do k=flo(2),fhi(2),ratio(2)
               do j=flo(1),fhi(1)-ratio(dir),ratio(1)
                  do i=flo(0),fhi(0)
                     df = fine(i,j+ratio(dir),k,n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n) &
                             + df*dble(M)/dble(ratio(dir))
                        do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                           do L=MAX(k,flo(2)),MIN(k+ratio(2)-1,fhi(2))
                              fine(P,j+M,L,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do k=flo(2),fhi(2)-ratio(dir),ratio(2)
               do j=flo(1),fhi(1),ratio(1)
                  do i=flo(0),fhi(0),ratio(0)
                     df = fine(i,j,k+ratio(dir),n)-fine(i,j,k,n)
                     do M=1,ratio(dir)-1
                        val = fine(i,j,k,n) &
                             + df*dble(M)/dble(ratio(dir))
                        do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                           do L=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                              fine(P,L,k+M,n) = val
                           enddo
                        enddo
                     enddo                     
                  enddo
               enddo
            enddo
         enddo
      endif

      end subroutine fort_edge_interp

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine fort_pc_edge_interp(lo, hi, nc, ratio, dir, &
           crse, crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2, &
           fine, fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2)

      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer lo(3),hi(3), nc, ratio(0:3-1), dir
      integer crse_l0, crse_l1, crse_l2, crse_h0, crse_h1, crse_h2
      integer fine_l0, fine_l1, fine_l2, fine_h0, fine_h1, fine_h2
      real(rt) &
           crse(crse_l0:crse_h0,crse_l1:crse_h1,crse_l2:crse_h2,nc)
      real(rt) &
           fine(fine_l0:fine_h0,fine_l1:fine_h1,fine_l2:fine_h2,nc)
      integer i,j,k,ii,jj,kk,n,L, P
      !
      ! For edge-based data, fill fine values with piecewise-constant interp of coarse data.
      ! Operate only on faces that overlap--ie, only fill the fine faces that make up each
      ! coarse face, leave the in-between faces alone.
      !
      if (dir.eq.0) then
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(2)-1
                        do L=0,ratio(1)-1
                           fine(ii,jj+L,kk+P,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(2)-1
                        do L=0,ratio(0)-1
                           fine(ii+L,jj,kk+P,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do k=lo(3),hi(3)
               kk = ratio(2)*k
               do j=lo(2),hi(2)
                  jj = ratio(1)*j
                  do i=lo(1),hi(1)
                     ii = ratio(0)*i
                     do P=0,ratio(1)-1
                        do L=0,ratio(0)-1
                           fine(ii+L,jj+P,kk,n) = crse(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

      end subroutine fort_pc_edge_interp

