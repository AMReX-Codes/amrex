module slope_module
  
  implicit none

  private

  public slope

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine slope( q, qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3, &
                       dq,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                       ilo,jlo,klo,ihi,jhi,khi,nv,idir)

      use meth_params_module
      use bl_constants_module

      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo,jlo,klo,ihi,jhi,khi,nv,idir

      double precision q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,nv)
      double precision dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,nv)

      integer i,j,k,n
      integer lo,hi

      double precision dlft, drgt, slop, dq1
      double precision dm, dp, dc, ds, sl, dl, dfm, dfp

      double precision, allocatable::dsgn(:),dlim(:),df(:),dcen(:)

      if (idir .eq. 1) then
          lo = ilo
          hi = ihi
      else if (idir .eq. 2) then
          lo = jlo
          hi = jhi
      else if (idir .eq. 3) then
          lo = klo
          hi = khi
      end if

      allocate(dsgn(lo-2:hi+2))
      allocate(dlim(lo-2:hi+2))
      allocate(  df(lo-2:hi+2))
      allocate(dcen(lo-2:hi+2))

      do n = 1, nv 

          ! Compute slopes in first coordinate direction
          if (idir .eq. 1) then
            do k = klo-1, khi+1
            do j = jlo-1, jhi+1

               ! First compute Fromm slopes
               do i = ilo-2, ihi+2
                  dlft = TWO*(q(i  ,j,k,n) - q(i-1,j,k,n))
                  drgt = TWO*(q(i+1,j,k,n) - q(i  ,j,k,n))
                  dcen(i) = FOURTH * (dlft+drgt)
                  dsgn(i) = sign(ONE, dcen(i))
                  slop = min( abs(dlft), abs(drgt) )
                  if (dlft*drgt .ge. ZERO) then
                     dlim(i) = slop
                  else
                     dlim(i) = ZERO
                  endif
                  df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
               enddo

               ! Now compute limited fourth order slopes
               do i = ilo-1, ihi+1
                  dq1         = FOUR3RD*dcen(i) - SIXTH*(df(i+1) + df(i-1))
                  dq(i,j,k,n) = dsgn(i)*min(dlim(i),abs(dq1))
               enddo

            enddo
            enddo

          else if (idir .eq. 2) then

            ! Compute slopes in second coordinate direction
            do k = klo-1, khi+1
            do i = ilo-1, ihi+1

               ! First compute Fromm slopes for this column
               do j = jlo-2, jhi+2
                  dlft = TWO*(q(i,j  ,k,n) - q(i,j-1,k,n))
                  drgt = TWO*(q(i,j+1,k,n) - q(i,j  ,k,n))
                  dcen(j) = FOURTH * (dlft+drgt)
                  dsgn(j) = sign( ONE, dcen(j) )
                  slop = min( abs(dlft), abs(drgt) )
                  if (dlft*drgt .ge. ZERO) then
                     dlim(j) = slop
                  else
                     dlim(j) = ZERO
                  endif
                  df(j) = dsgn(j)*min( dlim(j),abs(dcen(j)) )
               enddo

               ! Now compute limited fourth order slopes
               do j = jlo-1, jhi+1
                  dq1 = FOUR3RD*dcen(j) - SIXTH*( df(j+1) + df(j-1) )
                  dq(i,j,k,n) = dsgn(j)*min(dlim(j),abs(dq1))
               enddo
            enddo
            enddo

          else if (idir .eq. 3) then

            ! Compute slopes in third coordinate direction
            do j = jlo-1, jhi+1
            do i = ilo-1, ihi+1

               ! First compute Fromm slopes for this column
               do k = klo-2, khi+2
                  dlft = TWO*(q(i,j,k  ,n) - q(i,j,k-1,n))
                  drgt = TWO*(q(i,j,k+1,n) - q(i,j,k  ,n))
                  dcen(k) = FOURTH * (dlft+drgt)
                  dsgn(k) = sign( ONE, dcen(k) )
                  slop = min( abs(dlft), abs(drgt) )
                  if (dlft*drgt .ge. ZERO) then
                     dlim(k) = slop
                  else
                     dlim(k) = ZERO
                  endif
                  df(k) = dsgn(k)*min( dlim(k),abs(dcen(k)) )
               enddo

               ! Now compute limited fourth order slopes
               do k = klo-1, khi+1
                  dq1 = FOUR3RD*dcen(k) - SIXTH*( df(k+1) + df(k-1) )
                  dq(i,j,k,n) = dsgn(k)*min(dlim(k),abs(dq1))
               enddo
            enddo
            enddo

          endif ! idir

      end do  ! n = 1, nv

      deallocate(dsgn,dlim,df,dcen)

      end subroutine slope

end module slope_module
