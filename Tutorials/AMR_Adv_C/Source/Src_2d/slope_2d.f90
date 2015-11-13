module slope_module
 
  implicit none
 
  private
 
  public slope
 
contains
 
! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine slope(q,qd_l1,qd_l2,qd_h1,qd_h2, &
                       dq,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                       ilo1,ilo2,ihi1,ihi2,nv,idir)

      implicit none

      integer ilo,ihi
      integer qd_l1,qd_l2,qd_h1,qd_h2
      integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
      integer ilo1,ilo2,ihi1,ihi2,nv,idir

      double precision     q( qd_l1: qd_h1, qd_l2: qd_h2,nv)
      double precision    dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,nv)

!     Local arrays
      double precision, allocatable::dsgn(:),dlim(:),df(:),dcen(:)

      integer i, j, n
      double precision dlft, drgt, slop, dq1
      double precision four3rd, sixth

      four3rd = 4.d0/3.d0
      sixth = 1.d0/6.d0

      ilo = MIN(ilo1,ilo2)
      ihi = MAX(ihi1,ihi2)

      allocate (dsgn(ilo-2:ihi+2))
      allocate (dlim(ilo-2:ihi+2))
      allocate (  df(ilo-2:ihi+2))
      allocate (dcen(ilo-2:ihi+2))

      do n = 1, nv 
          if (idir .eq. 1) then

             ! slopes in first coordinate direction
             do j = ilo2-1, ihi2+1

                ! first compute Fromm slopes
                do i = ilo1-2, ihi1+2
                      dlft = 2.d0*(q(i  ,j,n) - q(i-1,j,n))
                      drgt = 2.d0*(q(i+1,j,n) - q(i  ,j,n))
                      dcen(i) = .25d0 * (dlft+drgt)
                      dsgn(i) = sign(1.d0, dcen(i))
                      slop = min( abs(dlft), abs(drgt) )
!                      dlim(i) = cvmgp( slop, 0.d0, dlft*drgt )
                      if (dlft*drgt .ge. 0.d0) then
                         dlim(i) = slop
                      else
                         dlim(i) = 0.d0
                      endif
                      df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
                  enddo

!                 Now limited fourth order slopes
                  do i = ilo1-1, ihi1+1
                      dq1 = four3rd*dcen(i) - sixth*(df(i+1) + df(i-1))
                      dq(i,j,n) = dsgn(i)*min(dlim(i),abs(dq1))
                  enddo
              enddo

          else


!            Compute slopes in second coordinate direction
             do i = ilo1-1, ihi1+1

!               First compute Fromm slopes for this column
                do j = ilo2-2, ihi2+2
                      dlft = 2.d0*(q(i,j  ,n) - q(i,j-1,n))
                      drgt = 2.d0*(q(i,j+1,n) - q(i,j  ,n))
                      dcen(j) = .25d0 * (dlft+drgt)
                      dsgn(j) = sign( 1.d0, dcen(j) )
                      slop = min( abs(dlft), abs(drgt) )
                      if (dlft*drgt .ge. 0.d0) then
                         dlim(j) = slop
                      else
                         dlim(j) = 0.d0
                      endif
                      df(j) = dsgn(j)*min( dlim(j),abs(dcen(j)) )
                  enddo

!                 Now compute limited fourth order slopes
                  do j = ilo2-1, ihi2+1
                      dq1 = four3rd*dcen(j) - &
                           sixth*( df(j+1) + df(j-1) )
                      dq(i,j,n) = dsgn(j)*min(dlim(j),abs(dq1))
                  enddo
              enddo

          endif
      enddo

      deallocate(dsgn,dlim,df,dcen)

      end subroutine slope

end module slope_module 
