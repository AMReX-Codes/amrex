module slope_module

  use amrex_base_module
 
  implicit none

  real(amrex_real), parameter:: four3rd=4._amrex_real/3._amrex_real, &
       sixth=1._amrex_real/6._amrex_real
  
  private
 
  public :: slopex, slopey
 
contains
 
  subroutine slopex(lo, hi, &
                    q, qlo, qhi, &
                    dq, dqlo, dqhi)

    implicit none

    integer, intent(in) :: lo(2), hi(2), qlo(2), qhi(2), dqlo(2), dqhi(2)
    real(amrex_real), intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2))
    real(amrex_real), intent(out) :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2))

    integer :: i, j
    real(amrex_real), dimension(lo(1)-1:hi(1)+1) :: dsgn, dlim, df, dcen
    real(amrex_real) :: dlft, drgt, dq1

    do j = lo(2), hi(2)

       ! first compute Fromm slopes
       do i = lo(1)-1, hi(1)+1
          dlft = q(i  ,j) - q(i-1,j)
          drgt = q(i+1,j) - q(i  ,j)
          dcen(i) = .5d0 * (dlft+drgt)
          dsgn(i) = sign(1.d0, dcen(i))
          if (dlft*drgt .ge. 0.d0) then
             dlim(i) = 2.d0 * min( abs(dlft), abs(drgt) )
          else
             dlim(i) = 0.d0
          endif
          df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
       end do

       ! Now limited fourth order slopes
       do i = lo(1), hi(1)
          dq1 = four3rd*dcen(i) - sixth*(df(i+1) + df(i-1))
          dq(i,j) = dsgn(i)*min(dlim(i),abs(dq1))
       end do
    enddo

  end subroutine slopex


  subroutine slopey(lo, hi, &
                    q, qlo, qhi, &
                    dq, dqlo, dqhi)

    integer, intent(in) :: lo(2), hi(2), qlo(2), qhi(2), dqlo(2), dqhi(2)
    real(amrex_real), intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2))
    real(amrex_real), intent(out) :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2))

    real(amrex_real), dimension(:,:), pointer, contiguous :: dsgn, dlim, df, dcen

    call amrex_allocate(dsgn, lo(1), hi(1), lo(2)-1, hi(2)+1)
    call amrex_allocate(dlim, lo(1), hi(1), lo(2)-1, hi(2)+1)
    call amrex_allocate(df  , lo(1), hi(1), lo(2)-1, hi(2)+1)
    call amrex_allocate(dcen, lo(1), hi(1), lo(2)-1, hi(2)+1)

    call slopey_doit(lo, hi, &
                     q, qlo, qhi, &
                     dq, dqlo, dqhi, &
                     dsgn, dlim, df, dcen, (/lo(1),lo(2)-1/), (/hi(1),hi(2)+1/))

    call amrex_deallocate(dsgn)
    call amrex_deallocate(dlim)
    call amrex_deallocate(df)
    call amrex_deallocate(dcen)

  end subroutine slopey

  subroutine slopey_doit(lo, hi, &
                         q, qlo, qhi, &
                         dq, dqlo, dqhi, &
                         dsgn, dlim, df, dcen, ddlo, ddhi)

    integer, intent(in) :: lo(2), hi(2), qlo(2), qhi(2), dqlo(2), dqhi(2), &
         ddlo(2), ddhi(2)
    real(amrex_real), intent(in ) ::  q  ( qlo(1): qhi(1), qlo(2): qhi(2))
    real(amrex_real), intent(out) :: dq  (dqlo(1):dqhi(1),dqlo(2):dqhi(2))
    real(amrex_real)              :: dsgn(ddlo(1):ddhi(1),ddlo(2):ddhi(2))
    real(amrex_real)              :: dlim(ddlo(1):ddhi(1),ddlo(2):ddhi(2))
    real(amrex_real)              :: df  (ddlo(1):ddhi(1),ddlo(2):ddhi(2))
    real(amrex_real)              :: dcen(ddlo(1):ddhi(1),ddlo(2):ddhi(2))

    integer :: i, j
    real(amrex_real) :: dlft, drgt, dq1

    ! first compute Fromm slopes
    do j    = lo(2)-1, hi(2)+1
       do i = lo(1)  , hi(1)
          dlft = q(i,j  ) - q(i,j-1)
          drgt = q(i,j+1) - q(i,j  )
          dcen(i,j) = .5d0 * (dlft+drgt)
          dsgn(i,j) = sign( 1.d0, dcen(i,j) )
          if (dlft*drgt .ge. 0.d0) then
             dlim(i,j) = 2.d0 * min( abs(dlft), abs(drgt) )
          else
             dlim(i,j) = 0.d0
          endif
          df(i,j) = dsgn(i,j)*min( dlim(i,j),abs(dcen(i,j)) )
       end do
    end do

    ! Now compute limited fourth order slopes
    do j    = lo(2), hi(2)
       do i = lo(1), hi(1)
          dq1 = four3rd*dcen(i,j) - sixth*( df(i,j+1) + df(i,j-1) )
          dq(i,j) = dsgn(i,j)*min(dlim(i,j),abs(dq1))
       end do
    end do

  end subroutine slopey_doit

end module slope_module 
