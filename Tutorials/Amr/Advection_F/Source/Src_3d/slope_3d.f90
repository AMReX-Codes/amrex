module slope_module
 
  implicit none

  double precision, parameter:: four3rd=4.d0/3.d0, sixth=1.d0/6.d0
  
  private
 
  public :: slopex, slopey, slopez
 
contains
 
  subroutine slopex(lo, hi, &
                    q, qlo, qhi, &
                    dq, dqlo, dqhi)

    implicit none

    integer, intent(in) :: lo(3), hi(3), qlo(3), qhi(3), dqlo(3), dqhi(3)
    double precision, intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision, intent(out) :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2),dqlo(3):dqhi(3))

    integer :: i, j, k
    double precision, dimension(lo(1)-1:hi(1)+1) :: dsgn, dlim, df, dcen
    double precision :: dlft, drgt, dq1

    do    k = lo(3), hi(3)
       do j = lo(2), hi(2)

          ! first compute Fromm slopes
          do i = lo(1)-1, hi(1)+1
             dlft = q(i  ,j,k) - q(i-1,j,k)
             drgt = q(i+1,j,k) - q(i  ,j,k)
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
             dq(i,j,k) = dsgn(i)*min(dlim(i),abs(dq1))
          end do
       end do
    end do

  end subroutine slopex


  subroutine slopey(lo, hi, &
                    q, qlo, qhi, &
                    dq, dqlo, dqhi)

    use amrex_mempool_module, only : bl_allocate, bl_deallocate

    integer, intent(in) :: lo(3), hi(3), qlo(3), qhi(3), dqlo(3), dqhi(3)
    double precision, intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision, intent(out) :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2),dqlo(3):dqhi(3))

    ! Some compiler may not support 'contiguous'.  Remove it in that case.
    double precision, dimension(:,:), pointer, contiguous :: dsgn, dlim, df, dcen

    call bl_allocate(dsgn, lo(1), hi(1), lo(2)-1, hi(2)+1)
    call bl_allocate(dlim, lo(1), hi(1), lo(2)-1, hi(2)+1)
    call bl_allocate(df  , lo(1), hi(1), lo(2)-1, hi(2)+1)
    call bl_allocate(dcen, lo(1), hi(1), lo(2)-1, hi(2)+1)

    call slopey_doit(lo, hi, &
                     q, qlo, qhi, &
                     dq, dqlo, dqhi, &
                     dsgn, dlim, df, dcen, (/lo(1),lo(2)-1/), (/hi(1),hi(2)+1/))

    call bl_deallocate(dsgn)
    call bl_deallocate(dlim)
    call bl_deallocate(df)
    call bl_deallocate(dcen)

  end subroutine slopey

  subroutine slopey_doit(lo, hi, &
                         q, qlo, qhi, &
                         dq, dqlo, dqhi, &
                         dsgn, dlim, df, dcen, ddlo, ddhi)

    integer, intent(in) :: lo(3), hi(3), qlo(3), qhi(3), dqlo(3), dqhi(3), &
         ddlo(2), ddhi(2)
    double precision, intent(in ) ::  q  ( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision, intent(out) :: dq  (dqlo(1):dqhi(1),dqlo(2):dqhi(2),dqlo(3):dqhi(3))
    double precision              :: dsgn(ddlo(1):ddhi(1),ddlo(2):ddhi(2))
    double precision              :: dlim(ddlo(1):ddhi(1),ddlo(2):ddhi(2))
    double precision              :: df  (ddlo(1):ddhi(1),ddlo(2):ddhi(2))
    double precision              :: dcen(ddlo(1):ddhi(1),ddlo(2):ddhi(2))

    integer :: i, j, k
    double precision :: dlft, drgt, dq1

    do k = lo(3), hi(3)

       ! first compute Fromm slopes
       do j    = lo(2)-1, hi(2)+1
          do i = lo(1)  , hi(1)
             dlft = q(i,j  ,k) - q(i,j-1,k)
             drgt = q(i,j+1,k) - q(i,j  ,k)
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
             dq(i,j,k) = dsgn(i,j)*min(dlim(i,j),abs(dq1))
          end do
       end do

    end do

  end subroutine slopey_doit


  subroutine slopez(lo, hi, &
                    q, qlo, qhi, &
                    dq, dqlo, dqhi)

    use amrex_mempool_module, only : bl_allocate, bl_deallocate

    integer, intent(in) :: lo(3), hi(3), qlo(3), qhi(3), dqlo(3), dqhi(3)
    double precision, intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision, intent(out) :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2),dqlo(3):dqhi(3))

    ! Some compiler may not support 'contiguous'.  Remove it in that case.
    double precision, dimension(:,:,:), pointer, contiguous :: dsgn, dlim, df, dcen

    call bl_allocate(dsgn, lo(1), hi(1), lo(2), hi(2), lo(3)-1, hi(3)+1)
    call bl_allocate(dlim, lo(1), hi(1), lo(2), hi(2), lo(3)-1, hi(3)+1)
    call bl_allocate(df  , lo(1), hi(1), lo(2), hi(2), lo(3)-1, hi(3)+1)
    call bl_allocate(dcen, lo(1), hi(1), lo(2), hi(2), lo(3)-1, hi(3)+1)

    call slopez_doit(lo, hi, &
                     q, qlo, qhi, &
                     dq, dqlo, dqhi, &
                     dsgn, dlim, df, dcen, &
                     (/lo(1),lo(2),lo(3)-1/), (/hi(1),hi(2),hi(3)+1/))

    call bl_deallocate(dsgn)
    call bl_deallocate(dlim)
    call bl_deallocate(df)
    call bl_deallocate(dcen)

  end subroutine slopez

  subroutine slopez_doit(lo, hi, &
                         q, qlo, qhi, &
                         dq, dqlo, dqhi, &
                         dsgn, dlim, df, dcen, ddlo, ddhi)

    integer, intent(in) :: lo(3), hi(3), qlo(3), qhi(3), dqlo(3), dqhi(3), &
         ddlo(3), ddhi(3)
    double precision, intent(in ) ::  q  ( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision, intent(out) :: dq  (dqlo(1):dqhi(1),dqlo(2):dqhi(2),dqlo(3):dqhi(3))
    double precision              :: dsgn(ddlo(1):ddhi(1),ddlo(2):ddhi(2),ddlo(3):ddhi(3))
    double precision              :: dlim(ddlo(1):ddhi(1),ddlo(2):ddhi(2),ddlo(3):ddhi(3))
    double precision              :: df  (ddlo(1):ddhi(1),ddlo(2):ddhi(2),ddlo(3):ddhi(3))
    double precision              :: dcen(ddlo(1):ddhi(1),ddlo(2):ddhi(2),ddlo(3):ddhi(3))

    integer :: i, j, k
    double precision :: dlft, drgt, dq1

    ! first compute Fromm slopes
    do k       = lo(3)-1, hi(3)+1
       do j    = lo(2)  , hi(2)
          do i = lo(1)  , hi(1)
             dlft = q(i,j,k  ) - q(i,j,k-1)
             drgt = q(i,j,k+1) - q(i,j,k  )
             dcen(i,j,k) = .5d0 * (dlft+drgt)
             dsgn(i,j,k) = sign( 1.d0, dcen(i,j,k) )
             if (dlft*drgt .ge. 0.d0) then
                dlim(i,j,k) = 2.d0 * min( abs(dlft), abs(drgt) )
             else
                dlim(i,j,k) = 0.d0
             endif
             df(i,j,k) = dsgn(i,j,k)*min( dlim(i,j,k),abs(dcen(i,j,k)) )
          end do
       end do
    end do
       
    ! Now compute limited fourth order slopes
    do k       = lo(3), hi(3)
       do j    = lo(2), hi(2)
          do i = lo(1), hi(1)
             dq1 = four3rd*dcen(i,j,k) - sixth*( df(i,j,k+1) + df(i,j,k-1) )
             dq(i,j,k) = dsgn(i,j,k)*min(dlim(i,j,k),abs(dq1))
          end do
       end do
    end do

  end subroutine slopez_doit

end module slope_module 
