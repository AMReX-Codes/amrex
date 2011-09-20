  module X

    implicit none

    contains

      subroutine nsold(omega, ss, uu, ff, lo, ng)

        integer, intent(in)             :: ng, lo(:)
        double precision, intent(in)    :: omega
        double precision, intent(in   ) :: ff(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
        double precision, intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
        double precision, intent(in)    :: ss(lo(1):, lo(2):, lo(3):, 0:)
        integer                         :: i, j, k, ioff, hi(size(lo))
        double precision                :: dd
        logical                         :: jface, kface, doit

        hi = ubound(ff)

        if ((size(ss,dim=4) .eq. 21) .or. (size(ss,dim=4) .eq. 27)) then

           do k = lo(3),hi(3)
              kface = .false. ; if ( (k.eq.lo(3)) .or. (k.eq.hi(3)) ) kface = .true.

              do j = lo(2),hi(2)
                 jface = .false. ; if ( (j.eq.lo(2)) .or. (j.eq.hi(2)) ) jface = .true.

                 do i = lo(1),hi(1)

                    doit = .true.

                    if ( jface .or. kface .or. (i.eq.lo(1)) .or. (i.eq.hi(1)) ) then
                       doit = .false.
                    end if

                    if (doit) then
                       dd = ss(i,j,k,0)*uu(i,j,k) &
                            + ss(i,j,k, 1) * uu(i-1,j-1,k-1) + ss(i,j,k, 2) * uu(i  ,j-1,k-1) &
                            + ss(i,j,k, 3) * uu(i+1,j-1,k-1) + ss(i,j,k, 4) * uu(i-1,j  ,k-1) &
                            + ss(i,j,k, 5) * uu(i+1,j  ,k-1) + ss(i,j,k, 6) * uu(i-1,j+1,k-1) &
                            + ss(i,j,k, 7) * uu(i  ,j+1,k-1) + ss(i,j,k, 8) * uu(i+1,j+1,k-1) &
                            + ss(i,j,k, 9) * uu(i-1,j-1,k  ) + ss(i,j,k,10) * uu(i+1,j-1,k  ) &
                            + ss(i,j,k,11) * uu(i-1,j+1,k  ) + ss(i,j,k,12) * uu(i+1,j+1,k  ) &
                            + ss(i,j,k,13) * uu(i-1,j-1,k+1) + ss(i,j,k,14) * uu(i  ,j-1,k+1) &
                            + ss(i,j,k,15) * uu(i+1,j-1,k+1) + ss(i,j,k,16) * uu(i-1,j  ,k+1) &
                            + ss(i,j,k,17) * uu(i+1,j  ,k+1) + ss(i,j,k,18) * uu(i-1,j+1,k+1) &
                            + ss(i,j,k,19) * uu(i  ,j+1,k+1) + ss(i,j,k,20) * uu(i+1,j+1,k+1)
             
                       uu(i,j,k) = uu(i,j,k) + omega/ss(i,j,k,0)*(ff(i,j,k) - dd)
                    end if
                 end do
              end do
           end do
        end if

      end subroutine nsold

      subroutine nsnew(omega, ss, uu, ff, lo, ng)

        integer, intent(in)             :: ng, lo(:)
        double precision, intent(in)    :: omega
        double precision, intent(in   ) :: ff(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
        double precision, intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
        double precision, intent(in)    :: ss(0:, lo(1):, lo(2):, lo(3):)
        integer                         :: i, j, k, ioff, hi(size(lo))
        double precision                :: dd
        logical                         :: jface, kface, doit

        hi = ubound(ff)

        if ((size(ss,dim=1) .eq. 21) .or. (size(ss,dim=1) .eq. 27)) then

           do k = lo(3),hi(3)
              kface = .false. ; if ( (k.eq.lo(3)) .or. (k.eq.hi(3)) ) kface = .true.

              do j = lo(2),hi(2)
                 jface = .false. ; if ( (j.eq.lo(2)) .or. (j.eq.hi(2)) ) jface = .true.

                 do i = lo(1),hi(1)

                    doit = .true.

                    if ( jface .or. kface .or. (i.eq.lo(1)) .or. (i.eq.hi(1)) ) then
                       doit = .false.
                    end if

                    if (doit) then
                       dd = ss(0,i,j,k)*uu(i,j,k) &
                            + ss(1,i,j,k ) * uu(i-1,j-1,k-1) + ss(2,i,j,k ) * uu(i  ,j-1,k-1) &
                            + ss(3,i,j,k ) * uu(i+1,j-1,k-1) + ss(4,i,j,k ) * uu(i-1,j  ,k-1) &
                            + ss(5,i,j,k ) * uu(i+1,j  ,k-1) + ss(6,i,j,k ) * uu(i-1,j+1,k-1) &
                            + ss(7,i,j,k ) * uu(i  ,j+1,k-1) + ss(8,i,j,k ) * uu(i+1,j+1,k-1) &
                            + ss(9,i,j,k ) * uu(i-1,j-1,k  ) + ss(10,i,j,k) * uu(i+1,j-1,k  ) &
                            + ss(11,i,j,k) * uu(i-1,j+1,k  ) + ss(12,i,j,k) * uu(i+1,j+1,k  ) &
                            + ss(13,i,j,k) * uu(i-1,j-1,k+1) + ss(14,i,j,k) * uu(i  ,j-1,k+1) &
                            + ss(15,i,j,k) * uu(i+1,j-1,k+1) + ss(16,i,j,k) * uu(i-1,j  ,k+1) &
                            + ss(17,i,j,k) * uu(i+1,j  ,k+1) + ss(18,i,j,k) * uu(i-1,j+1,k+1) &
                            + ss(19,i,j,k) * uu(i  ,j+1,k+1) + ss(20,i,j,k) * uu(i+1,j+1,k+1)
             
                       uu(i,j,k) = uu(i,j,k) + omega/ss(0,i,j,k)*(ff(i,j,k) - dd)
                    end if
                 end do
              end do
           end do
        end if

      end subroutine nsnew

  end module X

  program main

    use X

    implicit none

    integer, parameter :: N = 64, NG = 1, lo(3) = (/1,1,1/)

    double precision, parameter :: omega = 1.3d0

    double precision :: ssold(N,N,N,0:20)
    double precision :: ssnew(0:20,N,N,N)
    double precision :: ff(1-NG:N+NG,1-NG:N,1-NG:N+NG)
    double precision :: uu(1-NG:N+NG,1-NG:N+NG,1-NG:N+NG)

!    logical, parameter :: oldway = .false.
    logical, parameter :: oldway = .true.

    integer i

    do i = 0,20
       ssold(:,:,:,i) = 1.0d0 + i
       ssnew(i,:,:,:) = 1.0d0 + i
    end do

    uu = 0.0d0

    ff = 1.0d0

    if ( oldway ) then
       do i = 1,250
          call nsold(omega,ssold,uu,ff,lo,NG)
       end do
    else
       do i = 1,250
          call nsnew(omega,ssnew,uu,ff,lo,NG)
       end do
    end if

  end program main
