
module amrex_extrapolater

  use amrex_fort_module, only : amrex_real

  implicit none
  integer, parameter :: finecell = 1 ! must be consistent with Extrapolater.H
  integer, parameter :: crsecell = 0

  ! The value of msk is either 0 or 1.

contains

  subroutine amrex_first_order_extrap (u, ulo, uhi, nu, msk, mlo, mhi, lo, hi, sc, nc) &
       bind(c,name='amrex_first_order_extrap')

    integer, intent(in) :: ulo(3), uhi(3), nu, mlo(3), mhi(3), lo(3), hi(3), sc, nc
    real(amrex_real), intent(inout) ::   u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),0:nu-1)
    integer     , intent(in)    :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i, j, k, n

    do n = sc, sc+nc-1
       ! set all crse cells to zero first
       do       k = lo(3)-1, hi(3)+1
          do    j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1
                if (msk(i,j,k) .eq. crsecell) then
                   u(i,j,k,n) = 0.d0
                end if
             end do
          end do
       end do

       ! z-lo, y-lo, x-lo
       i = lo(1)-1
       j = lo(2)-1
       k = lo(3)-1
       if (msk(i,j,k) .eq. crsecell) then
          if (      msk(i+1,j,k) .eq. finecell &
               .or. msk(i,j+1,k) .eq. finecell &
               .or. msk(i,j,k+1) .eq. finecell) then

             u(i,j,k,n) = (msk(i+1,j,k) * u(i+1,j,k,n) &
                  +        msk(i,j+1,k) * u(i,j+1,k,n) &
                  +        msk(i,j,k+1) * u(i,j,k+1,n)) &
                  / (msk(i+1,j,k) + msk(i,j+1,k) + msk(i,j,k+1))

          else if ( msk(i+1,j+1,k) .eq. finecell &
               .or. msk(i+1,j,k+1) .eq. finecell &
               .or. msk(i,j+1,k+1) .eq. finecell) then

             u(i,j,k,n) = (msk(i+1,j+1,k) * u(i+1,j+1,k,n) &
                  +        msk(i+1,j,k+1) * u(i+1,j,k+1,n) &
                  +        msk(i,j+1,k+1) * u(i,j+1,k+1,n)) &
                  / (msk(i+1,j+1,k) + msk(i+1,j,k+1) + msk(i,j+1,k+1))
          else
             u(i,j,k,n) = u(i+1,j+1,k+1,n)
          end if
       end if

       ! z-lo, y-lo, x-valid
       j = lo(2)-1
       k = lo(3)-1
       do i = lo(1), hi(1)
          if (msk(i,j,k) .eq. crsecell) then
             if (      msk(i-1,j,k) .eq. finecell &
                  .or. msk(i+1,j,k) .eq. finecell &
                  .or. msk(i,j+1,k) .eq. finecell &
                  .or. msk(i,j,k+1) .eq. finecell) then

                u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                     +        msk(i+1,j,k) * u(i+1,j,k,n) &
                     +        msk(i,j+1,k) * u(i,j+1,k,n) &
                     +        msk(i,j,k+1) * u(i,j,k+1,n)) &
                     / (msk(i-1,j,k) + msk(i+1,j,k) + msk(i,j+1,k) + msk(i,j,k+1))
             else
                u(i,j,k,n) = (msk(i-1,j+1,k) * u(i-1,j+1,k,n) &
                     +        msk(i+1,j+1,k) * u(i+1,j+1,k,n) &
                     +        msk(i-1,j,k+1) * u(i-1,j,k+1,n) &
                     +        msk(i+1,j,k+1) * u(i+1,j,k+1,n) &
                     +                         u(i,j+1,k+1,n)) / &
                     (msk(i-1,j+1,k) + msk(i+1,j+1,k) + msk(i-1,j,k+1) + msk(i+1,j,k+1) + 1)
             end if
          end if
       end do

       ! z-lo, y-lo, x-hi
       i = hi(1)+1
       j = lo(2)-1
       k = lo(3)-1
       if (msk(i,j,k) .eq. crsecell) then
          if (      msk(i-1,j,k) .eq. finecell &
               .or. msk(i,j+1,k) .eq. finecell &
               .or. msk(i,j,k+1) .eq. finecell) then

             u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                  +        msk(i,j+1,k) * u(i,j+1,k,n) &
                  +        msk(i,j,k+1) * u(i,j,k+1,n)) &
                  / (msk(i-1,j,k) + msk(i,j+1,k) + msk(i,j,k+1))

          else if ( msk(i-1,j+1,k) .eq. finecell &
               .or. msk(i-1,j,k+1) .eq. finecell &
               .or. msk(i,j+1,k+1) .eq. finecell) then

             u(i,j,k,n) = (msk(i-1,j+1,k) * u(i-1,j+1,k,n) &
                  +        msk(i-1,j,k+1) * u(i-1,j,k+1,n) &
                  +        msk(i,j+1,k+1) * u(i,j+1,k+1,n)) &
                  / (msk(i-1,j+1,k) + msk(i-1,j,k+1) + msk(i,j+1,k+1))
          else
             u(i,j,k,n) = u(i-1,j+1,k+1,n)
          end if
       end if

       ! z-lo, y-valid, x-lo
       i = lo(1)-1
       k = lo(3)-1
       do j = lo(2), hi(2)
          if (msk(i,j,k) .eq. crsecell) then
             if (      msk(i+1,j,k) .eq. finecell &
                  .or. msk(i,j-1,k) .eq. finecell &
                  .or. msk(i,j+1,k) .eq. finecell &
                  .or. msk(i,j,k+1) .eq. finecell) then

                u(i,j,k,n) = (msk(i+1,j,k) * u(i+1,j,k,n) &
                     +        msk(i,j-1,k) * u(i,j-1,k,n) &
                     +        msk(i,j+1,k) * u(i,j+1,k,n) &
                     +        msk(i,j,k+1) * u(i,j,k+1,n)) / &
                     (msk(i+1,j,k) + msk(i,j-1,k) + msk(i,j+1,k) + msk(i,j,k+1))
             else
                u(i,j,k,n) = (msk(i+1,j-1,k) * u(i+1,j-1,k,n) &
                     +        msk(i+1,j+1,k) * u(i+1,j+1,k,n) &
                     +                         u(i+1,j,k+1,n) &
                     +        msk(i,j-1,k+1) * u(i,j-1,k+1,n) &
                     +        msk(i,j+1,k+1) * u(i,j+1,k+1,n)) / &
                     (msk(i+1,j-1,k) + msk(i+1,j+1,k) + 1 + msk(i,j-1,k+1) + msk(i,j+1,k+1))
             end if
          end if
       end do

       ! z-lo, y-valid, x-valid
       k = lo(3)-1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .eq. crsecell) then
                u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                     +        msk(i+1,j,k) * u(i+1,j,k,n) &
                     +        msk(i,j-1,k) * u(i,j-1,k,n) &
                     +        msk(i,j+1,k) * u(i,j+1,k,n) &
                     +                       u(i,j,k+1,n)) / &
                (msk(i-1,j,k) + msk(i+1,j,k) + msk(i,j-1,k) + msk(i,j+1,k) + 1)
             end if
          end do
       end do

       ! z-lo, y-valid, x-hi
       i = hi(1)+1
       k = lo(3)-1
       do j = lo(2), hi(2)
          if (msk(i,j,k) .eq. crsecell) then
             if (      msk(i-1,j,k) .eq. finecell &
                  .or. msk(i,j-1,k) .eq. finecell &
                  .or. msk(i,j+1,k) .eq. finecell &
                  .or. msk(i,j,k+1) .eq. finecell) then

                u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                     +        msk(i,j-1,k) * u(i,j-1,k,n) &
                     +        msk(i,j+1,k) * u(i,j+1,k,n) &
                     +        msk(i,j,k+1) * u(i,j,k+1,n)) &
                     / (msk(i-1,j,k) + msk(i,j-1,k) + msk(i,j+1,k) + msk(i,j,k+1))
             else
                u(i,j,k,n) = (msk(i-1,j-1,k)*u(i-1,j-1,k,n) &
                     +        msk(i-1,j+1,k)*u(i-1,j+1,k,n) &
                     +                       u(i-1,j,k+1,n) &
                     +        msk(i,j-1,k+1)*u(i,j-1,k+1,n) &
                     +        msk(i,j+1,k+1)*u(i,j+1,k+1,n)) / &
                     (msk(i-1,j-1,k) + msk(i-1,j+1,k) + 1 + msk(i,j-1,k+1) + msk(i,j+1,k+1))
             end if
          end if
       end do

       ! z-lo, y-hi, x-lo
       i = lo(1)-1
       j = hi(2)+1
       k = lo(3)-1
       if (msk(i,j,k) .eq. crsecell) then
          if (      msk(i+1,j,k) .eq. finecell &
               .or. msk(i,j-1,k) .eq. finecell &
               .or. msk(i,j,k+1) .eq. finecell) then

             u(i,j,k,n) = (msk(i+1,j,k) * u(i+1,j,k,n) &
                  +        msk(i,j-1,k) * u(i,j-1,k,n) &
                  +        msk(i,j,k+1) * u(i,j,k+1,n)) &
                  / (msk(i+1,j,k) + msk(i,j-1,k) + msk(i,j,k+1))

          else if ( msk(i+1,j-1,k) .eq. finecell &
               .or. msk(i+1,j,k+1) .eq. finecell &
               .or. msk(i,j-1,k+1) .eq. finecell) then

             u(i,j,k,n) = (msk(i+1,j-1,k) * u(i+1,j-1,k,n) &
                  +        msk(i+1,j,k+1) * u(i+1,j,k+1,n) &
                  +        msk(i,j-1,k+1) * u(i,j-1,k+1,n)) &
                  / (msk(i+1,j-1,k) + msk(i+1,j,k+1) + msk(i,j-1,k+1))
          else
             u(i,j,k,n) = u(i+1,j-1,k+1,n)
          end if
       end if

       ! z-lo, y-hi, x-valid
       j = hi(2)+1
       k = lo(3)-1
       do i = lo(1), hi(1)
          if (msk(i,j,k) .eq. crsecell) then
             if (      msk(i-1,j,k) .eq. finecell &
                  .or. msk(i+1,j,k) .eq. finecell &
                  .or. msk(i,j-1,k) .eq. finecell &
                  .or. msk(i,j,k+1) .eq. finecell) then

                u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                     +        msk(i+1,j,k) * u(i+1,j,k,n) &
                     +        msk(i,j-1,k) * u(i,j-1,k,n) &
                     +        msk(i,j,k+1) * u(i,j,k+1,n)) &
                     / (msk(i-1,j,k) + msk(i+1,j,k) + msk(i,j-1,k) + msk(i,j,k+1))
             else
                u(i,j,k,n) = (msk(i-1,j-1,k) * u(i-1,j-1,k,n) &
                     +        msk(i+1,j-1,k) * u(i+1,j-1,k,n) &
                     +        msk(i-1,j,k+1) * u(i-1,j,k+1,n) &
                     +        msk(i+1,j,k+1) * u(i+1,j,k+1,n) &
                     +                         u(i,j-1,k+1,n)) / &
                     (msk(i-1,j-1,k) + msk(i+1,j-1,k) + msk(i-1,j,k+1) + msk(i+1,j,k+1) + 1)
             end if
          end if
       end do

       ! z-lo, y-hi, x-hi
       i = hi(1)+1
       j = hi(2)+1
       k = lo(3)-1
       if (msk(i,j,k) .eq. crsecell) then
          if (      msk(i-1,j,k) .eq. finecell &
               .or. msk(i,j-1,k) .eq. finecell &
               .or. msk(i,j,k+1) .eq. finecell) then

             u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                  +        msk(i,j-1,k) * u(i,j-1,k,n) &
                  +        msk(i,j,k+1) * u(i,j,k+1,n)) &
                  / (msk(i-1,j,k) + msk(i,j-1,k) + msk(i,j,k+1))

          else if ( msk(i-1,j-1,k) .eq. finecell &
               .or. msk(i-1,j,k+1) .eq. finecell &
               .or. msk(i,j-1,k+1) .eq. finecell) then

             u(i,j,k,n) = (msk(i-1,j-1,k) * u(i-1,j-1,k,n) &
                  +        msk(i-1,j,k+1) * u(i-1,j,k+1,n) &
                  +        msk(i,j-1,k+1) * u(i,j-1,k+1,n)) &
                  / (msk(i-1,j-1,k) + msk(i-1,j,k+1) + msk(i,j-1,k+1))
          else
             u(i,j,k,n) = u(i-1,j-1,k+1,n)
          end if
       end if

       ! z-valid, y-lo, x-lo
       i = lo(1)-1
       j = lo(2)-1
       do k = lo(3), hi(3)
          if (msk(i,j,k) .eq. crsecell) then
             if (      msk(i+1,j,k) .eq. finecell &
                  .or. msk(i,j+1,k) .eq. finecell &
                  .or. msk(i,j,k-1) .eq. finecell &
                  .or. msk(i,j,k+1) .eq. finecell) then

                u(i,j,k,n) = (msk(i+1,j,k) * u(i+1,j,k,n) &
                     +        msk(i,j+1,k) * u(i,j+1,k,n) &
                     +        msk(i,j,k-1) * u(i,j,k-1,n) &
                     +        msk(i,j,k+1) * u(i,j,k+1,n)) &
                     / (msk(i+1,j,k) + msk(i,j+1,k) + msk(i,j,k-1) + msk(i,j,k+1))
             else
                u(i,j,k,n) = (                 u(i+1,j+1,k,n) &
                     +        msk(i+1,j,k-1) * u(i+1,j,k-1,n) &
                     +        msk(i+1,j,k+1) * u(i+1,j,k+1,n) &
                     +        msk(i,j+1,k-1) * u(i,j+1,k-1,n) &
                     +        msk(i,j+1,k+1) * u(i,j+1,k+1,n)) / &
                     (1 + msk(i+1,j,k-1) + msk(i+1,j,k+1) + msk(i,j+1,k-1) + msk(i,j+1,k+1))
             end if
          end if
       end do

       ! z-valid, y-lo, x-valid
       j = lo(2)-1
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .eq. crsecell) then
                u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                     +        msk(i+1,j,k) * u(i+1,j,k,n) &
                     +                       u(i,j+1,k,n) &
                     +        msk(i,j,k-1) * u(i,j,k-1,n) &
                     +        msk(i,j,k+1) * u(i,j,k+1,n)) / &
                     (msk(i-1,j,k) + msk(i+1,j,k) + 1 + msk(i,j,k-1) + msk(i,j,k+1))
             end if
          end do
       end do

       ! z-valid, y-lo, x-hi
       i = hi(1)+1
       j = lo(2)-1
       do k = lo(3), hi(3)
          if (msk(i,j,k) .eq. crsecell) then
             if (      msk(i-1,j,k) .eq. finecell &
                  .or. msk(i,j+1,k) .eq. finecell &
                  .or. msk(i,j,k-1) .eq. finecell &
                  .or. msk(i,j,k+1) .eq. finecell) then

                u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                     +        msk(i,j+1,k) * u(i,j+1,k,n) &
                     +        msk(i,j,k-1) * u(i,j,k-1,n) &
                     +        msk(i,j,k+1) * u(i,j,k+1,n)) &
                     / (msk(i-1,j,k) + msk(i,j+1,k) + msk(i,j,k-1) + msk(i,j,k+1))
             else
                u(i,j,k,n) = (                 u(i-1,j+1,k,n) &
                     +        msk(i-1,j,k-1) * u(i-1,j,k-1,n) &
                     +        msk(i-1,j,k+1) * u(i-1,j,k+1,n) &
                     +        msk(i,j+1,k-1) * u(i,j+1,k-1,n) &
                     +        msk(i,j+1,k+1) * u(i,j+1,k+1,n)) / &
                     (1 + msk(i-1,j,k-1) + msk(i-1,j,k+1) + msk(i,j+1,k-1) + msk(i,j+1,k+1))
             end if
          end if
       end do

       ! z-valid, y-valid, x-lo
       i = lo(1)-1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             if (msk(i,j,k) .eq. crsecell) then
                u(i,j,k,n) = (               u(i+1,j,k,n) &
                     +        msk(i,j-1,k) * u(i,j-1,k,n) &
                     +        msk(i,j+1,k) * u(i,j+1,k,n) &
                     +        msk(i,j,k-1) * u(i,j,k-1,n) &
                     +        msk(i,j,k+1) * u(i,j,k+1,n)) / &
                     (1 + msk(i,j-1,k) + msk(i,j+1,k) + msk(i,j,k-1) + msk(i,j,k+1))
             end if
          end do
       end do
       
       ! z-valid, y-valid, x-hi
       i = hi(1)+1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             if (msk(i,j,k) .eq. crsecell) then
                u(i,j,k,n) = (               u(i-1,j,k,n) &
                     +        msk(i,j-1,k) * u(i,j-1,k,n) &
                     +        msk(i,j+1,k) * u(i,j+1,k,n) &
                     +        msk(i,j,k-1) * u(i,j,k-1,n) &
                     +        msk(i,j,k+1) * u(i,j,k+1,n)) / &
                (1 + msk(i,j-1,k) + msk(i,j+1,k) + msk(i,j,k-1) + msk(i,j,k+1))
             end if
          end do
       end do
       
       ! z-valid, y-hi, x-lo
       i = lo(1)-1
       j = hi(2)+1
       do k = lo(3), hi(3)
          if (msk(i,j,k) .eq. crsecell) then
             if (      msk(i+1,j,k) .eq. finecell &
                  .or. msk(i,j-1,k) .eq. finecell &
                  .or. msk(i,j,k-1) .eq. finecell &
                  .or. msk(i,j,k+1) .eq. finecell) then

                u(i,j,k,n) = (msk(i+1,j,k) * u(i+1,j,k,n) &
                     +        msk(i,j-1,k) * u(i,j-1,k,n) &
                     +        msk(i,j,k-1) * u(i,j,k-1,n) &
                     +        msk(i,j,k+1) * u(i,j,k+1,n)) &
                     / (msk(i+1,j,k) + msk(i,j-1,k) + msk(i,j,k-1) + msk(i,j,k+1))
             else
                u(i,j,k,n) = (                 u(i+1,j-1,k,n) &
                     +        msk(i+1,j,k-1) * u(i+1,j,k-1,n) &
                     +        msk(i+1,j,k+1) * u(i+1,j,k+1,n) &
                     +        msk(i,j-1,k-1) * u(i,j-1,k-1,n) &
                     +        msk(i,j-1,k+1) * u(i,j-1,k+1,n)) &
                     / (1 + msk(i+1,j,k-1) + msk(i+1,j,k+1) + msk(i,j-1,k-1) + msk(i,j-1,k+1))
             end if
          end if
       end do

       ! z-valid, y-hi, x-valid
       j = hi(2)+1
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .eq. crsecell) then
                u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                     +        msk(i+1,j,k) * u(i+1,j,k,n) &
                     +                       u(i,j-1,k,n) &
                     +        msk(i,j,k-1) * u(i,j,k-1,n) &
                     +        msk(i,j,k+1) * u(i,j,k+1,n)) / &
                     (msk(i-1,j,k) + msk(i+1,j,k) + 1 + msk(i,j,k-1) + msk(i,j,k+1))
             end if
          end do
       end do

       ! z-valid, y-hi, x-hi
       i = hi(1)+1
       j = hi(2)+1
       do k = lo(3), hi(3)
          if (msk(i,j,k) .eq. crsecell) then
             if (      msk(i-1,j,k) .eq. finecell &
                  .or. msk(i,j-1,k) .eq. finecell &
                  .or. msk(i,j,k-1) .eq. finecell &
                  .or. msk(i,j,k+1) .eq. finecell) then

                u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                     +        msk(i,j-1,k) * u(i,j-1,k,n) &
                     +        msk(i,j,k-1) * u(i,j,k-1,n) &
                     +        msk(i,j,k+1) * u(i,j,k+1,n)) &
                     / (msk(i-1,j,k) + msk(i,j-1,k) + msk(i,j,k-1) + msk(i,j,k+1))
             else
                u(i,j,k,n) = (                 u(i-1,j-1,k,n) &
                     +        msk(i-1,j,k-1) * u(i-1,j,k-1,n) &
                     +        msk(i-1,j,k+1) * u(i-1,j,k+1,n) &
                     +        msk(i,j-1,k-1) * u(i,j-1,k-1,n) &
                     +        msk(i,j-1,k+1) * u(i,j-1,k+1,n)) / &
                     (1 + msk(i-1,j,k-1) + msk(i-1,j,k+1) + msk(i,j-1,k-1) + msk(i,j-1,k+1))
             end if
          end if
       end do

       ! z-hi, y-lo, x-lo
       i = lo(1)-1
       j = lo(2)-1
       k = hi(3)+1
       if (msk(i,j,k) .eq. crsecell) then
          if (      msk(i+1,j,k) .eq. finecell &
               .or. msk(i,j+1,k) .eq. finecell &
               .or. msk(i,j,k-1) .eq. finecell) then

             u(i,j,k,n) = (msk(i+1,j,k) * u(i+1,j,k,n) &
                  +        msk(i,j+1,k) * u(i,j+1,k,n) &
                  +        msk(i,j,k-1) * u(i,j,k-1,n)) &
                  / (msk(i+1,j,k) + msk(i,j+1,k) + msk(i,j,k-1))

          else if ( msk(i+1,j+1,k) .eq. finecell &
               .or. msk(i+1,j,k-1) .eq. finecell &
               .or. msk(i,j+1,k-1) .eq. finecell) then

             u(i,j,k,n) = (msk(i+1,j+1,k) * u(i+1,j+1,k,n) &
                  +        msk(i+1,j,k-1) * u(i+1,j,k-1,n) &
                  +        msk(i,j+1,k-1) * u(i,j+1,k-1,n)) &
                  / (msk(i+1,j+1,k) + msk(i+1,j,k-1) + msk(i,j+1,k-1))
          else
             u(i,j,k,n) = u(i+1,j+1,k-1,n)
          end if
       end if

       ! z-hi, y-lo, x-valid
       j = lo(2)-1
       k = hi(3)+1
       do i = lo(1), hi(1)
          if (msk(i,j,k) .eq. crsecell) then
             if (      msk(i-1,j,k) .eq. finecell &
                  .or. msk(i+1,j,k) .eq. finecell &
                  .or. msk(i,j+1,k) .eq. finecell &
                  .or. msk(i,j,k-1) .eq. finecell) then

                u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                     +        msk(i+1,j,k) * u(i+1,j,k,n) &
                     +        msk(i,j+1,k) * u(i,j+1,k,n) &
                     +        msk(i,j,k-1) * u(i,j,k-1,n)) &
                     / (msk(i-1,j,k) + msk(i+1,j,k) + msk(i,j+1,k) + msk(i,j,k-1))
             else
                u(i,j,k,n) = (msk(i-1,j+1,k) * u(i-1,j+1,k,n) &
                     +        msk(i+1,j+1,k) * u(i+1,j+1,k,n) &
                     +        msk(i-1,j,k-1) * u(i-1,j,k-1,n) &
                     +        msk(i+1,j,k-1) * u(i+1,j,k-1,n) &
                     +                         u(i,j+1,k-1,n)) / &
                     (msk(i-1,j+1,k) + msk(i+1,j+1,k) + msk(i-1,j,k-1) + msk(i+1,j,k-1) + 1)
             end if                
          end if
       end do

       ! z-hi, y-lo, x-hi
       i = hi(1)+1
       j = lo(2)-1
       k = hi(3)+1
       if (msk(i,j,k) .eq. crsecell) then
          if (      msk(i-1,j,k) .eq. finecell &
               .or. msk(i,j+1,k) .eq. finecell &
               .or. msk(i,j,k-1) .eq. finecell) then

             u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                  +        msk(i,j+1,k) * u(i,j+1,k,n) &
                  +        msk(i,j,k-1) * u(i,j,k-1,n)) &
                  / (msk(i-1,j,k) + msk(i,j+1,k) + msk(i,j,k-1))

          else if ( msk(i-1,j+1,k) .eq. finecell &
               .or. msk(i-1,j,k-1) .eq. finecell &
               .or. msk(i,j+1,k-1) .eq. finecell) then

             u(i,j,k,n) = (msk(i-1,j+1,k) * u(i-1,j+1,k,n) &
                  +        msk(i-1,j,k-1) * u(i-1,j,k-1,n) &
                  +        msk(i,j+1,k-1) * u(i,j+1,k-1,n)) &
                  / (msk(i-1,j+1,k) + msk(i-1,j,k-1) + msk(i,j+1,k-1))
          else
             u(i,j,k,n) = u(i-1,j+1,k-1,n)
          end if
       end if

       ! z-hi, y-valid, x-lo
       i = lo(1)-1
       k = hi(3)+1
       do j = lo(2), hi(2)
          if (msk(i,j,k) .eq. crsecell) then
             if (      msk(i+1,j,k) .eq. finecell &
                  .or. msk(i,j-1,k) .eq. finecell &
                  .or. msk(i,j+1,k) .eq. finecell &
                  .or. msk(i,j,k-1) .eq. finecell) then

                u(i,j,k,n) = (msk(i+1,j,k) * u(i+1,j,k,n) &
                     +        msk(i,j-1,k) * u(i,j-1,k,n) &
                     +        msk(i,j+1,k) * u(i,j+1,k,n) &
                     +        msk(i,j,k-1) * u(i,j,k-1,n)) / &
                     (msk(i+1,j,k) + msk(i,j-1,k) + msk(i,j+1,k) + msk(i,j,k-1))
             else
                u(i,j,k,n) = (msk(i+1,j-1,k) * u(i+1,j-1,k,n) &
                     +        msk(i+1,j+1,k) * u(i+1,j+1,k,n) &
                     +                         u(i+1,j,k-1,n) &
                     +        msk(i,j-1,k-1) * u(i,j-1,k-1,n) &
                     +        msk(i,j+1,k-1) * u(i,j+1,k-1,n)) / &
                     (msk(i+1,j-1,k) + msk(i+1,j+1,k) + 1 + msk(i,j-1,k-1) + msk(i,j+1,k-1))
             end if                
          end if
       end do

       ! z-hi, y-valid, x-valid
       k = hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .eq. crsecell) then
                u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                     +        msk(i+1,j,k) * u(i+1,j,k,n) &
                     +        msk(i,j-1,k) * u(i,j-1,k,n) &
                     +        msk(i,j+1,k) * u(i,j+1,k,n) &
                     +                       u(i,j,k-1,n)) / &
                     (msk(i-1,j,k) + msk(i+1,j,k) + msk(i,j-1,k) + msk(i,j+1,k) + 1)
             end if
          end do
       end do

       ! z-hi, y-valid, x-hi
       i = hi(1)+1
       k = hi(3)+1
       do j = lo(2), hi(2)
          if (msk(i,j,k) .eq. crsecell) then
             if (      msk(i-1,j,k) .eq. finecell &
                  .or. msk(i,j-1,k) .eq. finecell &
                  .or. msk(i,j+1,k) .eq. finecell &
                  .or. msk(i,j,k-1) .eq. finecell) then

                u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                     +        msk(i,j-1,k) * u(i,j-1,k,n) &
                     +        msk(i,j+1,k) * u(i,j+1,k,n) &
                     +        msk(i,j,k-1) * u(i,j,k-1,n)) / &
                     (msk(i-1,j,k) + msk(i,j-1,k) + msk(i,j+1,k) + msk(i,j,k-1))
             else
                u(i,j,k,n) = (msk(i-1,j-1,k) * u(i-1,j-1,k,n) &
                     +        msk(i-1,j+1,k) * u(i-1,j+1,k,n) &
                     +                         u(i-1,j,k-1,n) &
                     +        msk(i,j-1,k-1) * u(i,j-1,k-1,n) &
                     +        msk(i,j+1,k-1) * u(i,j+1,k-1,n)) / &
                     (msk(i-1,j-1,k) + msk(i-1,j+1,k) + 1 + msk(i,j-1,k-1) + msk(i,j+1,k-1))
             end if                                
          end if
       end do

       ! z-hi, y-hi, x-lo
       i = lo(1)-1
       j = hi(2)+1
       k = hi(3)+1
       if (msk(i,j,k) .eq. crsecell) then
          if (      msk(i+1,j,k) .eq. finecell &
               .or. msk(i,j-1,k) .eq. finecell &
               .or. msk(i,j,k-1) .eq. finecell) then

             u(i,j,k,n) = (msk(i+1,j,k) * u(i+1,j,k,n) &
                  +        msk(i,j-1,k) * u(i,j-1,k,n) &
                  +        msk(i,j,k-1) * u(i,j,k-1,n)) &
                  / (msk(i+1,j,k) + msk(i,j-1,k) + msk(i,j,k-1))

          else if ( msk(i+1,j-1,k) .eq. finecell &
               .or. msk(i+1,j,k-1) .eq. finecell &
               .or. msk(i,j-1,k-1) .eq. finecell) then

             u(i,j,k,n) = (msk(i+1,j-1,k) * u(i+1,j-1,k,n) &
                  +        msk(i+1,j,k-1) * u(i+1,j,k-1,n) &
                  +        msk(i,j-1,k-1) * u(i,j-1,k-1,n)) &
                  / (msk(i+1,j-1,k) + msk(i+1,j,k-1) + msk(i,j-1,k-1))
          else
             u(i,j,k,n) = u(i+1,j-1,k-1,n)
          end if
       end if

       ! z-hi, y-hi, x-valid
       j = hi(2)+1
       k = hi(3)+1
       do i = lo(1), hi(1)
          if (msk(i,j,k) .eq. crsecell) then
             if (      msk(i-1,j,k) .eq. finecell &
                  .or. msk(i+1,j,k) .eq. finecell &
                  .or. msk(i,j-1,k) .eq. finecell &
                  .or. msk(i,j,k-1) .eq. finecell) then

                u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                     +        msk(i+1,j,k) * u(i+1,j,k,n) &
                     +        msk(i,j-1,k) * u(i,j-1,k,n) &
                     +        msk(i,j,k-1) * u(i,j,k-1,n)) &
                     / (msk(i-1,j,k) + msk(i+1,j,k) + msk(i,j-1,k) + msk(i,j,k-1))
             else
                u(i,j,k,n) = (msk(i-1,j-1,k) * u(i-1,j-1,k,n) &
                     +        msk(i+1,j-1,k) * u(i+1,j-1,k,n) &
                     +        msk(i-1,j,k-1) * u(i-1,j,k-1,n) &
                     +        msk(i+1,j,k-1) * u(i+1,j,k-1,n) &
                     +                         u(i,j-1,k-1,n)) / &
                     (msk(i-1,j-1,k) + msk(i+1,j-1,k) + msk(i-1,j,k-1) + msk(i+1,j,k-1) + 1)
             end if                
          end if
       end do

       ! z-hi, y-hi, x-hi
       i = hi(1)+1
       j = hi(2)+1
       k = hi(3)+1
       if (msk(i,j,k) .eq. crsecell) then
          if (      msk(i-1,j,k) .eq. finecell &
               .or. msk(i,j-1,k) .eq. finecell &
               .or. msk(i,j,k-1) .eq. finecell) then

             u(i,j,k,n) = (msk(i-1,j,k) * u(i-1,j,k,n) &
                  +        msk(i,j-1,k) * u(i,j-1,k,n) &
                  +        msk(i,j,k-1) * u(i,j,k-1,n)) &
                  / (msk(i-1,j,k) + msk(i,j-1,k) + msk(i,j,k-1))

          else if ( msk(i-1,j-1,k) .eq. finecell &
               .or. msk(i-1,j,k-1) .eq. finecell &
               .or. msk(i,j-1,k-1) .eq. finecell) then

             u(i,j,k,n) = (msk(i-1,j-1,k) * u(i-1,j-1,k,n) &
                  +        msk(i-1,j,k-1) * u(i-1,j,k-1,n) &
                  +        msk(i,j-1,k-1) * u(i,j-1,k-1,n)) &
                  / (msk(i-1,j-1,k) + msk(i-1,j,k-1) + msk(i,j-1,k-1))
          else
             u(i,j,k,n) = u(i-1,j-1,k-1,n)
          end if
       end if
    end do

  end subroutine amrex_first_order_extrap

end module amrex_extrapolater
