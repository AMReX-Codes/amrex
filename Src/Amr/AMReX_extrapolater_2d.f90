BL
module bl_extrapolater

  use amrex_fort_module, only : amrex_real

  implicit none
  integer, parameter :: finecell = 1 ! must be consistent with Extrapolater.H
  integer, parameter :: crsecell = 0

  ! The value of msk is either 0 or 1.

contains

  subroutine first_order_extrap (u, ulo, uhi, nu, msk, mlo, mhi, lo, hi, sc, nc) &
       bind(c,name='first_order_extrap')

    integer, intent(in) :: ulo(2), uhi(2), nu, mlo(2), mhi(2), lo(2), hi(2), sc, nc
    real(amrex_real), intent(inout) ::   u(ulo(1):uhi(1),ulo(2):uhi(2),0:nu-1)
    integer     , intent(in)    :: msk(mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: i, j, n

    do n = sc, sc+nc-1
       ! set all crse cells to zero first
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             if (msk(i,j) .eq. crsecell) then
                u(i,j,n) = 0.d0
             end if
          end do
       end do

       ! ylo, xlo
       j = lo(2)-1
       i = lo(1)-1
       if (msk(i,j) .eq. crsecell) then
          if (msk(i,j+1) .eq. finecell .or. msk(i+1,j) .eq. finecell) then
             u(i,j,n) = (msk(i,j+1)*u(i,j+1,n) + msk(i+1,j)*u(i+1,j,n)) &
                  &   / (msk(i,j+1)            + msk(i+1,j)           )
          else
             u(i,j,n) = u(i+1,j+1,n)
          end if
       end if

       ! ylo, x-valid
       j = lo(2)-1
       do i = lo(1), hi(1)
          if (msk(i,j) .eq. crsecell) then
             u(i,j,n) = (msk(i-1,j)*u(i-1,j,n)+msk(i+1,j)*u(i+1,j,n)+u(i,j+1,n)) &
                  &   / (msk(i-1,j)           +msk(i+1,j)           +1         )
          end if
       end do

       ! ylo, xhi
       j = lo(2)-1
       i = hi(1)+1
       if (msk(i,j) .eq. crsecell) then
          if (msk(i-1,j).eq.finecell .or. msk(i,j+1).eq.finecell) then
             u(i,j,n) = (msk(i-1,j)*u(i-1,j,n)+msk(i,j+1)*u(i,j+1,n)) &
                  &   / (msk(i-1,j)           +msk(i,j+1))
          else
             u(i,j,n) = u(i-1,j+1,n)
          end if
       end if

       ! y-valid, xlo
       i = lo(1)-1
       do j = lo(2), hi(2)
          if (msk(i,j) .eq. crsecell) then
             u(i,j,n) = (msk(i,j-1)*u(i,j-1,n)+u(i+1,j,n)+msk(i,j+1)*u(i,j+1,n)) &
                  &   / (msk(i,j-1)           +1         +msk(i,j+1)           )
          end if
       end do

       ! y-valid, xhi
       i = hi(1)+1
       do j = lo(2), hi(2)
          if (msk(i,j) .eq. crsecell) then
             u(i,j,n) = (msk(i,j-1)*u(i,j-1,n)+u(i-1,j,n)+msk(i,j+1)*u(i,j+1,n)) &
                  &   / (msk(i,j-1)           +1         +msk(i,j+1)           )
          end if
       end do
       
       ! yhi, xlo
       j = hi(2)+1
       i = lo(1)-1
       if (msk(i,j) .eq. crsecell) then
          if (msk(i,j-1).eq.finecell .or. msk(i+1,j).eq.finecell) then
             u(i,j,n) = (msk(i,j-1)*u(i,j-1,n)+msk(i+1,j)*u(i+1,j,n)) &
                  &   / (msk(i,j-1)           +msk(i+1,j)           )
          else
             u(i,j,n) = u(i+1,j-1,n)
          end if
       end if

       ! yhi, xvalid
       j = hi(2)+1
       do i = lo(1), hi(1)
          if (msk(i,j) .eq. crsecell) then
             u(i,j,n) = (u(i,j-1,n)+msk(i-1,j)*u(i-1,j,n)+msk(i+1,j)*u(i+1,j,n)) &
                  &   / (1         +msk(i-1,j)           +msk(i+1,j)           )
          end if
       end do

       ! yhi, xhi
       i = hi(1)+1
       j = hi(2)+1
       if (msk(i,j) .eq. crsecell) then
          if (msk(i-1,j).eq.finecell .or. msk(i,j-1).eq.finecell) then
             u(i,j,n) = (msk(i-1,j)*u(i-1,j,n)+msk(i,j-1)*u(i,j-1,n)) &
                  &   / (msk(i-1,j)           +msk(i,j-1)           )
          else
             u(i,j,n) = u(i-1,j-1,n)
          end if
       end if
    end do

  end subroutine first_order_extrap

end module bl_extrapolater
