
module amrex_extrapolater

  use amrex_fort_module, only : amrex_real

  implicit none
  integer, parameter :: finecell = 1 ! must be consistent with Extrapolater.H
  integer, parameter :: crsecell = 0

contains

  subroutine amrex_first_order_extrap (u, ulo, uhi, nu, msk, mlo, mhi, lo, hi, sc, nc) &
       bind(c,name='amrex_first_order_extrap')

    integer, intent(in) :: ulo(1), uhi(1), nu, mlo(1), mhi(1), lo(1), hi(1), sc, nc
    real(amrex_real), intent(inout) :: u(ulo(1):uhi(1),0:nu-1)
    integer     , intent(in)  :: msk(mlo(1):mhi(1))

    integer :: n

    do n = sc, sc+nc-1
       if (msk(lo(1)-1) .eq. crsecell) then
          u(lo(1)-1,n) = u(lo(1),n)
       end if

       if (msk(hi(1)+1) .eq. crsecell) then
          u(hi(1)+1,n) = u(hi(1),n)
       end if
    end do

  end subroutine amrex_first_order_extrap

end module amrex_extrapolater
