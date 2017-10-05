module amrex_multifabutil_nd_module

  use amrex_fort_module, only : amrex_real
  implicit none
  private

  public :: amrex_fort_int_to_real

contains

  subroutine amrex_fort_int_to_real (lo, hi, ncomp, rdata, rlo, rhi, idata, ilo, ihi) &
       bind(c, name='amrex_fort_int_to_real')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, ilo, ihi
    integer, intent(in) :: ncomp
    real(amrex_real), intent(inout) :: rdata(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3),ncomp)
    integer, intent(in) :: idata(ilo(1):ihi(1),ilo(2):ihi(2),ilo(3):ihi(3),ncomp)

    integer :: i,j,k,n

    do n = 1, ncomp
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rdata(i,j,k,n) = idata(i,j,k,n)
             end do
          end do
       end do
    end do
  end subroutine amrex_fort_int_to_real

end module amrex_multifabutil_nd_module
