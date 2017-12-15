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

  subroutine amrex_fill_slice(tlo, thi, & 
       full_data, flo, fhi, &
       slice_data, slo, shi, & 
       fstart, nfull, ncomp) &
       bind(C, name="amrex_fill_slice")

    use amrex_fort_module, only : amrex_real
    
    integer       ,   intent(in)    :: ncomp, fstart, nfull
    integer       ,   intent(in)    :: flo(3), fhi(3)
    integer       ,   intent(in)    :: slo(3), shi(3)
    integer       ,   intent(in)    :: tlo(3), thi(3)
    real(amrex_real), intent(inout) ::  full_data(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3), nfull)
    real(amrex_real), intent(inout) :: slice_data(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3), ncomp)
    
    integer n, i, j, k
    
    do n = 1, ncomp
       do k = tlo(3), thi(3)
          do j = tlo(2), thi(2)
             do i = tlo(1), thi(1)
                slice_data(i, j, k, n) = full_data(i, j, k, fstart+n)
             end do
          end do
       end do
    end do
    
  end subroutine amrex_fill_slice
  
end module amrex_multifabutil_nd_module
