
module amrex_ebmultifabutil_module

  use amrex_fort_module, only : amrex_real
  use amrex_ebcellflag_module, only : is_covered_cell
  implicit none

  private
  
  public :: amrex_eb_set_covered

contains

  subroutine amrex_eb_set_covered (lo, hi, d, dlo, dhi, f, flo, fhi, v, nc) &
       bind(c,name='amrex_eb_set_covered')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), flo(3), fhi(3), nc
    real(amrex_real), intent(inout) :: d(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nc)
    integer, intent(in) :: f(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(amrex_real), intent(in) :: v(nc)

    integer :: i, j, k, n

    do n = 1, nc
       do       k = lo(3),hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (is_covered_cell(f(i,j,k))) then
                   d(i,j,k,n) = v(n)
                end if
             end do
          end do
       end do
    end do

  end subroutine amrex_eb_set_covered

end module amrex_ebmultifabutil_module
