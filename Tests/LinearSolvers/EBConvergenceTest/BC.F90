module bc
private
public :: apply_bc
contains
  subroutine apply_bc(lo, hi, phi, phlo, phhi, &
                      dx, problo, probhi)&
                      bind(c, name='apply_bc')
  use amrex_fort_module, only: amrex_real
  implicit none
  integer, dimension(2), intent(in) :: lo, hi, phlo, phhi
  real(amrex_real), dimension(2), intent(in) :: dx, problo, probhi
  real(amrex_real), intent(inout) :: phi(phlo(1):phhi(1), phlo(2):phhi(2))
  
  real(amrex_real) :: xl, xh, x, yl, yh, y, denom
  integer          :: i, j

  xl = problo(1)
  xh = probhi(1)
  yl = problo(2)
  yh = probhi(2)
  do j = lo(2)+1, hi(2)-1
    y = yl + (dble(j) + 0.5d0)*dx(2)
    denom = sqrt((xl-0.5d0)*(xl-0.5d0) + (y-0.5d0)*(y-0.5d0))
    phi(lo(1),j) = (xl-0.5d0)/denom
    denom = sqrt((xh-0.5d0)*(xh-0.5d0) + (y-0.5d0)*(y-0.5d0))
    phi(hi(1),j) = (xh-0.5d0)/denom
  end do
 
  do i = lo(1)+1, hi(1)-1
    x = xl + (dble(i) + 0.5d0)*dx(1)
    denom = sqrt((x-0.5d0)*(x-0.5d0) + (yl-0.5d0)*(yl-0.5d0))
    phi(i,lo(2)) = (x-0.5d0)/denom
    denom = sqrt((x-0.5d0)*(x-0.5d0) + (yh-0.5d0)*(yh-0.5d0))
    phi(i,hi(2)) = (x-0.5d0)/denom
  end do

  denom = sqrt((xl-0.5d0)*(xl-0.5d0) + (yl-0.5d0)*(yl-0.5d0))
  phi(lo(1),lo(2)) = (xl-0.5d0)/denom

  denom = sqrt((xh-0.5d0)*(xh-0.5d0) + (yh-0.5d0)*(yh-0.5d0))
  phi(hi(1), hi(2)) = (xh-0.5d0)/denom

  denom = sqrt((xl-0.5d0)*(xl-0.5d0) + (yh-0.5d0)*(yh-0.5d0))
  phi(lo(1), hi(2)) = (xl-0.5d0)/denom

  denom = sqrt((xh-0.5d0)*(xh-0.5d0) + (yl-0.5d0)*(yl-0.5d0))
  phi(hi(1),lo(2)) = (xh-0.5d0)/denom
  end subroutine apply_bc
end module bc
