module fornberg_weights_module
  implicit none
  integer, parameter, private :: dp_t = kind(1.0d0)
contains
  subroutine weights(c, z, x)
    real(kind=dp_t), intent(in) :: z
    real(kind=dp_t), dimension(0:), intent(in) :: x
    real(kind=dp_t), dimension(0:,0:), intent(out) :: c
    real(kind=dp_t) :: c1, c2, c3, c4, c5
    integer :: i, j, k, mn
    c1 = 1.0_dp_t
    c4 = x(0) - z
    c(:,:) = 0.0_dp_t
    c(0,0) = 1.0_dp_t
    do i = 1, ubound(x,dim=1)
       mn = min(i, ubound(c,dim=2))
       c2 = 1.0_dp_t
       c5 = c4
       c4 = x(i) - z
       do j = 0, i - 1
          c3 = x(i) - x(j)
          c2 = c2*c3
          if ( j == i-1 ) then
             do k = mn, 1, -1
                c(i,k) = c1*(k*c(i-1,k-1) - c5*c(i-1,k))/c2
             end do
             c(i,0) = -c1*c5*c(i-1,0)/c2
          end if
          do k = mn, 1, -1
             c(j,k) = (c4*c(j,k) - k*c(j,k-1))/c3
          end do
          c(j,0) = c4*c(j,0)/c3
       end do
       c1 = c2
    end do
  end subroutine weights
end module fornberg_weights_module

program main
  use fornberg_weights_module
  implicit none
  integer :: i
  integer, parameter :: dp_t = kind(1.0d0)
  real(kind(0.0d0)) :: z, x(0:2), c(0:2,0:2)
  z = 0
  x = (/-0.5_dp_t, 0.0_dp_t, 1.0_dp_t/)
  call weights(c, z, x)
  do i = lbound(c,dim=2), ubound(c,dim=2)
     print *, '*', i, c(:,i)
  end do
end program main
