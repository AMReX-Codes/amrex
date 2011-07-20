module tridiag_module
 
  use bl_error_module
 
  implicit none

contains

  subroutine tridiag(a,b,c,r,u,n)

      integer           , intent(in   ) ::  n
      real (kind = dp_t), intent(in   ) :: a(n), b(n), c(n), r(n)
      real (kind = dp_t), intent(inout) :: u(n)

      integer j
      real (kind = dp_t), allocatable :: gam(:)
      real (kind = dp_t) :: bet

      allocate(gam(n))

      if (b(1) .eq. 0) call bl_error('tridiag: CANT HAVE B(1) = ZERO')

      bet = b(1)
      u(1) = r(1)/bet

      do j = 2,n
        gam(j) = c(j-1)/bet
        bet = b(j) - a(j)*gam(j)
        if (bet .eq. 0) call bl_error('tridiag: TRIDIAG FAILED')
        u(j) = (r(j)-a(j)*u(j-1))/bet
      end do

      do j = n-1,1,-1
        u(j) = u(j) - gam(j+1)*u(j+1)
      end do

  end subroutine tridiag

end module tridiag_module
