
module amrex_abec_util_module

  use amrex_fort_module
  implicit none

contains

!-----------------------------------------------------------------------
!>
!>     Tridiagonal solve
!>
    subroutine tridiag(a,b,c,r,u,n)

      integer n
      integer nmax

      real(amrex_real) a(n)
      real(amrex_real) b(n)
      real(amrex_real) c(n)
      real(amrex_real) r(n)
      real(amrex_real) u(n)

      parameter (nmax = 4098)

      integer j
      real(amrex_real) bet
      real(amrex_real) gam(nmax)
      if (n .gt. nmax ) call bl_error('tridiag: size exceeded')
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

end module amrex_abec_util_module
