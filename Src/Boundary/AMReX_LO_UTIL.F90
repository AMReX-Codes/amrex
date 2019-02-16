
module amrex_lo_util_module

  use amrex_fort_module
  use amrex_constants_module
  implicit none

contains

!>     polyInterpCoeff:
!!
!!     This routine returns the Lagrange interpolating coefficients for a
!!     polynomial through N points, evaluated at xInt (see Numerical Recipes,
!!     v2, p102, e.g.):
!! ``
!!            (x-x2)(x-x3)...(x-xN)              (x-x1)(x-x2)...(x-x(N-1))
!!    P(x) = ----------------------- y1  + ... + ------------------------  yN
!!           (x1-x2)(x1-x3)...(x1-xN)            (x1-x2)(x1-x3)...(x1-xN)
!!
!!     P(xInt) = sum_(i=1)^(N) y[i]*c[i]``
!!
    subroutine polyInterpCoeff(xInt, x, N, c)

      implicit none

      integer N, i, j
      real(amrex_real) xInt, x(N), c(N), num, den

      do j=1,N
         num = one
         den = one
         do i = 1,j-1
            num = num*(xInt - x(i))
            den = den*(x(j) - x(i))
         end do
         do i = j+1,N
            num = num*(xInt - x(i))
            den = den*(x(j) - x(i))
         end do
#ifdef AMREX_DEBUG
         if (den .eq. zero) STOP 'polyInterpCoeff::invalid data'
#endif
         c(j) = num/den
      end do

    end subroutine polyInterpCoeff

end module amrex_lo_util_module
