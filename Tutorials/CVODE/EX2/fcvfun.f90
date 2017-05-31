! Right hand side of ODE for CVODE to solve. Note that CVODE requires this function to have exactly this name.
subroutine fcvfun(t, y, ydot, ipar, rpar, ier)
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double) :: t
  real(c_double), dimension(*) :: y
  real(c_double), dimension(*) :: ydot
  integer(c_long), dimension(*) :: ipar
  real(c_double), dimension(*) :: rpar
  integer(c_int) :: ier

  YDOT(1) = -0.04D0 * Y(1) + 1.0D4 * Y(2) * Y(3)
  YDOT(3) = 3.0D7 * Y(2) * Y(2)
  YDOT(2) = -YDOT(1) - YDOT(3)

  ! CVODE looks for a zero return code, or else it will print warnings or errors.
  ier = 0
end subroutine fcvfun
