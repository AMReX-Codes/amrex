module ode_params
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_long), parameter :: neq = 1
end module ode_params

! Right hand side of ODE for CVODE to solve.

module rhs_mod
  implicit none

  contains

    integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
           result(ierr) bind(C,name='RhsFn')

      use, intrinsic :: iso_c_binding
      use cvode_interface
      use fsundials_nvector_mod
      use fnvector_serial_mod
      use ode_params

      implicit none

      real(c_double), value :: tn
      type(N_Vector)        :: sunvec_f, sunvec_y
      type(c_ptr), value    :: user_data

      ! pointers to data in SUNDAILS vectors
      real(c_double), pointer :: fvec(:)

      ! get data array from SUNDIALS vectors
      fvec => FN_VGetArrayPointer(sunvec_f)

      fvec(1) = 2.0*tn

    end function RhsFn

end module rhs_mod
