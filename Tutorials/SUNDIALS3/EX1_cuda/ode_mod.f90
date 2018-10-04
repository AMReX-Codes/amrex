module ode_params
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_long), parameter :: neq = 1
end module ode_params

! Wrapper for Const to test initial setup
module fnvector_cuda_mod

  !======= Interfaces =========
  interface

     ! -----------------------------------------------------------------
     ! N_VConst_Cuda
     ! -----------------------------------------------------------------

     subroutine FN_VConst_Cuda(c, z) &
          bind(C,name='N_VConst_Cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: c
       type(c_ptr),    value :: z
     end subroutine FN_VConst_Cuda

     end interface
end module

! Right hand side of ODE for CVODE to solve.

module rhs_mod
  implicit none
  contains

    integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
           result(ierr) bind(C,name='RhsFn')

      use, intrinsic :: iso_c_binding
      use fnvector_cuda_mod
      use ode_params

      implicit none

      real(c_double), value :: tn
      type(c_ptr), value    :: sunvec_y
      type(c_ptr), value    :: sunvec_f
      type(c_ptr), value    :: user_data

      ! pointers to data in SUNDAILS vectors
      real(c_double), pointer :: yvec(:)
      real(c_double), pointer :: fvec(:)

      ! get data arrays from SUNDIALS vectors
!      call N_VGetData_Serial(sunvec_f, neq, fvec)
      call FN_VConst_Cuda(2.0*tn, sunvec_f)

    end function RhsFn

end module rhs_mod
