! Right hand side of ODE for CVODE to solve.

module rhs_mod
  implicit none
  contains

    attributes(device) subroutine RhsFn(tn, yvec, fvec, neq) &
           bind(C,name='RhsFn')

      use, intrinsic :: iso_c_binding

      implicit none

      real(c_double), value :: tn
      integer(c_int), value :: neq
!      type(c_ptr), value    :: sunvec_y
!      type(c_ptr), value    :: sunvec_f
!      type(c_ptr), value    :: user_data

      ! pointers to data in SUNDAILS vectors
      real(c_double) :: yvec
      real(c_double) :: fvec

      fvec=2.0*tn
      ! get data arrays from SUNDIALS vectors
!      call N_VGetData_Serial(sunvec_f, neq, fvec)
!      call FN_VConst_Cuda(2.0*tn, sunvec_f)

!      ierr=0
    end subroutine RhsFn

end module rhs_mod
