module ode_params
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_long), parameter :: neq = 3
end module ode_params

! Right hand side of ODE for CVODE to solve.

module rhs_mod
  implicit none

  contains

    integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
           result(ierr) bind(C,name='RhsFn')

      use, intrinsic :: iso_c_binding
      use fnvector_serial
      use cvode_interface
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
      call N_VGetData_Serial(sunvec_f, neq, fvec)
      call N_VGetData_Serial(sunvec_y, neq, yvec)

      fvec(1) = -0.04D0 * yvec(1) + 1.0D4 * yvec(2) * yvec(3)
      fvec(3) = 3.0D7 * yvec(2) * yvec(2)
      fvec(2) = -fvec(1) - fvec(3)

      ierr = 0
      return
    end function RhsFn

end module rhs_mod

module jac_mod

  implicit none

  contains

    integer(c_int) function JacFn(N, tn, sunvec_y, sunvec_f, sunDlsMat_J, &
           user_data, tmp1, tmp2, tmp3) result(ierr) bind(C,name='JacFn')
        ! ----------------------------------------------------------------
        ! Description: JacFn provides the Jacobian function J(t,y) = df/dy
        !
        ! Inputs:
        !           N - problem size
        !          tn - current time
        !    sunvec_y - C pointer to current solution N_Vector
        !    sunvec_f - C pointer to right-hand side function N_Vector
        ! sunDlsMat_J - C pointer to dense Jacobian matrix
        !   user_data - C pointer to user-defined data
        !        tmp1 - temporary workspace N_Vector
        !        tmp2 - temporary workspace N_Vector
        !        tmp3 - temporary workspace N_Vector
        !
        ! Output:
        !        ierr - return flag: 
        !                0 = success, 
        !                1 = recoverable error, 
        !               -1 = non-recoverable error
        ! ----------------------------------------------------------------
    
        !======= Inclusions ===========
        use, intrinsic :: iso_c_binding
        use fnvector_serial
        use sundials_fdlsmat
        use ode_params
    
        !======= Declarations =========
        implicit none
        
        ! calling variables
        integer(c_long), value :: N
        real(c_double),  value :: tn
        type(c_ptr),     value :: sunvec_y
        type(c_ptr),     value :: sunvec_f
        type(sundlsmat)        :: sundlsmat_J
        type(c_ptr),     value :: user_data
        type(c_ptr),     value :: tmp1, tmp2, tmp3
    
        ! pointers to data in SUNDAILS vector and matrix
        real(c_double), pointer :: yvec(:)
        real(c_double), pointer :: Jmat(:,:)
    
        real(c_double) :: y1, y2, y3
    
        !======= Internals ============
        
        ! get data array from SUNDIALS vector
        call N_VGetData_Serial(sunvec_y, neq, yvec)
    
        ! get data array from SUNDIALS matrix
        call sundlsmat_GetData_Dense(sundlsmat_J, Jmat)
    
          y1 = yvec(1)
          y2 = yvec(2)
          y3 = yvec(3)
          jmat(1,1) = -0.04d0
          jmat(1,2) = 1.0d4 * y3
          jmat(1,3) = 1.0d4 * y2
          jmat(2,1) =  0.04d0
          jmat(2,2) = -1.0d4 * y3 - 6.0d7 * y2
          jmat(2,3) = -1.0d4 * y2
          jmat(3,2) = 6.0d7 * y2
    
          ierr = 0
          return
    end function JacFn

end module jac_mod
