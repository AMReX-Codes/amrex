! -----------------------------------------------------------------
! $Revision$
! $Date$
! -----------------------------------------------------------------
! Programmer(s): David J. Gardner @ LLNL
!                Daniel R. Reynolds @ SMU
! -----------------------------------------------------------------
! LLNS Copyright Start
! Copyright (c) 2014, Lawrence Livermore National Security
! This work was performed under the auspices of the U.S. Department
! of Energy by Lawrence Livermore National Laboratory in part under
! Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
! Produced at the Lawrence Livermore National Laboratory.
! All rights reserved.
! For details, see the LICENSE file.
! LLNS Copyright End
! -----------------------------------------------------------------
! This module implements the Fortran 2003 interface to the SUNDIALS
! serial NVECTOR structure.
! -----------------------------------------------------------------

module fnvector_serial

  !======= Interfaces =========
  interface

     ! -----------------------------------------------------------------
     ! Function : N_VNew_Serial
     ! -----------------------------------------------------------------
     ! This function creates and allocates memory for a serial vector.
     ! -----------------------------------------------------------------

     type(c_ptr) function N_VNew_Serial(vec_length) &
          bind(C,name='N_VNewSerial')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: vec_length
     end function N_VNew_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VNewEmpty_Serial
     ! -----------------------------------------------------------------
     ! This function creates a new serial N_Vector with an empty (NULL)
     ! data array.
     ! -----------------------------------------------------------------

     type(c_ptr) function N_VNewEmpty_Serial(vec_length) &
          bind(C,name='N_VNewEmpty_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: vec_length
     end function N_VNewEmpty_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VMake_Serial
     ! -----------------------------------------------------------------
     ! This function creates and allocates memory for a serial vector
     ! with a user-supplied data array.
     ! -----------------------------------------------------------------

     type(c_ptr) function N_VMake_Serial(length, v_data) &
          bind(C,name='N_VMake_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: length
       real(c_double)         :: v_data(length)
     end function N_VMake_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VCloneVectorArray_Serial
     ! -----------------------------------------------------------------
     ! This function creates an array of 'count' SERIAL vectors by
     ! cloning a given vector w.
     ! -----------------------------------------------------------------

     type(c_ptr) function N_VCloneVectorArray_Serial(count, w) &
          bind(C,name='N_VCloneVectorArray_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: count
       type(c_ptr),    value :: w
     end function N_VCloneVectorArray_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VCloneVectorArrayEmpty_Serial
     ! -----------------------------------------------------------------
     ! This function creates an array of 'count' SERIAL vectors each
     ! with an empty (NULL) data array by cloning w.
     ! -----------------------------------------------------------------

     type(c_ptr) function N_VCloneVectorArrayEmpty_Serial(count, w) &
          bind(C,name='N_VCloneVectorArrayEmpty_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: count
       type(c_ptr),    value :: w
     end function N_VCloneVectorArrayEmpty_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VDestroyVectorArray_Serial
     ! -----------------------------------------------------------------
     ! This function frees an array of SERIAL vectors created with
     ! N_VCloneVectorArray_Serial or N_VCloneVectorArrayEmpty_Serial.
     ! -----------------------------------------------------------------

     subroutine N_VDestroyVectorArray_Serial(vs, count) &
          bind(C,name='N_VDestroyVectorArray_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr)    :: vs
       integer(c_int) :: count
     end subroutine N_VDestroyVectorArray_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VGetLength_Serial
     ! -----------------------------------------------------------------
     ! This function returns number of vector elements.
     ! -----------------------------------------------------------------

     integer(c_long) function N_VGetLength_Serial(v) &
          bind(C,name='N_VGetLength_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
     end function N_VGetLength_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VPrint_Serial
     ! -----------------------------------------------------------------
     ! This function prints the content of a serial vector to stdout.
     ! -----------------------------------------------------------------

     subroutine N_VPrint_Serial(v) &
          bind(C,name='N_VPrint_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
     end subroutine N_VPrint_Serial

     ! =================================================================
     ! serial implementations of various useful vector operations
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function : N_VGetVectorID
     ! -----------------------------------------------------------------
     ! Returns an identifier for the vector type from enumeration
     ! N_Vector_ID.
     ! -----------------------------------------------------------------

     ! integer(c_int) function N_VGetVectorID_Serial(v) &
     !      bind(C,name'N_VGetVectorID_Serial')
     !   use, intrinsic :: iso_c_binding
     !   implicit none
     !   type(c_ptr), value :: vec
     ! end function N_VGetVectorID_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VCloneEmpty_Serial
     ! -----------------------------------------------------------------
     ! Creates a new vector of the same type as an existing vector,
     ! but does not allocate storage.
     ! -----------------------------------------------------------------

     type(c_ptr) function N_VCloneEmpty_Serial(w) &
          bind(C,name='N_VCloneEmpty_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: w
     end function N_VCloneEmpty_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VClone_Serial
     ! -----------------------------------------------------------------
     ! Creates a new vector of the same type as an existing vector.
     ! It does not copy the vector, but rather allocates storage for
     ! the new vector.
     ! -----------------------------------------------------------------

     type(c_ptr) function N_VClone_Serial(w) &
          bind(C,name='N_VClone_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: w
     end function N_VClone_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VDestroy_Serial
     ! -----------------------------------------------------------------
     ! Destroys a vector created with N_VClone_Serial
     ! -----------------------------------------------------------------

     subroutine N_VDestroy_Serial(v) &
          bind(C,name='N_VDestroy_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
     end subroutine N_VDestroy_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VSpace_Serial
     ! -----------------------------------------------------------------
     ! Returns space requirements for one N_Vector (type 'realtype' in
     ! lrw and type 'long int' in liw).
     ! -----------------------------------------------------------------

     subroutine N_VSpace_Serial(v, lrw, liw) &
          bind(C,name='N_VSpace_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
       integer(c_long)    :: lrw
       integer(c_long)    :: liw
     end subroutine N_VSpace_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VGetArrayPointer_Serial
     ! -----------------------------------------------------------------
     ! Returns a pointer to the data component of the given N_Vector.
     !
     ! NOTE: This function assumes that the internal data is stored
     ! as a contiguous 'realtype' array. This routine is only used in
     ! the solver-specific interfaces to the dense and banded linear
     ! solvers, as well as the interfaces to  the banded preconditioners
     ! distributed with SUNDIALS.
     ! -----------------------------------------------------------------

     type(c_ptr) function N_VGetArrayPointer_Serial(vec) &
          bind(C,name='N_VGetArrayPointer_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: vec
     end function N_VGetArrayPointer_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VSetArrayPointer_Serial
     ! -----------------------------------------------------------------
     ! Overwrites the data field in the given N_Vector with a user-supplied
     ! array of type 'realtype'.
     !
     ! NOTE: This function assumes that the internal data is stored
     ! as a contiguous 'realtype' array. This routine is only used in
     ! the interfaces to the dense linear solver.
     ! -----------------------------------------------------------------

     subroutine N_VSetArrayPointer_Serial(v_data, v) &
          bind(C,name='N_VSetArrayPointer_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double)     :: v_data
       type(c_ptr), value :: v
     end subroutine N_VSetArrayPointer_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VLinearSum_Serial
     ! -----------------------------------------------------------------
     ! Performs the operation z = a*x + b*y
     ! -----------------------------------------------------------------

     subroutine N_VLinearSum_Serial(a, x, b, y, z) &
          bind(C,name='N_VLinearSum_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: a
       type(c_ptr),    value :: x
       real(c_double), value :: b
       type(c_ptr),    value :: y
       type(c_ptr),    value :: z
     end subroutine N_VLinearSum_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VConst_Serial
     ! -----------------------------------------------------------------
     ! Performs the operation z[i] = c for i = 0, 1, ..., N-1
     ! -----------------------------------------------------------------

     subroutine N_VConst_Serial(c, z) &
          bind(C,name='N_VConst_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: c
       type(c_ptr),    value :: z
     end subroutine N_VConst_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VProd_Serial
     ! -----------------------------------------------------------------
     ! Performs the operation z[i] = x[i]*y[i] for i = 0, 1, ..., N-1
     ! -----------------------------------------------------------------

     subroutine N_VProd_Serial(x, y, z) &
          bind(C,name='N_VProd_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: y
       type(c_ptr), value :: z
     end subroutine N_VProd_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VDiv_Serial
     ! -----------------------------------------------------------------
     ! Performs the operation z[i] = x[i]/y[i] for i = 0, 1, ..., N-1
     ! -----------------------------------------------------------------

     subroutine N_VDiv_Serial(x, y, z) &
          bind(C,name='N_VDiv_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: y
       type(c_ptr), value :: z
     end subroutine N_VDiv_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VScale_Serial
     ! -----------------------------------------------------------------
     ! Performs the operation z = c*x
     ! -----------------------------------------------------------------

     subroutine N_VScale_Serial(c, x, z) &
          bind(C,name='N_VScale_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: c
       type(c_ptr),    value :: x
       type(c_ptr),    value :: z
     end subroutine N_VScale_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VAbs_Serial
     ! -----------------------------------------------------------------
     ! Performs the operation z[i] = |x[i]| for i = 0, 1, ..., N-1
     ! -----------------------------------------------------------------

     subroutine N_VAbs_Serial(x, z) &
          bind(C,name='N_VAbs_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: z
     end subroutine N_VAbs_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VInv_Serial
     ! -----------------------------------------------------------------
     ! Performs the operation z[i] = 1/x[i] for i = 0, 1, ..., N-1
     !
     ! This routine does not check for division by 0. It should be
     ! called only with an N_Vector x which is guaranteed to have
     ! all non-zero components.
     ! -----------------------------------------------------------------

     subroutine N_VInv_Serial(x, z) &
          bind(C,name='N_VInv_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: z
     end subroutine N_VInv_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VAddConst
     ! -----------------------------------------------------------------
     ! Performs the operation z[i] = x[i] + b   for i = 0, 1, ..., N-1
     ! -----------------------------------------------------------------

     subroutine N_VAddConst_Serial(x, b, z) &
          bind(C,name='N_VAddConst_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: x
       real(c_double), value :: b
       type(c_ptr),    value :: z
     end subroutine N_VAddConst_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VDotProd_Serial
     ! -----------------------------------------------------------------
     ! Returns the dot product of two vectors:
     !     sum (i = 0 to N-1) {x[i]*y[i]}
     ! -----------------------------------------------------------------

     real(c_double) function N_VDotProd_Serial(x, y) &
          bind(C,name='N_VDotProd_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: y
     end function N_VDotProd_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VMaxNorm_Serial
     ! -----------------------------------------------------------------
     ! Returns the maximum norm of x:
     !     max (i = 0 to N-1) ABS(x[i])
     ! -----------------------------------------------------------------

     real(c_double) function N_VMaxNorm_Serial(x) &
          bind(C,name='N_VMaxNorm_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
     end function N_VMaxNorm_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VWrmsNorm_Serial
     ! -----------------------------------------------------------------
     ! Returns the weighted root mean square norm of x with weight
     ! vector w:
     !     sqrt [(sum (i = 0 to N-1) {(x[i]*w[i])^2})/N]
     ! -----------------------------------------------------------------

     real(c_double) function N_VWrmsNorm_Serial(x, w) &
          bind(C,name='N_VWrmsNorm_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: w
     end function N_VWrmsNorm_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VWrmsNormMask_Serial
     ! -----------------------------------------------------------------
     ! Returns the weighted root mean square norm of x with weight
     ! vector w, masked by the elements of id:
     !     sqrt [(sum (i = 0 to N-1) {(x[i]*w[i]*msk[i])^2})/N]
     ! where msk[i] = 1.0 if id[i] > 0 and
     !     msk[i] = 0.0 if id[i] < 0
     ! -----------------------------------------------------------------

     real(c_double) function N_VWrmsNormMask_Serial(x, w, id) &
          bind(C,name='N_VWrmsNormMask_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: w
       type(c_ptr), value :: id
     end function N_VWrmsNormMask_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VMin_Serial
     ! -----------------------------------------------------------------
     ! Returns the smallest element of x:
     !     min (i = 0 to N-1) x[i]
     ! -----------------------------------------------------------------

     real(c_double) function N_VMin_Serial(x) &
          bind(C,name='N_VMin_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
     end function N_VMin_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VWL2Norm_Serial
     ! -----------------------------------------------------------------
     ! Returns the weighted Euclidean L2 norm of x with weight
     ! vector w:
     !     sqrt [(sum (i = 0 to N-1) {(x[i]*w[i])^2})]
     ! -----------------------------------------------------------------

     real(c_double) function N_VWL2Norm_Serial(x, w) &
          bind(C,name='N_VWL2Norm_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: w
     end function N_VWL2Norm_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VL1Norm_Serial
     ! -----------------------------------------------------------------
     ! Returns the L1 norm of x:
     !     sum (i = 0 to N-1) {ABS(x[i])}
     ! -----------------------------------------------------------------

     real(c_double) function N_VL1Norm_Serial(x) &
          bind(C,name='N_VL1Norm_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
     end function N_VL1Norm_Serial

     ! -----------------------------------------------------------------
     ! Subroutine : N_VCompare_Serial
     ! -----------------------------------------------------------------
     ! Performs the operation
     !     z[i] = 1.0 if ABS(x[i]) >= c   i = 0, 1, ..., N-1
     !     0.0 otherwise
     ! -----------------------------------------------------------------

     subroutine N_VCompare_Serial(c, x, z) &
          bind(C,name='N_VCompare_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: c
       type(c_ptr),    value :: x
       type(c_ptr),    value :: z
     end subroutine N_VCompare_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VInvTest_Serial
     ! -----------------------------------------------------------------
     ! Performs the operation z[i] = 1/x[i] with a test for
     !     x[i] == 0.0 before inverting x[i].
     !
     ! This routine returns TRUE if all components of x are non-zero
     ! (successful inversion) and returns FALSE otherwise.
     ! -----------------------------------------------------------------

     integer(c_int) function N_VInvTest_Serial(x, z) &
          bind(C,name='N_VInvTest_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: z
     end function N_VInvTest_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VConstrMask_Serial
     ! -----------------------------------------------------------------
     ! Performs the operation :
     !     m[i] = 1.0 if constraint test fails for x[i]
     !     m[i] = 0.0 if constraint test passes for x[i]
     ! where the constraint tests are as follows:
     !     If c[i] = +2.0, then x[i] must be >  0.0.
     !     If c[i] = +1.0, then x[i] must be >= 0.0.
     !     If c[i] = -1.0, then x[i] must be <= 0.0.
     !     If c[i] = -2.0, then x[i] must be <  0.0.
     ! This routine returns a boolean FALSE if any element failed
     ! the constraint test, TRUE if all passed. It also sets a
     ! mask vector m, with elements equal to 1.0 where the
     ! corresponding constraint test failed, and equal to 0.0
     ! where the constraint test passed.
     !
     ! This routine is specialized in that it is used only for
     ! constraint checking.
     ! -----------------------------------------------------------------

     integer(c_int) function N_VConstrMask_Serial(c, x, m) &
          bind(C,name='N_VConstrMask_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: c
       type(c_ptr), value :: x
       type(c_ptr), value :: m
     end function N_VConstrMask_Serial

     ! -----------------------------------------------------------------
     ! Function : N_VMinQuotient_Serial
     ! -----------------------------------------------------------------
     ! Performs the operation :
     !     minq  = min ( num[i]/denom[i]) over all i such that
     !     denom[i] != 0.
     ! This routine returns the minimum of the quotients obtained
     ! by term-wise dividing num[i] by denom[i]. A zero element
     ! in denom will be skipped. If no such quotients are found,
     ! then the large value BIG_REAL is returned.
     ! -----------------------------------------------------------------

     real(c_double) function N_VMinQuotient_Serial(num, denom) &
          bind(C,name='N_VMinQuotient_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: num
       type(c_ptr), value :: denom
     end function N_VMinQuotient_Serial

  end interface

contains

  ! ================================================================
  ! Helpful N_Vector_Serial Functions / Subroutines
  ! ================================================================

  subroutine N_VGetData_Serial(SUNVec, length, f_array)
    ! ----------------------------------------------------------------
    ! Description: Extracts data array from serial SUNDIALS N_Vector
    ! ----------------------------------------------------------------

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    type(c_ptr)             :: SUNVec
    integer(c_long)         :: length
    real(c_double), pointer :: f_array(:)

    ! C pointer for N_Vector interal data array
    type(c_ptr) :: c_array

    !======= Internals ============

    ! get data pointer from N_Vector
    c_array = N_VGetArrayPointer_Serial(SUNVec)

    ! convert c pointer to f pointer
    call c_f_pointer(c_array, f_array, (/length/))

  end subroutine N_VGetData_Serial

end module fnvector_serial
