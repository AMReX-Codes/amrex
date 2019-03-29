module arkode_interface
  use farkode_mod
  use farkode_arkstep_mod
  use fsunmatrix_dense_mod
  use fsunlinsol_dense_mod
  contains
  integer(c_int) function FARKDense(arkode_mem, N) 
    use, intrinsic :: iso_c_binding
    use fnvector_serial_mod
    implicit none 
    type(c_ptr),     value :: arkode_mem
    integer(c_long), value :: N
    type(c_ptr)            :: sunmat_A
    type(c_ptr)            :: sunlinsol_LS
    type(c_ptr)            :: sunvec_y
    integer(c_int)         :: ierr

    sunvec_y = FN_VNewEmpty_Serial(N)
    sunmat_A = FSUNDenseMatrix(N, N)
    sunlinsol_LS = FSUNDenseLinearSolver(sunvec_y, sunmat_A)
    ierr = FARKStepSetLinearSolver(arkode_mem, sunlinsol_LS, sunmat_A)
  end function FARKDense
end module arkode_interface
