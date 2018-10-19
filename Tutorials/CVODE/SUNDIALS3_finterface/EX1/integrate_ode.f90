subroutine integrate_ode(mf, lo, hi, cvode_meth, cvode_itmeth) bind(C, name="integrate_ode")

  use amrex_error_module
  use rhs_mod
  use ode_params
  use fcvode_mod            ! Fortran interface to CVODE
  use fnvector_serial_mod   ! Fortran interface to serial N_Vector
  use fsunmat_dense_mod     ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod  ! Fortran interface to dense SUNLinearSolver
!  use fnvector_serial
!  use cvode_interface
  use, intrinsic :: iso_c_binding

  implicit none

  integer, intent(in) :: lo(3),hi(3)
  real(c_double), intent(inout) :: mf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
  integer(c_int), intent(in) :: cvode_meth
  integer(c_int), intent(in) :: cvode_itmeth

  ! local variables
  integer i,j,k

  integer(c_int) :: ierr ! CVODE return status
  real(c_double) :: atol, rtol
  real(c_double) :: t0, t1
  real(c_double), pointer :: yvec(:)
  type(c_ptr) :: sunvec_y      ! sundials vector
  type(c_ptr) :: sunmat_A      ! sundials matrix
  type(c_ptr) :: sunlinsol_LS  ! sundials linear solver
  type(c_ptr) :: CVmem

  allocate(yvec(neq))

  ! Allocate a CVODE C struct from the array of variables to be integrated. The resulting C struct points to the same memory as the
  ! Fortran pointer array.
  sunvec_y = FN_VMake_Serial(neq, yvec)
  if (.not. c_associated(sunvec_y)) call amrex_abort("integrate_ode: failed in N_VMake_Serial()")

  CVmem = FCVodeCreate(CV_BDF, CV_NEWTON)
  if (.not. c_associated(CVmem)) call amrex_abort("integrate_ode: failed in FCVodeCreate()")

  ! create a dense matrix
  sunmat_A = FSUNDenseMatrix(neq, neq);
  if (.not. c_associated(sunmat_A)) print *,'ERROR: sunmat = NULL'

  ! create a dense linear solver
  sunlinsol_LS = FSUNDenseLinearSolver(sunvec_y, sunmat_A);
  if (.not. c_associated(sunlinsol_LS)) print *,'ERROR: sunlinsol = NULL'

  t0 = 0.0d0 ! initial time for integration
  t1 = 2.0d0 ! final time for integration

  ! Set up CVODE solver machinery. This allocates a lot of memory. Since we will be solving the same system of ODEs in each cell and
  ! only changing initial conditions among cells, we don't want to allocate and deallocate all this stuff in each cell, because that
  ! will be slow (and unnecessary). So instead, we allocate the solver stuff once outside the (i,j,k) loop, and then within the
  ! (i,j,k) loop we just "re-initialize" the solver, which does not allocate any new data, but rather just sets up new initial
  ! conditions.
  ierr = FCVodeInit(CVmem, c_funloc(RhsFn), t0, sunvec_y)
  if (ierr /= 0) call amrex_abort("integrate_ode: failed in FCVodeInit()")

  ! Set error tolerances tolerances
  atol = 1.0d-10
  rtol = 1.0d-6
  ierr = FCVodeSStolerances(CVmem, rtol, atol)
  if (ierr /= 0) call amrex_abort("integrate_ode: failed in FCVodeSStolerances()")

  ! attach linear solver
  ierr = FCVDlsSetLinearSolver(CVmem, sunlinsol_LS, sunmat_A);
  if (ierr /= 0) then
     write(*,*) 'Error in FCVDlsSetLinearSolver, ierr = ', ierr, '; halting'
     stop
  end if

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           ! Set initial conditions for the ODE for this cell.
           yvec(1) = real(i+j+k, c_double)

           t0 = 0.0d0 ! initial time for integration
           ierr = FCVodeReInit(cvmem, t0, sunvec_y)
           if (ierr /= 0) call amrex_abort("integrate_ode: failed in FCVodeReInit()")

           ierr = FCVode(CVmem, t1, sunvec_y, t0, CV_NORMAL)
           if (ierr /= 0) call amrex_abort("integrate_ode: failed in FCVode()")

           ! Copy the solution of the ODE to the MultiFab.
           mf(i,j,k) = yvec(1)

        end do
     end do
  end do

  ! Free memory
  call FN_VDestroy_Serial(sunvec_y)
  call FCVodeFree(cvmem)

  deallocate(yvec)

end subroutine integrate_ode
