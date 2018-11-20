subroutine integrate_ode_with_jac(mf, lo, hi, cvode_meth, cvode_itmeth) bind(C, name="integrate_ode_with_jac")

  use amrex_error_module
  use rhs_mod
  use jac_mod
  use ode_params
  use fnvector_serial
  use cvode_interface
  use, intrinsic :: iso_c_binding

  implicit none

  integer, intent(in) :: lo(3),hi(3)
  real(c_double), intent(inout) :: mf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), neq)
  integer(c_int), intent(in) :: cvode_meth
  integer(c_int), intent(in) :: cvode_itmeth

  ! local variables
  integer i,j,k

  integer(c_int) :: ierr ! CVODE return status
  real(c_double) :: atol(neq), rtol
  real(c_double) :: t0, t1
  real(c_double), pointer :: yvec(:)
  type(c_ptr) :: sunvec_y
  type(c_ptr) :: CVmem
  type(c_ptr) :: atol_cptr
  integer(c_long), parameter :: mxsteps = 2000

  allocate(yvec(neq))

  ! Allocate a CVODE C struct from the array of variables to be integrated. The resulting C struct points to the same memory as the
  ! Fortran pointer array.
  sunvec_y = N_VMake_Serial(neq, yvec)
  if (.not. c_associated(sunvec_y)) call amrex_abort("integrate_ode_with_jac: failed in N_VMake_Serial()")

  CVmem = FCVodeCreate(CV_BDF, CV_NEWTON)
  if (.not. c_associated(CVmem)) call amrex_abort("integrate_ode_with_jac: failed in FCVodeCreate()")

  t0 = 0.0d0 ! initial time for integration
  t1 = 4.0d10 ! final time for integration

  ! Set up CVODE solver machinery. This allocates a lot of memory. Since we will be solving the same system of ODEs in each cell and
  ! only changing initial conditions among cells, we don't want to allocate and deallocate all this stuff in each cell, because that
  ! will be slow (and unnecessary). So instead, we allocate the solver stuff once outside the (i,j,k) loop, and then within the
  ! (i,j,k) loop we just "re-initialize" the solver, which does not allocate any new data, but rather just sets up new initial
  ! conditions.
  ierr = FCVodeInit(CVmem, c_funloc(RhsFn), t0, sunvec_y)
  if (ierr /= 0) call amrex_abort("integrate_ode: failed in FCVodeInit()")

  ! Set error tolerances tolerances
  atol(1) = 1.0d-8
  atol(2) = 1.0d-14
  atol(3) = 1.0d-6
  atol_cptr = N_VMake_Serial(neq, atol)
  rtol = 1.0d-4
  ierr = FCVodeSVtolerances(CVmem, rtol, atol_cptr)
  if (ierr /= 0) call amrex_abort("integrate_ode: failed in FCVodeSVtolerances()")

  ! Tell CVODE to use a dense linear solver.
  ierr = FCVDense(CVmem, neq)
  if (ierr /= 0) call amrex_abort("integrate_ode: failed in FCVDense()")

  ! set Jacobian routine
  ierr = FCVDlsSetDenseJacFn(CVmem, c_funloc(JacFn))
  if (ierr /= 0) call amrex_abort("integrate_ode_with_jac: failed in FCVDlsSetDenseJacFn()")

  ! By default, CVode will quit if it cannot reached the final time within 500 steps. For this particular problem, we need more than
  ! 500 steps, so increase the default to 2000.
  ierr = FCVodeSetMaxNumSteps(CVmem, mxsteps)
  if (ierr /= 0) call amrex_abort("integrate_ode_with_jac: failed in FCVodeSetMaxNumSteps()")

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           call N_VGetData_Serial(sunvec_y, neq, yvec)
           ! Set initial conditions for the ODE for this cell. We will solve the same system of ODEs with the same initial
           ! conditions in every cell.
           yvec(1) = 1.0d0
           yvec(2) = 0.0d0
           yvec(3) = 0.0d0

           t0 = 0.0d0 ! initial time for integration
           ierr = FCVodeReInit(cvmem, t0, sunvec_y)
           if (ierr /= 0) call amrex_abort("integrate_ode_with_jac: failed in FCVodeReInit()")

           ierr = FCVode(CVmem, t1, sunvec_y, t0, CV_NORMAL)
           if (ierr /= 0) call amrex_abort("integrate_ode_with_jac: failed in FCVode()")

           ! Copy the solution of the ODE to the MultiFab.
           mf(i,j,k,1:neq) = yvec(1:neq)

        end do
     end do
  end do

  ! Free memory
  call N_VDestroy_Serial(sunvec_y)
  call FCVodeFree(cvmem)

  deallocate(yvec)

end subroutine integrate_ode_with_jac
