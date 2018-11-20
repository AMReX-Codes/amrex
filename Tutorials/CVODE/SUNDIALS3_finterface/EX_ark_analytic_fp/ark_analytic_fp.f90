! ------------------------------------------------------------------
! Programmer(s): David J. Gardner @ LLNL
!                modified by Jean M. Sexton @ LBL
! ------------------------------------------------------------------
! LLNS Copyright Start
! Copyright (c) 2014, Lawrence Livermore National Security
! This work was performed under the auspices of the U.S. Department
! of Energy by Lawrence Livermore National Laboratory in part under
! Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
! Produced at the Lawrence Livermore National Laboratory.
! All rights reserved.
! For details, see the LICENSE file.
! LLNS Copyright End
! ------------------------------------------------------------------
! The following is a simple example problem with an analytical
! solution.
!
!   dy/dt = lamda*y + 1/(1+t^2) - lamda*atan(t)
!
! for t in the interval [0.0, 10.0], with initial condition: y=0.
!
! The stiffness of the problem is directly proportional to the
! value of lamda. The value of lamda should be negative to
! result in a well-posed ODE; for values with magnitude larger
! than 100 the problem becomes quite stiff.
! ------------------------------------------------------------------

module ode_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none

  ! number of equations
  integer(c_long), parameter :: neq = 1

  ! ODE parameters
  double precision, parameter :: lamda = -100.0d0

contains

  ! ----------------------------------------------------------------
  ! RhsFn provides the right hand side function for the
  ! ODE: dy/dt = f(t,y)
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C,name='RhsFn')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fnvector_serial_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn        ! current time
    type(c_ptr), value    :: sunvec_y  ! solution N_Vector
    type(c_ptr), value    :: sunvec_f  ! rhs N_Vector
    type(c_ptr), value    :: user_data ! user-defined data

    ! pointers to data in SUNDAILS vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: fvec(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    call FN_VGetData_Serial(sunvec_y, yvec)
    call FN_VGetData_Serial(sunvec_f, fvec)

    ! fill RHS vector
    fvec(1) = lamda*yvec(1) + 1.0/(1.0+tn*tn) - lamda*atan(tn);

    ! return success
    ierr = 0
    return

  end function RhsFn

end module ode_mod


program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use farkode_mod         ! Fortran interface to ARKODE
  use fnvector_serial_mod ! Fortran interface to serial N_Vector
  use ode_mod             ! ODE functions
  use arkode_interface

  !======= Declarations =========
  implicit none

  ! local variables
  real(c_double) :: tstart     ! initial time
  real(c_double) :: tend       ! final time
  real(c_double) :: rtol, atol ! relative and absolute tolerance
  real(c_double) :: dtout      ! output time interval
  real(c_double) :: tout       ! output time
  real(c_double) :: tcur       ! current time
  integer(c_int) :: imethod, idefault, pq

  integer(c_int) :: ierr       ! error flag from C functions
  integer(c_int) :: nout       ! number of outputs

  integer :: outstep           ! output loop counter

  type(c_ptr) :: sunvec_y      ! sundials vector
  type(c_ptr) :: arkode_mem     ! ARKODE memory

  ! solution vector, neq is set in the ode_functions module
  real(c_double) :: yvec(neq)

  !======= Internals ============

  ! initialize ODE
  tstart = 0.0d0
  tend   = 10.0d0
  tcur   = tstart
  tout   = tstart
  dtout  = 1.0d0
  nout   = ceiling(tend/dtout)

  ! initialize solution vector
  yvec(1) = 0.0d0

  ! create SUNDIALS N_Vector
  sunvec_y = FN_VMake_Serial(neq, yvec)
  if (.not. c_associated(sunvec_y)) print *,'ERROR: sunvec = NULL'

  ! create ARKode memory
  arkode_mem = FARKodeCreate()
  if (.not. c_associated(arkode_mem)) print *,'ERROR: arkode_mem = NULL'

  ! initialize ARKode
!  ierr = FARKodeInit(arkode_mem, c_funloc(RhsFn), c_null_ptr, tstart, sunvec_y)
  ierr = FARKodeInit(arkode_mem, c_null_ptr, c_funloc(RhsFn), tstart, sunvec_y)
  if (ierr /= 0) then
     write(*,*) 'Error in FARKodeInit, ierr = ', ierr, '; halting'
     stop
  end if

  ! Tell ARKODE to use a dense linear solver.
  ierr = FARKDense(arkode_mem, neq)
  if (ierr /= 0) then
     write(*,*) 'Error in FARKDense, ierr = ', ierr, '; halting'
     stop
  end if

  ! set relative and absolute tolerances
  rtol = 1.0d-6
  atol = 1.0d-10

  ierr = FARKodeSStolerances(arkode_mem, rtol, atol)
  if (ierr /= 0) then
     write(*,*) 'Error in FARKodeSStolerances, ierr = ', ierr, '; halting'
     stop
  end if

  imethod = 4
  idefault = 1
  pq = 0
  ierr = FARKodeSetAdaptivityMethod(arkode_mem, imethod, idefault, pq, c_null_ptr)
  
  ! Start time stepping
  print *, '   '
  print *, 'Finished initialization, starting time steps'
  print *, '   '
  print *, '       t           y        '
  print *, '----------------------------'
  print '(2x,2(es12.5,1x))', tcur, yvec(1)
  do outstep = 1,nout

     ! call ARKode
     tout = min(tout + dtout, tend)
     ierr = FARKode(arkode_mem, tout, sunvec_y, tcur, ARK_NORMAL)
     if (ierr /= 0) then
        write(*,*) 'Error in FARKODE, ierr = ', ierr, '; halting'
        stop
     endif

     ! output current solution
     print '(2x,2(es12.5,1x))', tcur, yvec(1)

  enddo

  ! diagnostics output
  call ARKodeStats(arkode_mem)
  call FN_VDestroy_Serial(sunvec_y)

  ! clean up
  call FARKodeFree(arkode_mem)

end program main


! ----------------------------------------------------------------
! ARKodeStats
!
! Print ARKODE statstics to stdandard out
! ----------------------------------------------------------------
subroutine ARKodeStats(arkode_mem)

  !======= Inclusions ===========
  use iso_c_binding
  use farkode_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), intent(in) :: arkode_mem ! solver memory structure

  integer(c_int)  :: ierr ! error flag

  integer(c_long) :: nsteps     ! num steps
  integer(c_long) :: nst_a      ! num steps attempted
  integer(c_long) :: nfe        ! num explicit function evals
  integer(c_long) :: nfi        ! num implicit function evals
  integer(c_long) :: nfevals    ! num function evals
  integer(c_long) :: nlinsetups ! num linear solver setups
  integer(c_long) :: netfails   ! num error test fails

  integer(c_int) :: qlast ! method order in last step
  integer(c_int) :: qcur  ! method order for next step

  real(c_double) :: hinused ! initial step size
  real(c_double) :: hlast   ! last step size
  real(c_double) :: hcur    ! step size for next step
  real(c_double) :: tcur    ! internal time reached

  integer(c_long) :: nniters  ! nonlinear solver iterations
  integer(c_long) :: nncfails ! nonlinear solver fails
  integer(c_long) :: njacevals! number of Jacobian evaluations
  integer(c_long) :: nfeLS    ! number of nonlinear solver rhs evaluations

  !======= Internals ============

  ierr = FARKodeGetNumSteps(arkode_mem, nsteps)
  ierr = FARKodeGetNumStepAttempts(arkode_mem, nst_a)
  ierr = FARKodeGetNumRhsEvals(arkode_mem, nfe, nfi)
  nfevals=nfe+nfi
  ierr = FARKodeGetActualInitStep(arkode_mem, hinused)
  ierr = FARKodeGetLastStep(arkode_mem, hlast)
  ierr = FARKodeGetCurrentStep(arkode_mem, hcur)
  ierr = FARKodeGetCurrentTime(arkode_mem, tcur)
  ierr = FARKodeGetNumLinSolvSetups(arkode_mem, nlinsetups)
  ierr = FARKodeGetNumErrTestFails(arkode_mem, netfails)
  ierr = FARKodeGetNumNonlinSolvIters(arkode_mem, nniters)
  ierr = FARKodeGetNumNonlinSolvConvFails(arkode_mem, nncfails)
  ierr = FARKDlsGetNumJacEvals(arkode_mem, njacevals)
  ierr = FARKDlsGetNumRhsEvals(arkode_mem, nfeLS)

  print *, ' '
  print *, ' General Solver Stats:'
  print '(4x,A,i9)'    ,'Total internal steps taken    =',nsteps
  print '(4x,A,i9)'    ,'Total internal steps attempts =',nst_a
  print '(4x,A,i9)'    ,'Total rhs function calls      =',nfevals
  print '(4x,A,i9)'    ,'Num lin solver setup calls    =',nlinsetups
  print '(4x,A,i9)'    ,'Num error test failures       =',netfails
  print '(4x,A,es12.5)','First internal step size      =',hinused
  print '(4x,A,es12.5)','Last internal step size       =',hlast
  print '(4x,A,es12.5)','Next internal step size       =',hcur
  print '(4x,A,es12.5)','Current internal time         =',tcur
  print '(4x,A,i9)'    ,'Num nonlinear solver iters    =',nniters
  print '(4x,A,i9)'    ,'Num nonlinear solver fails    =',nncfails
  print *, ' '

  return

end subroutine ARKodeStats
