! ------------------------------------------------------------------
! $Revision$
! $Date$
! ------------------------------------------------------------------
! Programmer(s): David J. Gardner @ LLNL
!                Daniel R. Reynolds @ SMU
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
! This file contains Fortran modules for interfacing with the main
! CVODE integrator using the ISO_C_BINDING module.
! ------------------------------------------------------------------
! CVODE is used to solve numerically the ordinary initial value
! problem:
!
!                 y' = f(t,y),
!                 y(t0) = y0,
!
! where t0, y0 in R^N, and f: R x R^N -> R^N are given.
! ------------------------------------------------------------------

module cvode_interface

  ! =================================================================
  !              C V O D E     C O N S T A N T S
  ! =================================================================
  use, intrinsic :: iso_c_binding, only : c_int

  ! -----------------------------------------------------------------
  ! Enumerations for inputs to CVodeCreate and CVode.
  ! -----------------------------------------------------------------

  ! lmm
  integer(c_int), parameter :: CV_ADAMS = 1
  integer(c_int), parameter :: CV_BDF   = 2

  ! iter
  integer(c_int), parameter :: CV_FUNCTIONAL = 1
  integer(c_int), parameter :: CV_NEWTON     = 2

  ! itask
  integer(c_int), parameter :: CV_NORMAL   = 1
  integer(c_int), parameter :: CV_ONE_STEP = 2

  ! -----------------------------------------------------------------
  ! CVODE return flags
  ! -----------------------------------------------------------------

  integer(c_int), parameter :: CV_SUCCESS           =   0
  integer(c_int), parameter :: CV_TSTOP_RETURN      =   1
  integer(c_int), parameter :: CV_ROOT_RETURN       =   2

  integer(c_int), parameter :: CV_WARNING           =  99

  integer(c_int), parameter :: CV_TOO_MUCH_WORK     =  -1
  integer(c_int), parameter :: CV_TOO_MUCH_ACC      =  -2
  integer(c_int), parameter :: CV_ERR_FAILURE       =  -3
  integer(c_int), parameter :: CV_CONV_FAILURE      =  -4

  integer(c_int), parameter :: CV_LINIT_FAIL        =  -5
  integer(c_int), parameter :: CV_LSETUP_FAIL       =  -6
  integer(c_int), parameter :: CV_LSOLVE_FAIL       =  -7
  integer(c_int), parameter :: CV_RHSFUNC_FAIL      =  -8
  integer(c_int), parameter :: CV_FIRST_RHSFUNC_ERR =  -9
  integer(c_int), parameter :: CV_REPTD_RHSFUNC_ERR =  -10
  integer(c_int), parameter :: CV_UNREC_RHSFUNC_ERR =  -11
  integer(c_int), parameter :: CV_RTFUNC_FAIL       =  -12

  integer(c_int), parameter :: CV_MEM_FAIL          =  -20
  integer(c_int), parameter :: CV_MEM_NULL          =  -21
  integer(c_int), parameter :: CV_ILL_INPUT         =  -22
  integer(c_int), parameter :: CV_NO_MALLOC         =  -23
  integer(c_int), parameter :: CV_BAD_K             =  -24
  integer(c_int), parameter :: CV_BAD_T             =  -25
  integer(c_int), parameter :: CV_BAD_DKY           =  -26
  integer(c_int), parameter :: CV_TOO_CLOSE         =  -27

  ! =================================================================
  !              C V D I A G     C O N S T A N T S
  ! =================================================================

  ! -----------------------------------------------------------------
  ! CVDIAG return values
  ! -----------------------------------------------------------------

  integer(c_int), parameter :: CVDIAG_SUCCESS         =  0
  integer(c_int), parameter :: CVDIAG_MEM_NULL        = -1
  integer(c_int), parameter :: CVDIAG_LMEM_NULL       = -2
  integer(c_int), parameter :: CVDIAG_ILL_INPUT       = -3
  integer(c_int), parameter :: CVDIAG_MEM_FAIL        = -4

  ! Additional last_flag values
  integer(c_int), parameter :: CVDIAG_INV_FAIL        = -5
  integer(c_int), parameter :: CVDIAG_RHSFUNC_UNRECVR = -6
  integer(c_int), parameter :: CVDIAG_RHSFUNC_RECVR   = -7

  ! =================================================================
  !              C V D I R E C T     C O N S T A N T S
  ! =================================================================

  ! -----------------------------------------------------------------
  ! CVDLS return values
  ! -----------------------------------------------------------------

  integer(c_int), parameter :: CVDLS_SUCCESS         =  0
  integer(c_int), parameter :: CVDLS_MEM_NULL        = -1
  integer(c_int), parameter :: CVDLS_LMEM_NULL       = -2
  integer(c_int), parameter :: CVDLS_ILL_INPUT       = -3
  integer(c_int), parameter :: CVDLS_MEM_FAIL        = -4

  ! Additional last_flag values
  integer(c_int), parameter :: CVDLS_JACFUNC_UNRECVR = -5
  integer(c_int), parameter :: CVDLS_JACFUNC_RECVR   = -6

  ! =================================================================
  !              C V S P A R S E    C O N S T A N T S
  ! =================================================================

  ! -----------------------------------------------------------------
  ! CVSLS return values
  ! -----------------------------------------------------------------

  integer(c_int), parameter :: CVSLS_SUCCESS           =  0
  integer(c_int), parameter :: CVSLS_MEM_NULL          = -1
  integer(c_int), parameter :: CVSLS_LMEM_NULL         = -2
  integer(c_int), parameter :: CVSLS_ILL_INPUT         = -3
  integer(c_int), parameter :: CVSLS_MEM_FAIL          = -4
  integer(c_int), parameter :: CVSLS_JAC_NOSET         = -5
  integer(c_int), parameter :: CVSLS_PACKAGE_FAIL      = -6

  ! Additional last_flag values
  integer(c_int), parameter :: CVSLS_JACFUNC_UNRECVR   = -7
  integer(c_int), parameter :: CVSLS_JACFUNC_RECVR     = -8

  ! Return values for the adjoint module
  integer(c_int), parameter :: CVSLS_NO_ADJ     = -101
  integer(c_int), parameter :: CVSLS_LMEMB_NULL = -102

  ! =================================================================
  !              C V S P I L S    C O N S T A N T S
  ! =================================================================

  ! -----------------------------------------------------------------
  ! CVSPILS return values
  ! -----------------------------------------------------------------

  integer(c_int), parameter :: CVSPILS_SUCCESS   =  0
  integer(c_int), parameter :: CVSPILS_MEM_NULL  = -1
  integer(c_int), parameter :: CVSPILS_LMEM_NULL = -2
  integer(c_int), parameter :: CVSPILS_ILL_INPUT = -3
  integer(c_int), parameter :: CVSPILS_MEM_FAIL  = -4
  integer(c_int), parameter :: CVSPILS_PMEM_NULL = -5

  ! =================================================================
  !          U S E R - C A L L A B L E   R O U T I N E S
  ! =================================================================

  interface
     ! =================================================================
     ! Interfaces from cvode.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function : FCVodeCreate
     ! -----------------------------------------------------------------
     ! FCVodeCreate creates an internal memory block for a problem to
     ! be solved by CVODE.
     !
     ! lmm   is the type of linear multistep method to be used.
     !       The legal values are CV_ADAMS and CV_BDF (see previous
     !       description).
     !
     ! iter  is the type of iteration used to solve the nonlinear
     !       system that arises during each internal time step.
     !       The legal values are CV_FUNCTIONAL and CV_NEWTON.
     !
     ! If successful, FCVodeCreate returns a pointer to initialized
     ! problem memory. This pointer should be passed to FCVodeInit.
     ! If an initialization error occurs, FCVodeCreate prints an error
     ! message to standard err and returns NULL.
     ! -----------------------------------------------------------------

     type(c_ptr) function FCVodeCreate(lmm, iter) &
          bind(C,name='CVodeCreate')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: lmm
       integer(c_int), value :: iter
     end function FCVodeCreate

     ! -----------------------------------------------------------------
     ! Integrator optional input specification functions
     ! -----------------------------------------------------------------
     ! The following functions can be called to set optional inputs
     ! to values other than the defaults given below:
     !
     ! Function                 |  Optional input / [ default value ]
     ! -----------------------------------------------------------------
     !                          |
     ! FCVodeSetErrHandlerFn    | user-provided ErrHandler function.
     !                          | [internal]
     !                          |
     ! FCVodeSetErrFile         | the file pointer for an error file
     !                          | where all CVODE warning and error
     !                          | messages will be written if the default
     !                          | internal error handling function is used.
     !                          | This parameter can be stdout (standard
     !                          | output), stderr (standard error), or a
     !                          | file pointer (corresponding to a user
     !                          | error file opened for writing) returned
     !                          | by fopen.
     !                          | If not called, then all messages will
     !                          | be written to the standard error stream.
     !                          | [stderr]
     !                          |
     ! FCVodeSetUserData        | a pointer to user data that will be
     !                          | passed to the user's f function every
     !                          | time f is called.
     !                          | [NULL]
     !                          |
     ! FCVodeSetMaxOrd          | maximum lmm order to be used by the
     !                          | solver.
     !                          | [12 for Adams , 5 for BDF]
     !                          |
     ! FCVodeSetMaxNumSteps     | maximum number of internal steps to be
     !                          | taken by the solver in its attempt to
     !                          | reach tout.
     !                          | [500]
     !                          |
     ! FCVodeSetMaxHnilWarns    | maximum number of warning messages
     !                          | issued by the solver that t+h==t on the
     !                          | next internal step. A value of -1 means
     !                          | no such messages are issued.
     !                          | [10]
     !                          |
     ! FCVodeSetStabLimDet      | flag to turn on/off stability limit
     !                          | detection (TRUE = on, FALSE = off).
     !                          | When BDF is used and order is 3 or
     !                          | greater, CVsldet is called to detect
     !                          | stability limit.  If limit is detected,
     !                          | the order is reduced.
     !                          | [FALSE]
     !                          |
     ! FCVodeSetInitStep        | initial step size.
     !                          | [estimated by CVODE]
     !                          |
     ! FCVodeSetMinStep         | minimum absolute value of step size
     !                          | allowed.
     !                          | [0.0]
     !                          |
     ! FCVodeSetMaxStep         | maximum absolute value of step size
     !                          | allowed.
     !                          | [infinity]
     !                          |
     ! FCVodeSetStopTime        | the independent variable value past
     !                          | which the solution is not to proceed.
     !                          | [infinity]
     !                          |
     ! FCVodeSetMaxErrTestFails | Maximum number of error test failures
     !                          | in attempting one step.
     !                          | [7]
     !                          |
     ! FCVodeSetMaxNonlinIters  | Maximum number of nonlinear solver
     !                          | iterations at one solution.
     !                          | [3]
     !                          |
     ! FCVodeSetMaxConvFails    | Maximum number of convergence failures
     !                          | allowed in attempting one step.
     !                          | [10]
     !                          |
     ! FCVodeSetNonlinConvCoef  | Coefficient in the nonlinear
     !                          | convergence test.
     !                          | [0.1]
     !                          |
     ! -----------------------------------------------------------------
     !                          |
     ! FCVodeSetIterType        | Changes the current nonlinear iteration
     !                          | type.
     !                          | [set by FCVodecreate]
     !                          |
     ! -----------------------------------------------------------------
     !                             |
     ! FCVodeSetRootDirection      | Specifies the direction of zero
     !                             | crossings to be monitored
     !                             | [both directions]
     !                             |
     ! FCVodeSetNoInactiveRootWarn | disable warning about possible
     !                             | g==0 at beginning of integration
     !                             |
     ! -----------------------------------------------------------------

     ! -----------------------------------------------------------------
     ! Return flag:
     !   CV_SUCCESS   if successful
     !   CV_MEM_NULL  if the cvode memory is NULL
     !   CV_ILL_INPUT if an argument has an illegal value
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeSetErrHandlerFn(cvode_mem, ehfun, eh_data) &
          bind(C,name='CVodeSetErrHandlerFn')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: ehfun
       type(c_ptr),    value :: eh_data
     end function FCVodeSetErrHandlerFn

     ! >>> NOT CURRENTLY IMPLEMENTED IN FORTRAN INTERFACE
     ! int CVodeSetErrFile(void *cvode_mem, FILE *errfp);

     integer(c_int) function FCVodeSetUserData(cvode_mem, user_data) &
          bind(C,name='CVodeSetUserData')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       type(c_ptr), value :: user_data
     end function FCVodeSetUserData

     integer(c_int) function FCVodeSetMaxOrd(cvode_mem, maxord) &
          bind(C,name='CVodeSetMaxOrd')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: maxord
     end function FCVodeSetMaxOrd

     integer(c_int) function FCVodeSetMaxNumSteps(cvode_mem, mxsteps) &
          bind(C,name='CVodeSetMaxNumSteps')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: cvode_mem
       integer(c_long), value :: mxsteps
     end function FCVodeSetMaxNumSteps

     integer(c_int) function FCVodeSetMaxHnilWarns(cvode_mem, mxhnil) &
          bind(C,name='CVodeSetMaxHnilWarns')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: mxhnil
     end function FCVodeSetMaxHnilWarns

     integer(c_int) function FCVodeSetStabLimDet(cvode_mem, stldet) &
          bind(C,name='CVodeSetStabLimDet')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: stldet
     end function FCVodeSetStabLimDet

     integer(c_int) function FCVodeSetInitStep(cvode_mem, hin) &
          bind(C,name='CVodeSetInitStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: hin
     end function FCVodeSetInitStep

     integer(c_int) function FCVodeSetMinStep(cvode_mem, hmin) &
          bind(C,name='CVodeSetMinStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: hmin
     end function FCVodeSetMinStep

     integer(c_int) function FCVodeSetMaxStep(cvode_mem, hmax) &
          bind(C,name='CVodeSetMaxStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: hmax
     end function FCVodeSetMaxStep

     integer(c_int) function FCVodeSetStopTime(cvode_mem, tstop) &
          bind(C,name='CVodeSetStopTime')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: tstop
     end function FCVodeSetStopTime

     integer(c_int) function FCVodeSetMaxErrTestFails(cvode_mem, maxnef) &
          bind(C,name='CVodeSetMaxErrTestFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: maxnef
     end function FCVodeSetMaxErrTestFails

     integer(c_int) function FCVodeSetMaxNonlinIters(cvode_mem, maxcor) &
          bind(C,name='CVodeSetMaxNonlinIters')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: maxcor
     end function FCVodeSetMaxNonlinIters

     integer(c_int) function FCVodeSetMaxConvFails(cvode_mem, maxncf) &
          bind(C,name='CVodeSetMaxConvFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: maxncf
     end function FCVodeSetMaxConvFails

     integer(c_int) function FCVodeSetNonlinConvCoef(cvode_mem, nlscoef) &
          bind(C,name='CVodeSetNonlinConvCoef')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: nlscoef
     end function FCVodeSetNonlinConvCoef

     integer(c_int) function FCVodeSetIterType(cvode_mem, iter) &
          bind(C,name='CVodeSetIterType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: iter
     end function FCVodeSetIterType

     integer(c_int) function FCVodeSetRootDirection(cvode_mem, rootdir) &
          bind(C,name='CVodeSetRootDirection')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       real(c_double)     :: rootdir
     end function FCVodeSetRootDirection

     integer(c_int) function FCVodeSetNoInactiveRootWarn(cvode_mem) &
       bind(C,name='CVodeSetNoInactiveRootWarn')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
     end function FCVodeSetNoInactiveRootWarn

     ! -----------------------------------------------------------------
     ! Function : FCVodeInit
     ! -----------------------------------------------------------------
     ! FCVodeInit allocates and initializes memory for a problem to
     ! to be solved by CVODE.
     !
     ! cvode_mem is pointer to CVODE memory returned by FCVodeCreate.
     !
     ! f       is the name of the C function defining the right-hand
     !         side function in y' = f(t,y).
     !
     ! t0      is the initial value of t.
     !
     ! y0      is the initial condition vector y(t0).
     !
     ! Return flag:
     !  CV_SUCCESS if successful
     !  CV_MEM_NULL if the cvode memory was NULL
     !  CV_MEM_FAIL if a memory allocation failed
     !  CV_ILL_INPUT f an argument has an illegal value.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeInit(cvode_mem, f, t0, y0) &
          bind(C,name='CVodeInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: f
       real(c_double), value :: t0
       type(c_ptr),    value :: y0
     end function FCVodeInit

     ! -----------------------------------------------------------------
     ! Function : FCVodeReInit
     ! -----------------------------------------------------------------
     ! FCVodeReInit re-initializes CVode for the solution of a problem,
     ! where a prior call to FCVodeInit has been made with the same
     ! problem size N. FCVodeReInit performs the same input checking
     ! and initializations that FCVodeInit does.
     ! But it does no memory allocation, assuming that the existing
     ! internal memory is sufficient for the new problem.
     !
     ! The use of FCVodeReInit requires that the maximum method order,
     ! maxord, is no larger for the new problem than for the problem
     ! specified in the last call to FCVodeInit.  This condition is
     ! automatically fulfilled if the multistep method parameter lmm
     ! is unchanged (or changed from CV_ADAMS to CV_BDF) and the default
     ! value for maxord is specified.
     !
     ! All of the arguments to FCVodeReInit have names and meanings
     ! identical to those of FCVodeInit.
     !
     ! The return value of FCVodeReInit is equal to CV_SUCCESS = 0 if
     ! there were no errors; otherwise it is a negative int equal to:
     !   CV_MEM_NULL      indicating cvode_mem was NULL (i.e.,
     !                    FCVodeCreate has not been called).
     !   CV_NO_MALLOC     indicating that cvode_mem has not been
     !                    allocated (i.e., FCVodeInit has not been
     !                    called).
     !   CV_ILL_INPUT     indicating an input argument was illegal
     !                    (including an attempt to increase maxord).
     ! In case of an error return, an error message is also printed.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeReInit(cvode_mem, t0, y0) &
          bind(C,name='CVodeReInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: t0
       type(c_ptr),    value :: y0
     end function FCVodeReInit

     ! -----------------------------------------------------------------
     ! Functions : FCVodeSStolerances
     !             FCVodeSVtolerances
     !             FCVodeWFtolerances
     ! -----------------------------------------------------------------
     !
     ! These functions specify the integration tolerances. One of them
     ! MUST be called before the first call to CVode.
     !
     ! FCVodeSStolerances specifies scalar relative and absolute tolerances.
     ! FCVodeSVtolerances specifies scalar relative tolerance and a vector
     !   absolute tolerance (a potentially different absolute tolerance
     !   for each vector component).
     ! FCVodeWFtolerances specifies a user-provides function (of type CVEwtFn)
     !   which will be called to set the error weight vector.
     !
     ! The tolerances reltol and abstol define a vector of error weights,
     ! ewt, with components
     !   ewt[i] = 1/(reltol*abs(y[i]) + abstol)      (in the SS case), or
     !   ewt[i] = 1/(reltol*abs(y[i]) + abstol[i])   (in the SV case).
     ! This vector is used in all error and convergence tests, which
     ! use a weighted RMS norm on all error-like vectors v:
     !    WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*ewt[i])^2 ),
     ! where N is the problem dimension.
     !
     ! The return value of these functions is equal to CV_SUCCESS = 0 if
     ! there were no errors; otherwise it is a negative int equal to:
     !   CV_MEM_NULL      indicating cvode_mem was NULL (i.e.,
     !                    FCVodeCreate has not been called).
     !   CV_NO_MALLOC     indicating that cvode_mem has not been
     !                    allocated (i.e., FCVodeInit has not been
     !                    called).
     !   CV_ILL_INPUT     indicating an input argument was illegal
     !                    (e.g. a negative tolerance)
     ! In case of an error return, an error message is also printed.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeSStolerances(cvode_mem, reltol, abstol) &
          bind(C,name='CVodeSStolerances')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: reltol
       real(c_double), value :: abstol
     end function FCVodeSStolerances

     integer(c_int) function FCVodeSVtolerances(cvode_mem, reltol, abstol) &
          bind(C,name='CVodeSVtolerances')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: reltol
       type(c_ptr),    value :: abstol
     end function FCVodeSVtolerances

     integer(c_int) function FCVodeWFtolerances(cvode_mem, efun) &
          bind(C,name='CVodeWFtolerances')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: efun
     end function FCVodeWFtolerances

     ! -----------------------------------------------------------------
     ! Function : FCVodeRootInit
     ! -----------------------------------------------------------------
     ! FCVodeRootInit initializes a rootfinding problem to be solved
     ! during the integration of the ODE system.  It must be called
     ! after FCVodeCreate, and before CVode.  The arguments are:
     !
     ! cvode_mem = pointer to CVODE memory returned by FCVodeCreate.
     !
     ! nrtfn     = number of functions g_i, an int >= 0.
     !
     ! g         = name of user-supplied function, of type CVRootFn,
     !             defining the functions g_i whose roots are sought.
     !
     ! If a new problem is to be solved with a call to FCVodeReInit,
     ! where the new problem has no root functions but the prior one
     ! did, then call FCVodeRootInit with nrtfn = 0.
     !
     ! The return value of FCVodeRootInit is CV_SUCCESS = 0 if there were
     ! no errors; otherwise it is a negative int equal to:
     !   CV_MEM_NULL     indicating cvode_mem was NULL, or
     !   CV_MEM_FAIL     indicating a memory allocation failed.
     !                    (including an attempt to increase maxord).
     !   CV_ILL_INPUT   indicating nrtfn > 0 but g = NULL.
     ! In case of an error return, an error message is also printed.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeRootInit(cvode_mem, nrtfn, g) &
          bind(C,name='CVodeRootInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: nrtfn
       type(c_funptr), value :: g
     end function FCVodeRootInit

     ! -----------------------------------------------------------------
     ! Function : FCVode
     ! -----------------------------------------------------------------
     ! FCVode integrates the ODE over an interval in t.
     ! If itask is CV_NORMAL, then the solver integrates from its
     ! current internal t value to a point at or beyond tout, then
     ! interpolates to t = tout and returns y(tout) in the user-
     ! allocated vector yout. If itask is CV_ONE_STEP, then the solver
     ! takes one internal time step and returns in yout the value of
     ! y at the new internal time. In this case, tout is used only
     ! during the first call to CVode to determine the direction of
     ! integration and the rough scale of the t variable. If tstop is
     ! enabled (through a call to FCVodeSetStopTime), then CVode returns
     ! the solution at tstop. Once the integrator returns at a tstop
     ! time, any future testing for tstop is disabled (and can be
     ! reenabled only though a new call to FCVodeSetStopTime).
     ! The time reached by the solver is placed in (*tret). The
     ! user is responsible for allocating the memory for this value.
     !
     ! cvode_mem is the pointer to CVODE memory returned by
     !           FCVodeCreate.
     !
     ! tout  is the next time at which a computed solution is desired.
     !
     ! yout  is the computed solution vector. In CV_NORMAL mode with no
     !       errors and no roots found, yout=y(tout).
     !
     ! tret  is a pointer to a real location. CVode sets (*tret) to
     !       the time reached by the solver and returns
     !       yout=y(*tret).
     !
     ! itask is CV_NORMAL or CV_ONE_STEP. These two modes are described above.
     !
     ! Here is a brief description of each return value:
     !
     ! CV_SUCCESS:      CVode succeeded and no roots were found.
     !
     ! CV_ROOT_RETURN:  CVode succeeded, and found one or more roots.
     !                  If nrtfn > 1, call FCVodeGetRootInfo to see
     !                  which g_i were found to have a root at (*tret).
     !
     ! CV_TSTOP_RETURN: CVode succeeded and returned at tstop.
     !
     ! CV_MEM_NULL:     The cvode_mem argument was NULL.
     !
     ! CV_NO_MALLOC:    cvode_mem was not allocated.
     !
     ! CV_ILL_INPUT:    One of the inputs to FCVode is illegal. This
     !                  includes the situation when a component of the
     !                  error weight vectors becomes < 0 during
     !                  internal time-stepping.  It also includes the
     !                  situation where a root of one of the root
     !                  functions was found both at t0 and very near t0.
     !                  The ILL_INPUT flag will also be returned if the
     !                  linear solver routine CV--- (called by the user
     !                  after calling FCVodeCreate) failed to set one of
     !                  the linear solver-related fields in cvode_mem or
     !                  if the linear solver's init routine failed. In
     !                  any case, the user should see the printed
     !                  error message for more details.
     !
     ! CV_TOO_MUCH_WORK: The solver took mxstep internal steps but
     !                  could not reach tout. The default value for
     !                  mxstep is MXSTEP_DEFAULT = 500.
     !
     ! CV_TOO_MUCH_ACC: The solver could not satisfy the accuracy
     !                  demanded by the user for some internal step.
     !
     ! CV_ERR_FAILURE:  Error test failures occurred too many times
     !                  (= MXNEF = 7) during one internal time step or
     !                  occurred with |h| = hmin.
     !
     ! CV_CONV_FAILURE: Convergence test failures occurred too many
     !                  times (= MXNCF = 10) during one internal time
     !                  step or occurred with |h| = hmin.
     !
     ! CV_LINIT_FAIL:   The linear solver's initialization function
     !                  failed.
     !
     ! CV_LSETUP_FAIL:  The linear solver's setup routine failed in an
     !                  unrecoverable manner.
     !
     ! CV_LSOLVE_FAIL:  The linear solver's solve routine failed in an
     !                  unrecoverable manner.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVode(cvode_mem, tout, yout, tret, itask) &
          bind(C,name='CVode')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: tout
       type(c_ptr),    value :: yout
       real(c_double)        :: tret
       integer(c_int), value :: itask
     end function FCVode

     ! -----------------------------------------------------------------
     ! Function : FCVodeGetDky
     ! -----------------------------------------------------------------
     ! FCVodeGetDky computes the kth derivative of the y function at
     ! time t, where tn-hu <= t <= tn, tn denotes the current
     ! internal time reached, and hu is the last internal step size
     ! successfully used by the solver. The user may request
     ! k=0, 1, ..., qu, where qu is the order last used. The
     ! derivative vector is returned in dky. This vector must be
     ! allocated by the caller. It is only legal to call this
     ! function after a successful return from CVode.
     !
     ! cvode_mem is the pointer to CVODE memory returned by
     !           FCVodeCreate.
     !
     ! t   is the time at which the kth derivative of y is evaluated.
     !     The legal range for t is [tn-hu,tn] as described above.
     !
     ! k   is the order of the derivative of y to be computed. The
     !     legal range for k is [0,qu] as described above.
     !
     ! dky is the output derivative vector [((d/dy)^k)y](t).
     !
     ! The return value for FCVodeGetDky is one of:
     !
     !   CV_SUCCESS:  FCVodeGetDky succeeded.
     !
     !   CV_BAD_K:    k is not in the range 0, 1, ..., qu.
     !
     !   CV_BAD_T:    t is not in the interval [tn-hu,tn].
     !
     !   CV_BAD_DKY:  The dky argument was NULL.
     !
     !   CV_MEM_NULL: The cvode_mem argument was NULL.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeGetDky(cvode_mem, t, k, dky) &
          bind(C,name='CVodeGetDky')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: t
       integer(c_int), value :: k
       type(c_ptr),    value :: dky
     end function FCVodeGetDky

     ! -----------------------------------------------------------------
     ! Integrator optional output extraction functions
     ! -----------------------------------------------------------------
     ! The following functions can be called to get optional outputs
     ! and statistics related to the main integrator.
     ! -----------------------------------------------------------------
     ! FCVodeGetWorkSpace returns the CVODE real and integer workspaces
     ! FCVodeGetNumSteps returns the cumulative number of internal
     !                  steps taken by the solver
     ! FCVodeGetNumRhsEvals returns the number of calls to the user's
     !                     f function
     ! FCVodeGetNumLinSolvSetups returns the number of calls made to
     !                          the linear solver's setup routine
     ! FCVodeGetNumErrTestFails returns the number of local error test
     !                         failures that have occured
     ! FCVodeGetLastOrder returns the order used during the last
     !                   internal step
     ! FCVodeGetCurrentOrder returns the order to be used on the next
     !                      internal step
     ! FCVodeGetNumStabLimOrderReds returns the number of order
     !                             reductions due to stability limit
     !                             detection
     ! FCVodeGetActualInitStep returns the actual initial step size
     !                        used by CVODE
     ! FCVodeGetLastStep returns the step size for the last internal
     !                  step
     ! FCVodeGetCurrentStep returns the step size to be attempted on
     !                     the next internal step
     ! FCVodeGetCurrentTime returns the current internal time reached
     !                     by the solver
     ! FCVodeGetTolScaleFactor returns a suggested factor by which the
     !                        user's tolerances should be scaled when
     !                        too much accuracy has been requested for
     !                        some internal step
     ! FCVodeGetErrWeights returns the current error weight vector.
     !                    The user must allocate space for eweight.
     ! FCVodeGetEstLocalErrors returns the vector of estimated local
     !                        errors. The user must allocate space
     !                        for ele.
     ! FCVodeGetNumGEvals returns the number of calls to the user's
     !                   g function (for rootfinding)
     ! FCVodeGetRootInfo returns the indices for which g_i was found to
     !                  have a root. The user must allocate space for
     !                  rootsfound. For i = 0 ... nrtfn-1,
     !                  rootsfound[i] = 1 if g_i has a root, and = 0 if not.
     !
     ! FCVodeGet* return values:
     !   CV_SUCCESS   if succesful
     !   CV_MEM_NULL  if the cvode memory was NULL
     !   CV_NO_SLDET  if stability limit was not turned on
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeGetWorkSpace(cvode_mem, lenrw, leniw) &
          bind(C,name='CVodeGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: lenrw
       integer(c_long)    :: leniw
     end function FCVodeGetWorkSpace

     integer(c_int) function FCVodeGetNumSteps(cvode_mem, nsteps) &
          bind(C,name='CVodeGetNumSteps')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nsteps
     end function FCVodeGetNumSteps

     integer(c_int) function FCVodeGetNumRhsEvals(cvode_mem, nfevals) &
          bind(C,name='CVodeGetNumRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nfevals
     end function FCVodeGetNumRhsEvals

     integer(c_int) function FCVodeGetNumLinSolvSetups(cvode_mem, nlinsetups) &
          bind(C,name='CVodeGetNumLinSolvSetups')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nlinsetups
     end function FCVodeGetNumLinSolvSetups

     integer(c_int) function FCVodeGetNumErrTestFails(cvode_mem, netfails) &
          bind(C,name='CVodeGetNumErrTestFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: netfails
     end function FCVodeGetNumErrTestFails

     integer(c_int) function FCVodeGetLastOrder(cvode_mem, qlast) &
          bind(C,name='CVodeGetLastOrder')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_int)    :: qlast
     end function FCVodeGetLastOrder

     integer(c_int) function FCVodeGetCurrentOrder(cvode_mem, qcur) &
          bind(C,name='CVodeGetCurrentOrder')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_int)     :: qcur
     end function FCVodeGetCurrentOrder

     integer(c_int) function FCVodeGetNumStabLimOrderReds(cvode_mem, nslred) &
          bind(C,name='CVodeGetNumStabLimOrderReds')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nslred
     end function FCVodeGetNumStabLimOrderReds

     integer(c_int) function FCVodeGetActualInitStep(cvode_mem, hinused) &
          bind(C,name='CVodeGetActualInitStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       real(c_double)     :: hinused
     end function FCVodeGetActualInitStep

     integer(c_int) function FCVodeGetLastStep(cvode_mem, hlast) &
          bind(C,name='CVodeGetLastStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       real(c_double)     :: hlast
     end function FCVodeGetLastStep

     integer(c_int) function FCVodeGetCurrentStep(cvode_mem, hcur) &
         bind(C,name='CVodeGetCurrentStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       real(c_double)     :: hcur
     end function FCVodeGetCurrentStep

     integer(c_int) function FCVodeGetCurrentTime(cvode_mem, tcur) &
          bind(C,name='CVodeGetCurrentTime')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       real(c_double)     :: tcur
     end function FCVodeGetCurrentTime

     integer(c_int) function FCVodeGetTolScaleFactor(cvode_mem, tolsfac) &
          bind(C,name='CVodeGetTolScaleFactor')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       real(c_double)     :: tolsfac
     end function FCVodeGetTolScaleFactor

     integer(c_int) function FCVodeGetErrWeights(cvode_mem, eweight) &
          bind(C,name='CVodeGetEffWeights')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       type(c_ptr), value :: eweight
     end function FCVodeGetErrWeights

     integer(c_int) function FCVodeGetEstLocalErrors(cvode_mem, ele) &
          bind(C,name='CVodeGetEstLocalErrors')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       type(c_ptr), value :: ele
     end function FCVodeGetEstLocalErrors

     integer(c_int) function FCVodeGetNumGEvals(cvode_mem, ngevals) &
          bind(C,name='CVodeGetNumGEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: ngevals
     end function FCVodeGetNumGEvals

     integer(c_int) function FCVodeGetRootInfo(cvode_mem, rootsfound) &
          bind(C,name='CVodeGetRootInfo')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_int)     :: rootsfound
     end function FCVodeGetRootInfo

     ! -----------------------------------------------------------------
     ! As a convenience, the following functions provides the
     ! optional outputs in one group.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeGetIntegratorStats(cvode_mem, nsteps, nfevals, &
          nlinsetups, netfails, qlast, qcur, hinused, hlast, hcur, tcur) &
          bind(C,name='CVodeGetIntegratorStats')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nsteps
       integer(c_long)    :: nfevals
       integer(c_long)    :: nlinsetups
       integer(c_long)    :: netfails
       integer(c_int)     :: qlast
       integer(c_int)     :: qcur
       real(c_double)     :: hinused
       real(c_double)     :: hlast
       real(c_double)     :: hcur
       real(c_double)     :: tcur
     end function FCVodeGetIntegratorStats

     ! -----------------------------------------------------------------
     ! Nonlinear solver optional output extraction functions
     ! -----------------------------------------------------------------
     ! The following functions can be called to get optional outputs
     ! and statistics related to the nonlinear solver.
     ! -----------------------------------------------------------------
     ! FCVodeGetNumNonlinSolvIters returns the number of nonlinear
     !                            solver iterations performed.
     ! FCVodeGetNumNonlinSolvConvFails returns the number of nonlinear
     !                                convergence failures.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeGetNumNonlinSolvIters(cvode_mem, nniters) &
          bind(C,name='CVodeGetNumNonlinSolvIters')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nniters
     end function FCVodeGetNumNonlinSolvIters

     integer(c_int) function FCVodeGetNumNonlinSolvConvFails(cvode_mem, nncfails) &
          bind(C,name='CVodeGetNumNonlinSolvConvFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nncfails
     end function FCVodeGetNumNonlinSolvConvFails

     ! -----------------------------------------------------------------
     ! As a convenience, the following function provides the
     ! nonlinear solver optional outputs in a group.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeGetNonlinSolvStats(cvode_mem, nniters, nncfails) &
          bind(C,name='CVodeGetNonlinSolvStats')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nniters
       integer(c_long)    :: nncfails
     end function FCVodeGetNonlinSolvStats

     ! -----------------------------------------------------------------
     ! The following function returns the name of the constant
     ! associated with a CVODE return flag
     ! -----------------------------------------------------------------

     ! >>> NOT CURRENTLY IMPLEMENTED IN FORTRAN INTERFACE
     ! char* CVodeGetReturnFlagName(long int flag);

     ! -----------------------------------------------------------------
     ! Subroutine : FCVodeFree
     ! -----------------------------------------------------------------
     ! FCVodeFree frees the problem memory cvode_mem allocated by
     ! FCVodeCreate and FCVodeInit. Its only argument is the pointer
     ! cvode_mem returned by FCVodeCreate.
     ! -----------------------------------------------------------------

     subroutine FCVodeFree(cvode_mem) &
          bind(C,name='CVodeFree')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: cvode_mem ! DO NOT use value attribute input is void**
     end subroutine FCVodeFree

     ! =================================================================
     ! Interfaces from cvode_band.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function : FCVBand
     ! -----------------------------------------------------------------
     ! A call to the FCVBand function links the main CVODE integrator
     ! with the CVBAND linear solver.
     !
     ! cvode_mem is the pointer to the integrator memory returned by
     !           FCVodeCreate.
     !
     ! N is the size of the ODE system.
     !
     ! mupper is the upper bandwidth of the band Jacobian
     !        approximation.
     !
     ! mlower is the lower bandwidth of the band Jacobian
     !        approximation.
     !
     ! The return value of FCVBand is one of:
     !    CVDLS_SUCCESS   if successful
     !    CVDLS_MEM_NULL  if the cvode memory was NULL
     !    CVDLS_MEM_FAIL  if there was a memory allocation failure
     !    CVDLS_ILL_INPUT if a required vector operation is missing or
     !                       if a bandwidth has an illegal value.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBand(cvode_mem, N, mupper, mlower) &
          bind(C,name='CVBand')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: cvode_mem
       integer(c_long), value :: N
       integer(c_long), value :: mupper
       integer(c_long), value :: mlower
     end function FCVBand

     ! =================================================================
     ! Interfaces from cvode_bandpre.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function : FCVBandPrecInit
     ! -----------------------------------------------------------------
     ! FCVBandPrecInit allocates and initializes the BANDPRE preconditioner
     ! module. This functino must be called AFTER one of the SPILS linear
     ! solver modules has been attached to the CVODE integrator.
     !
     ! The parameters of FCVBandPrecInit are as follows:
     !
     ! cvode_mem is the pointer to CVODE memory returned by FCVodeCreate.
     !
     ! N is the problem size.
     !
     ! mu is the upper half bandwidth.
     !
     ! ml is the lower half bandwidth.
     !
     ! The return value of FCVBandPrecInit is one of:
     !   CVSPILS_SUCCESS if no errors occurred
     !   CVSPILS_MEM_NULL if the integrator memory is NULL
     !   CVSPILS_LMEM_NULL if the linear solver memory is NULL
     !   CVSPILS_ILL_INPUT if an input has an illegal value
     !   CVSPILS_MEM_FAIL if a memory allocation request failed
     !
     ! NOTE: The band preconditioner assumes a serial implementation
     !       of the NVECTOR package. Therefore, FCVBandPrecInit will
     !       first test for a compatible N_Vector internal
     !       representation by checking for required functions.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBandPrecInit(cvode_mem, N, mu, ml) &
          bind(C,name='CVBandPrecInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: cvode_mem
       integer(c_long), value :: N
       integer(c_long), value :: mu
       integer(c_long), value :: ml
     end function FCVBandPrecInit

     ! -----------------------------------------------------------------
     ! Optional output functions : FCVBandPrecGet*
     ! -----------------------------------------------------------------
     ! FCVBandPrecGetWorkSpace returns the real and integer work space used
     !                        by FCVBANDPRE.
     ! FCVBandPrecGetNumRhsEvals returns the number of calls made from
     !                          FCVBANDPRE to the user's right-hand side
     !                          routine f.
     !
     ! The return value of FCVBandPrecGet* is one of:
     !   CVSPILS_SUCCESS if no errors occurred
     !   CVSPILS_MEM_NULL if the integrator memory is NULL
     !   CVSPILS_LMEM_NULL if the linear solver memory is NULL
     !   CVSPILS_PMEM_NULL if the preconditioner memory is NULL
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBandPrecGetWorkSpace(cvode_mem, lenrwLS, leniwLS) &
          bind(C,name='CVBandPrecGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FCVBandPrecGetWorkSpace

     integer(c_int) function FCVBandPrecGetNumRhsEvals(cvode_mem, nfevalsBP) &
          bind(C,name='CVBandPrecGetNumRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nfevalsBP
     end function FCVBandPrecGetNumRhsEvals

     ! =================================================================
     ! Interfaces from cvode_bbdpre.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function : FCVBBDPrecInit
     ! -----------------------------------------------------------------
     ! FCVBBDPrecInit allocates and initializes the BBD preconditioner.
     !
     ! The parameters of FCVBBDPrecInit are as follows:
     !
     ! cvode_mem is the pointer to the integrator memory.
     !
     ! Nlocal is the length of the local block of the vectors y etc.
     !        on the current processor.
     !
     ! mudq, mldq are the upper and lower half-bandwidths to be used
     !            in the difference quotient computation of the local
     !            Jacobian block.
     !
     ! mukeep, mlkeep are the upper and lower half-bandwidths of the
     !                retained banded approximation to the local Jacobian
     !                block.
     !
     ! dqrely is an optional input. It is the relative increment
     !        in components of y used in the difference quotient
     !        approximations. To specify the default, pass 0.
     !        The default is dqrely = sqrt(unit roundoff).
     !
     ! gloc is the name of the user-supplied function g(t,y) that
     !      approximates f and whose local Jacobian blocks are
     !      to form the preconditioner.
     !
     ! cfn is the name of the user-defined function that performs
     !     necessary interprocess communication for the
     !     execution of gloc.
     !
     ! The return value of FCVBBDPrecInit is one of:
     !   CVSPILS_SUCCESS if no errors occurred
     !   CVSPILS_MEM_NULL if the integrator memory is NULL
     !   CVSPILS_LMEM_NULL if the linear solver memory is NULL
     !   CVSPILS_ILL_INPUT if an input has an illegal value
     !   CVSPILS_MEM_FAIL if a memory allocation request failed
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBBDPrecInit(cvode_mem, Nlocal, mudq, mldq, &
          mukeep, mlkeep, dqrely, gloc, cfn) &
          bind(C,name='CVBBDPrecInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: cvode_mem
       integer(c_long), value :: Nlocal
       integer(c_long), value :: mudq
       integer(c_long), value :: mldq
       integer(c_long), value :: mukeep
       integer(c_long), value :: mlkeep
       real(c_double),  value :: dqrely
       type(c_funptr),  value :: gloc
       type(c_funptr),  value :: cfn
     end function FCVBBDPrecInit

     ! -----------------------------------------------------------------
     ! Function : FCVBBDPrecReInit
     ! -----------------------------------------------------------------
     ! FCVBBDPrecReInit re-initializes the BBDPRE module when solving a
     ! sequence of problems of the same size with CVSPGMR/CVBBDPRE or
     ! CVSPBCG/CVBBDPRE or CVSPTFQMR/CVBBDPRE provided there is no change
     ! in Nlocal, mukeep, or mlkeep. After solving one problem, and after
     ! calling FCVodeReInit to re-initialize the integrator for a subsequent
     ! problem, call FCVBBDPrecReInit.
     !
     ! All arguments have the same names and meanings as those
     ! of FCVBBDPrecInit.
     !
     ! The return value of FCVBBDPrecReInit is one of:
     !   CVSPILS_SUCCESS if no errors occurred
     !   CVSPILS_MEM_NULL if the integrator memory is NULL
     !   CVSPILS_LMEM_NULL if the linear solver memory is NULL
     !   CVSPILS_PMEM_NULL if the preconditioner memory is NULL
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBBDPrecReInit(cvode_mem, mudq, mldq, dqrely) &
          bind(C,name='CVBBNPrecReInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: cvode_mem
       integer(c_long), value :: mudq
       integer(c_long), value :: mldq
       real(c_double),  value :: dqrely
     end function FCVBBDPrecReInit

     ! -----------------------------------------------------------------
     ! BBDPRE optional output extraction routines
     ! -----------------------------------------------------------------
     ! FCVBBDPrecGetWorkSpace returns the BBDPRE real and integer work space
     !                       sizes.
     ! FCVBBDPrecGetNumGfnEvals returns the number of calls to gfn.
     !
     ! The return value of FCVBBDPrecGet* is one of:
     !   CVSPILS_SUCCESS if no errors occurred
     !   CVSPILS_MEM_NULL if the integrator memory is NULL
     !   CVSPILS_LMEM_NULL if the linear solver memory is NULL
     !   CVSPILS_PMEM_NULL if the preconditioner memory is NULL
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBBDPrecGetWorkSpace(cvode_mem, lenrwLS, leniwLS) &
          bind(C,name='CVBBDPrecGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FCVBBDPrecGetWorkSpace

     integer(c_int) function FCVBBDPrecGetNumGfnEvals(cvode_mem, ngevalsBBDP) &
          bind(C,name='CVBBDPrecGetNumGfnEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: ngevalsBBDP
     end function FCVBBDPrecGetNumGfnEvals

     ! =================================================================
     ! Interfaces from cvode_dense.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function: FCVDense
     ! -----------------------------------------------------------------
     ! A call to the FCVDense function links the main integrator with
     ! the CVDENSE linear solver.
     !
     ! cvode_mem is the pointer to the integrator memory returned by
     !           FCVodeCreate.
     !
     ! N is the size of the ODE system.
     !
     ! The return value of FCVDense is one of:
     !    CVDLS_SUCCESS   if successful
     !    CVDLS_MEM_NULL  if the cvode memory was NULL
     !    CVDLS_MEM_FAIL  if there was a memory allocation failure
     !    CVDLS_ILL_INPUT if a required vector operation is missing
     ! -----------------------------------------------------------------

     integer(c_int) function FCVDense(cvode_mem, N) &
          bind(C,name='CVDense')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: cvode_mem
       integer(c_long), value :: N
     end function FCVDense

     ! =================================================================
     ! Interfaces from cvode_diag.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function : FCVDiag
     ! -----------------------------------------------------------------
     ! A call to the FCVDiag function links the main integrator with
     ! the CVDIAG linear solver.
     !
     ! cvode_mem is the pointer to the integrator memory returned by
     !           FCVodeCreate.
     !
     ! The return value of FCVDiag is one of:
     !    CVDIAG_SUCCESS   if successful
     !    CVDIAG_MEM_NULL  if the cvode memory was NULL
     !    CVDIAG_MEM_FAIL  if there was a memory allocation failure
     !    CVDIAG_ILL_INPUT if a required vector operation is missing
     ! -----------------------------------------------------------------

     integer(c_int) function FCVDiag(cvode_mem) &
          bind(C,name='CVDiag')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
     end function FCVDiag

     ! -----------------------------------------------------------------
     ! Optional outputs from the CVDIAG linear solver
     ! -----------------------------------------------------------------
     !
     ! FCVDiagGetWorkSpace returns the real and integer workspace used
     !                    by CVDIAG.
     ! FCVDiagGetNumRhsEvals returns the number of calls to the user
     !                      f routine due to finite difference Jacobian
     !                      evaluation.
     !                      Note: The number of diagonal approximate
     !                      Jacobians formed is equal to the number of
     !                      CVDiagSetup calls. This number is available
     !                      through FCVodeGetNumLinSolvSetups.
     ! FCVDiagGetLastFlag returns the last error flag set by any of
     !                   the CVDIAG interface functions.
     !
     ! The return value of FCVDiagGet* is one of:
     !    CVDIAG_SUCCESS   if successful
     !    CVDIAG_MEM_NULL  if the cvode memory was NULL
     !    CVDIAG_LMEM_NULL if the cvdiag memory was NULL
     ! -----------------------------------------------------------------

     integer(c_int) function FCVDiagGetWorkSpace(cvode_mem, lenrwLS, leniwLS) &
          bind(C,name='CVDiagGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FCVDiagGetWorkSpace

     integer(c_int) function FCVDiagGetNumRhsEvals(cvode_mem, nfevalsLS) &
          bind(C,name='CVDiagGetNumRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nfevalsLS
     end function FCVDiagGetNumRhsEvals

     integer(c_int) function FCVDiagGetLastFlag(cvode_mem, flag) &
          bind(C,name='CVDiagGetLastFlag')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: flag
     end function FCVDiagGetLastFlag

     ! -----------------------------------------------------------------
     ! The following function returns the name of the constant
     ! associated with a CVDIAG return flag
     ! -----------------------------------------------------------------

     ! >>> NOT CURRENTLY IMPLEMENTED IN FORTRAN INTERFACE
     ! char* CVDiagGetReturnFlagName(long int flag);

     ! =================================================================
     ! Interfaces from cvode_direct.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Optional inputs to the CVDLS linear solver
     ! -----------------------------------------------------------------
     !
     ! FCVDlsSetDenseJacFn specifies the dense Jacobian approximation
     ! routine to be used for a direct dense linear solver.
     !
     ! FCVDlsSetBandJacFn specifies the band Jacobian approximation
     ! routine to be used for a direct band linear solver.
     !
     ! By default, a difference quotient approximation, supplied with
     ! the solver is used.
     !
     ! The return value is one of:
     !    CVDLS_SUCCESS   if successful
     !    CVDLS_MEM_NULL  if the CVODE memory was NULL
     !    CVDLS_LMEM_NULL if the linear solver memory was NULL
     ! -----------------------------------------------------------------

     integer(c_int) function FCVDlsSetDenseJacFn(cvode_mem, jac) &
          bind(C,name='CVDlsSetDenseJacFn')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: jac
     end function FCVDlsSetDenseJacFn

     integer(c_int) function FCVDlsSetBandJacFn(cvode_mem, jac) &
          bind(C,name='CVDlsSetBandJacFn')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: jac
     end function FCVDlsSetBandJacFn

     ! -----------------------------------------------------------------
     ! Optional outputs from the CVDLS linear solver
     ! -----------------------------------------------------------------
     !
     ! FCVDlsGetWorkSpace   returns the real and integer workspace used
     !                     by the direct linear solver.
     ! FCVDlsGetNumJacEvals returns the number of calls made to the
     !                     Jacobian evaluation routine jac.
     ! FCVDlsGetNumRhsEvals returns the number of calls to the user
     !                     f routine due to finite difference Jacobian
     !                     evaluation.
     ! FCVDlsGetLastFlag    returns the last error flag set by any of
     !                     the CVDLS interface functions.
     !
     ! The return value of FCVDlsGet* is one of:
     !    CVDLS_SUCCESS   if successful
     !    CVDLS_MEM_NULL  if the CVODE memory was NULL
     !    CVDLS_LMEM_NULL if the linear solver memory was NULL
     ! -----------------------------------------------------------------

     integer(c_int) function FCVDlsGetWorkSpace(cvode_mem, lenrwLS, leniwLS) &
          bind(C,name='CVDlsGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FCVDlsGetWorkSpace

     integer(c_int) function FCVDlsGetNumJacEvals(cvode_mem, njevals) &
          bind(C,name='CVDlsGetNumJacEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: njevals
     end function FCVDlsGetNumJacEvals

     integer(c_int) function FCVDlsGetNumRhsEvals(cvode_mem, nfevalsLS) &
          bind(C,name='CVDlsGetNumRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nfevalsLS
     end function FCVDlsGetNumRhsEvals

     integer(c_int) function FCVDlsGetLastFlag(cvode_mem, flag) &
          bind(C,name='CVDlsGetLastFlag')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: flag
     end function FCVDlsGetLastFlag

     ! -----------------------------------------------------------------
     ! The following function returns the name of the constant
     ! associated with a CVDLS return flag
     ! -----------------------------------------------------------------

     ! >>> NOT CURRENTLY IMPLEMENTED IN FORTRAN INTERFACE
     ! char* CVDlsGetReturnFlagName(long int flag);

     ! =================================================================
     ! Interfaces from cvode_hypamgpre.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function : FCVBoomerAMGInit
     ! -----------------------------------------------------------------
     ! FCVBoomerAMGInit allocates and initializes the BBD preconditioner.
     !
     ! The parameters of FCVBoomerAMGInit are as follows:
     !
     ! cvode_mem is the pointer to the integrator memory.
     !
     ! Nlocal is the length of the local block of the vectors y etc.
     !        on the current processor.
     !
     ! mudq, mldq are the upper and lower half-bandwidths to be used
     !            in the difference quotient computation of the local
     !            Jacobian block.
     !
     ! mukeep, mlkeep are the upper and lower half-bandwidths of the
     !                retained banded approximation to the local Jacobian
     !                block.
     !
     ! dqrely is an optional input. It is the relative increment
     !        in components of y used in the difference quotient
     !        approximations. To specify the default, pass 0.
     !        The default is dqrely = sqrt(unit roundoff).
     !
     ! gloc is the name of the user-supplied function g(t,y) that
     !      approximates f and whose local Jacobian blocks are
     !      to form the preconditioner.
     !
     ! cfn is the name of the user-defined function that performs
     !     necessary interprocess communication for the
     !     execution of gloc.
     !
     ! The return value of FCVBoomerAMGInit is one of:
     !   CVSPILS_SUCCESS if no errors occurred
     !   CVSPILS_MEM_NULL if the integrator memory is NULL
     !   CVSPILS_LMEM_NULL if the linear solver memory is NULL
     !   CVSPILS_ILL_INPUT if an input has an illegal value
     !   CVSPILS_MEM_FAIL if a memory allocation request failed
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBoomerAMGInit(cvode_mem, ilower, iupper, &
          jlower, jupper, N) &
          bind(C,name='CVBoomerAMGInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: ilower
       integer(c_int), value :: iupper
       integer(c_int), value :: jlower
       integer(c_int), value :: jupper
       integer(c_int), value :: N
     end function FCVBoomerAMGInit

     ! -----------------------------------------------------------------
     ! Function : FCVBoomerAMGReInit
     ! -----------------------------------------------------------------
     ! FCVBoomerAMGReInit re-initializes the HYPRE_BOOMERAMG module when solving a
     ! sequence of problems of the same size with CVSPGMR/CVHYPRE_BOOMERAMG or
     ! CVSPBCG/CVHYPRE_BOOMERAMG or CVSPTFQMR/CVHYPRE_BOOMERAMG provided there is no change
     ! in Nlocal, mukeep, or mlkeep. After solving one problem, and after
     ! calling FCVodeReInit to re-initialize the integrator for a subsequent
     ! problem, call FCVBoomerAMGReInit.
     !
     ! All arguments have the same names and meanings as those
     ! of FCVBoomerAMGInit.
     !
     ! The return value of FCVBoomerAMGReInit is one of:
     !   CVSPILS_SUCCESS if no errors occurred
     !   CVSPILS_MEM_NULL if the integrator memory is NULL
     !   CVSPILS_LMEM_NULL if the linear solver memory is NULL
     !   CVSPILS_PMEM_NULL if the preconditioner memory is NULL
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBoomerAMGReInit(cvode_mem, mudq, mldq, dqrely) &
          bind(C,name='CVBoomerAMGReInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: cvode_mem
       integer(c_long), value :: mudq
       integer(c_long), value :: mldq
       real(c_double),  value :: dqrely
     end function FCVBoomerAMGReInit

     ! -----------------------------------------------------------------
     ! HYPRE_BOOMERAMG optional output extraction routines
     ! -----------------------------------------------------------------
     ! FCVBoomerAMGGetWorkSpace returns the HYPRE_BOOMERAMG real and integer work space
     !                       sizes.
     ! FCVBoomerAMGGetNumGfnEvals returns the number of calls to gfn.
     !
     ! The return value of FCVBoomerAMGGet* is one of:
     !   CVSPILS_SUCCESS if no errors occurred
     !   CVSPILS_MEM_NULL if the integrator memory is NULL
     !   CVSPILS_LMEM_NULL if the linear solver memory is NULL
     !   CVSPILS_PMEM_NULL if the preconditioner memory is NULL
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBoomerAMGGetWorkSpace(cvode_mem, lenrwLS, leniwLS) &
          bind(C,name='CVBoomerAMGGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FCVBoomerAMGGetWorkSpace

     integer(c_int) function FCVBoomerAMGGetNumGfnEvals(cvode_mem, ngevalsBBDP) &
          bind(C,name='CVBoomerAMGGetNumGfnEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: ngevalsBBDP
     end function FCVBoomerAMGGetNumGfnEvals

     ! =================================================================
     ! Interfaces from cvode_klu.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function : FCVKLU
     ! -----------------------------------------------------------------
     ! A call to the FCVKLU function links the main integrator
     ! with the CVKLU linear solver module.
     !
     ! cv_mem is the pointer to integrator memory returned by
     !     FCVCreate.
     !
     !
     ! FCVKLU returns:
     !     CVSLU_SUCCESS   = 0  if successful
     !     CVSLU_LMEM_FAIL = -1 if there was a memory allocation failure
     !     CVSLU_ILL_INPUT = -2 if NVECTOR found incompatible
     !
     ! NOTE: The KLU linear solver assumes a serial implementation
     !       of the NVECTOR package. Therefore, CVKLU will first
     !       test for a compatible N_Vector internal representation
     !       by checking that the functions N_VGetArrayPointer and
     !       N_VSetArrayPointer exist.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVKLU(cvode_mem, n, nnz, sparsetype) &
          bind(C,name='CVKLU')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: n
       integer(c_int), value :: nnz
       integer(c_int), value :: sparsetype
     end function FCVKLU

     ! -----------------------------------------------------------------
     ! FCVKLUReInit
     ! -----------------------------------------------------------------
     ! This routine reinitializes memory and flags for a new factorization
     ! (symbolic and numeric) to be conducted at the next solver setup
     ! call.  This routine is useful in the cases where the number of nonzeroes
     ! has changed or if the structure of the linear system has changed
     ! which would require a new symbolic (and numeric factorization).
     !
     ! The reinit_type argumenmt governs the level of reinitialization:
     !
     ! reinit_type = 1: The Jacobian matrix will be destroyed and
     !                  a new one will be allocated based on the nnz
     !                  value passed to this call. New symbolic and
     !                  numeric factorizations will be completed at the next
     !                  solver setup.
     !
     ! reinit_type = 2: Only symbolic and numeric factorizations will be
     !                  completed.  It is assumed that the Jacobian size
     !                  has not exceeded the size of nnz given in the prior
     !                  call to FCVKLU.
     !
     ! This routine assumes no other changes to solver use are necessary.
     !
     ! The return value is CVSLS_SUCCESS = 0, CVSLS_MEM_NULL = -1,
     ! CVSLS_LMEM_NULL = -2, CVSLS_ILL_INPUT = -3, or CVSLS_MEM_FAIL = -4.
     !
     ! -----------------------------------------------------------------

     integer(c_int) function FCVKLUReInit(cvode_mem, n, nnz, reinit_type) &
          bind(C,name='CVKLUReInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: n
       integer(c_int), value :: nnz
       integer(c_int), value :: reinit_type
     end function FCVKLUReInit

     ! -----------------------------------------------------------------
     ! Optional Input Specification Functions
     ! -----------------------------------------------------------------
     !
     ! FCVKLUSetOrdering sets the ordering used by KLU for reducing fill.
     ! Options are: 0 for AMD, 1 for COLAMD, and 2 for the natural ordering.
     ! The default used in CVODE is 1 for COLAMD.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVKLUSetOrdering(cvode_mem, ordering_choice) &
          bind(C,name='CVKLUSetOrdering')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: ordering_choice
     end function FCVKLUSetOrdering

     ! =================================================================
     ! Interfaces from cvode_lapack.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function : FCVLapackDense
     ! -----------------------------------------------------------------
     ! A call to the FCVLapackDense function links the main integrator
     ! with the CVLAPACK linear solver using dense Jacobians.
     !
     ! cvode_mem is the pointer to the integrator memory returned by
     !           FCVodeCreate.
     !
     ! N is the size of the ODE system.
     !
     ! The return value of FCVLapackDense is one of:
     !    CVLAPACK_SUCCESS   if successful
     !    CVLAPACK_MEM_NULL  if the CVODE memory was NULL
     !    CVLAPACK_MEM_FAIL  if there was a memory allocation failure
     !    CVLAPACK_ILL_INPUT if a required vector operation is missing
     ! -----------------------------------------------------------------

     integer(c_int) function FCVLapackDense(cvode_mem, N) &
          bind(C,name='CVLapackDense')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: N
     end function FCVLapackDense

     ! -----------------------------------------------------------------
     ! Function : CVLapackBand
     ! -----------------------------------------------------------------
     ! A call to the FCVLapackBand function links the main integrator
     ! with the CVLAPACK linear solver using banded Jacobians.
     !
     ! cvode_mem is the pointer to the integrator memory returned by
     !           FCVodeCreate.
     !
     ! N is the size of the ODE system.
     !
     ! mupper is the upper bandwidth of the band Jacobian approximation.
     !
     ! mlower is the lower bandwidth of the band Jacobian approximation.
     !
     ! The return value of FCVLapackBand is one of:
     !    CVLAPACK_SUCCESS   if successful
     !    CVLAPACK_MEM_NULL  if the CVODE memory was NULL
     !    CVLAPACK_MEM_FAIL  if there was a memory allocation failure
     !    CVLAPACK_ILL_INPUT if a required vector operation is missing or
     !                       if a bandwidth has an illegal value.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVLapackBand(cvode_mem, N, mupper, mlower) &
          bind(C,name='CVLapackBand')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: N
       integer(c_int), value :: mupper
       integer(c_int), value :: mlower
     end function FCVLapackBand

     ! =================================================================
     ! Interfaces from cvode_sparse.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Optional inputs to the CVSPARSE linear solver
     ! -----------------------------------------------------------------
     ! FCVSlsSetSparseJacFn specifies the Jacobian approximation
     ! routine to be used for a sparse direct linear solver.
     !
     ! The return value is one of:
     !    CVSLS_SUCCESS   if successful
     !    CVSLS_MEM_NULL  if the CVODE memory was NULL
     !    CVSLS_LMEM_NULL if the linear solver memory was NULL
     ! -----------------------------------------------------------------

     integer(c_int) function FCVSlsSetSparseJacFn(cvode_mem, jac) &
          bind(C,name='CVSlsSetSparseJacFn')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: jac
     end function FCVSlsSetSparseJacFn

     ! -----------------------------------------------------------------
     ! Optional outputs from the CVSLS linear solver
     ! -----------------------------------------------------------------
     !
     ! FCVSlsGetNumJacEvals returns the number of calls made to the
     !                      Jacobian evaluation routine jac.
     ! FCVSlsGetLastFlag    returns the last error flag set by any of
     !                      the IDADLS interface functions.
     !
     ! The return value of IDADlsGet* is one of:
     !    CVSLS_SUCCESS   if successful
     !    CVSLS_MEM_NULL  if the IDA memory was NULL
     !    CVSLS_LMEM_NULL if the linear solver memory was NULL
     ! -----------------------------------------------------------------

     integer(c_int) function FCVSlsGetNumJacEvals(cvode_mem, njevals) &
          bind(C,name='CVSlsGetNumJacEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: njevals
     end function FCVSlsGetNumJacEvals

     integer(c_int) function FCVSlsGetLastFlag(cvode_mem, flag) &
          bind(C,name='CVSlsGetLastFlag')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: flag
     end function FCVSlsGetLastFlag

     ! -----------------------------------------------------------------
     ! The following function returns the name of the constant
     ! associated with a CVSLS return flag
     ! -----------------------------------------------------------------

     ! >>> NOT CURRENTLY IMPLEMENTED IN FORTRAN INTERFACE
     ! char* CVSlsGetReturnFlagName(long int flag);

     ! =================================================================
     ! Interfaces from cvode_spbcgs.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function : FCVSpbcg
     ! -----------------------------------------------------------------
     ! A call to the FCVSpbcg function links the main CVODE integrator
     ! with the CVSPBCG linear solver.
     !
     ! cvode_mem is the pointer to the integrator memory returned by
     !           FCVodeCreate.
     !
     ! pretype   is the type of user preconditioning to be done.
     !           This must be one of the four enumeration constants
     !           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined
     !           in iterative.h. These correspond to no preconditioning,
     !           left preconditioning only, right preconditioning
     !           only, and both left and right preconditioning,
     !           respectively.
     !
     ! maxl      is the maximum Krylov dimension. This is an
     !           optional input to the CVSPBCG solver. Pass 0 to
     !           use the default value CVSPBCG_MAXL=5.
     !
     ! The return value of FCVSpbcg is one of:
     !    CVSPILS_SUCCESS   if successful
     !    CVSPILS_MEM_NULL  if the cvode memory was NULL
     !    CVSPILS_MEM_FAIL  if there was a memory allocation failure
     !    CVSPILS_ILL_INPUT if a required vector operation is missing
     ! The above constants are defined in cvode_spils.h
     !
     ! -----------------------------------------------------------------

     integer(c_int) function FCVSpbcg(cvode_mem, pretype, maxl) &
          bind(C,name='CVSpbcg')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: pretype
       integer(c_int), value :: maxl
     end function FCVSpbcg

     ! =================================================================
     ! Interfaces from cvode_spgmr.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function : FCVSpgmr
     ! -----------------------------------------------------------------
     ! A call to the FCVSpgmr function links the main CVODE integrator
     ! with the CVSPGMR linear solver.
     !
     ! cvode_mem is the pointer to the integrator memory returned by
     !           FCVodeCreate.
     !
     ! pretype   is the type of user preconditioning to be done.
     !           This must be one of the four enumeration constants
     !           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined
     !           in sundials_iterative.h.
     !           These correspond to no preconditioning,
     !           left preconditioning only, right preconditioning
     !           only, and both left and right preconditioning,
     !           respectively.
     !
     ! maxl      is the maximum Krylov dimension. This is an
     !           optional input to the CVSPGMR solver. Pass 0 to
     !           use the default value CVSPGMR_MAXL=5.
     !
     ! The return value of CVSpgmr is one of:
     !    CVSPILS_SUCCESS   if successful
     !    CVSPILS_MEM_NULL  if the cvode memory was NULL
     !    CVSPILS_MEM_FAIL  if there was a memory allocation failure
     !    CVSPILS_ILL_INPUT if a required vector operation is missing
     ! The above constants are defined in cvode_spils.h
     !
     ! -----------------------------------------------------------------

     integer(c_int) function FCVSpgmr(cvode_mem, pretype, maxl) &
          bind(C,name='CVSpgmr')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: pretype
       integer(c_int), value :: maxl
     end function FCVSpgmr

     ! =================================================================
     ! Interfaces from cvode_spils.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Optional inputs to the CVSPILS linear solver
     ! -----------------------------------------------------------------
     !
     ! FCVSpilsSetPrecType resets the type of preconditioner, pretype,
     !                from the value previously set.
     !                This must be one of PREC_NONE, PREC_LEFT,
     !                PREC_RIGHT, or PREC_BOTH.
     !
     ! FCVSpilsSetGSType specifies the type of Gram-Schmidt
     !                orthogonalization to be used. This must be one of
     !                the two enumeration constants MODIFIED_GS or
     !                CLASSICAL_GS defined in iterative.h. These correspond
     !                to using modified Gram-Schmidt and classical
     !                Gram-Schmidt, respectively.
     !                Default value is MODIFIED_GS.
     !
     ! FCVSpilsSetMaxl resets the maximum Krylov subspace size, maxl,
     !                from the value previously set.
     !                An input value <= 0, gives the default value.
     !
     ! FCVSpilsSetEpsLin specifies the factor by which the tolerance on
     !                the nonlinear iteration is multiplied to get a
     !                tolerance on the linear iteration.
     !                Default value is 0.05.
     !
     ! FCVSpilsSetPreconditioner specifies the PrecSetup and PrecSolve functions.
     !                Default is NULL for both arguments (no preconditioning)
     !
     ! FCVSpilsSetJacTimesVecFn specifies the jtimes function. Default is to
     !                use an internal finite difference approximation routine.
     !
     ! The return value of FCVSpilsSet* is one of:
     !    CVSPILS_SUCCESS   if successful
     !    CVSPILS_MEM_NULL  if the cvode memory was NULL
     !    CVSPILS_LMEM_NULL if the linear solver memory was NULL
     !    CVSPILS_ILL_INPUT if an input has an illegal value
     ! -----------------------------------------------------------------

     integer(c_int) function FCVSpilsSetPrecType(cvode_mem, pretype) &
          bind(C,name='CVSpilsSetPrecType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: pretype
     end function FCVSpilsSetPrecType

     integer(c_int) function FCVSpilsSetGSType(cvode_mem, gstype) &
          bind(C,name='CVspilsSetGSType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: gstype
     end function FCVSpilsSetGSType

     integer(c_int) function FCVSpilsSetMaxl(cvode_mem, maxl) &
          bind(C,name='CVSpilsSetMaxl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: maxl
     end function FCVSpilsSetMaxl

     integer(c_int) function FCVSpilsSetEpsLin(cvode_mem, eplifac) &
          bind(C,name='FCVSpilsSetEpsLin')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: eplifac
     end function FCVSpilsSetEpsLin

     integer(c_int) function FCVSpilsSetPreconditioner(cvode_mem, pset, psolve) &
          bind(C,name='CVSpilsSetPreconditioner')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: pset
       type(c_funptr), value :: psolve
     end function FCVSpilsSetPreconditioner

     integer(c_int) function FCVSpilsSetJacTimesVecFn(cvode_mem, jtv) &
          bind(C,name='CVSpilsSetJacTimesVecFn')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: jtv
     end function FCVSpilsSetJacTimesVecFn

     ! -----------------------------------------------------------------
     ! Optional outputs from the CVSPILS linear solver
     ! -----------------------------------------------------------------
     ! FCVSpilsGetWorkSpace returns the real and integer workspace used
     !                by the SPILS module.
     !
     ! FCVSpilsGetNumPrecEvals returns the number of preconditioner
     !                 evaluations, i.e. the number of calls made
     !                 to PrecSetup with jok==FALSE.
     !
     ! FCVSpilsGetNumPrecSolves returns the number of calls made to
     !                 PrecSolve.
     !
     ! FCVSpilsGetNumLinIters returns the number of linear iterations.
     !
     ! FCVSpilsGetNumConvFails returns the number of linear
     !                 convergence failures.
     !
     ! FCVSpilsGetNumJtimesEvals returns the number of calls to jtimes.
     !
     ! FCVSpilsGetNumRhsEvals returns the number of calls to the user
     !                 f routine due to finite difference Jacobian
     !                 times vector evaluation.
     !
     ! FCVSpilsGetLastFlag returns the last error flag set by any of
     !                 the CVSPILS interface functions.
     !
     ! The return value of FCVSpilsGet* is one of:
     !    CVSPILS_SUCCESS   if successful
     !    CVSPILS_MEM_NULL  if the cvode memory was NULL
     !    CVSPILS_LMEM_NULL if the linear solver memory was NULL
     ! -----------------------------------------------------------------

     integer(c_int) function FCVSpilsGetWorkSpace(cvode_mem, lenrwLS, leniwLS) &
          bind(C,name='CVSpilsGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FCVSpilsGetWorkSpace

     integer(c_int) function FCVSpilsGetNumPrecEvals(cvode_mem, npevals) &
          bind(C,name='CVSpilsGetNumPrecEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: npevals
     end function FCVSpilsGetNumPrecEvals

     integer(c_int) function FCVSpilsGetNumPrecSolves(cvode_mem, npsolves) &
          bind(C,name='CVSpilsGetNumPrecSolves')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: npsolves
     end function FCVSpilsGetNumPrecSolves

     integer(c_int) function FCVSpilsGetNumLinIters(cvode_mem, nliters) &
          bind(C,name='CVSpilsGetNumLinIters')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nliters
     end function FCVSpilsGetNumLinIters

     integer(c_int) function FCVSpilsGetNumConvFails(cvode_mem, nlcfails) &
          bind(C,name='CVSpilsGetNumConvFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nlcfails
     end function FCVSpilsGetNumConvFails

     integer(c_int) function FCVSpilsGetNumJtimesEvals(cvode_mem, njvevals) &
          bind(C,name='CVSpilsGetNumJtimesEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: njvevals
     end function FCVSpilsGetNumJtimesEvals

     integer(c_int) function FCVSpilsGetNumRhsEvals(cvode_mem, nfevalsLS) &
          bind(C,name='CVSpilsGetNumRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nfevalsLS
     end function FCVSpilsGetNumRhsEvals

     integer(c_int) function FCVSpilsGetLastFlag(cvode_mem, flag) &
          bind(C,name='CVSpilsGetLastFalg')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: flag
     end function FCVSpilsGetLastFlag

     ! -----------------------------------------------------------------
     ! The following function returns the name of the constant
     ! associated with a CVSPILS return flag
     ! -----------------------------------------------------------------

     ! >>> NOT CURRENTLY IMPLEMENTED IN FORTRAN INTERFACE
     ! char* CVSpilsGetReturnFlagName(long int flag);

     ! =================================================================
     ! Interfaces from cvode_sptfqmr.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function : FCVSptfqmr
     ! -----------------------------------------------------------------
     ! A call to the FCVSptfqmr function links the main CVODE integrator
     ! with the CVSPTFQMR linear solver.
     !
     ! cvode_mem is the pointer to the integrator memory returned by
     !           FCVodeCreate.
     !
     ! pretype   is the type of user preconditioning to be done.
     !           This must be one of the four enumeration constants
     !           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined
     !           in iterative.h. These correspond to no preconditioning,
     !           left preconditioning only, right preconditioning
     !           only, and both left and right preconditioning,
     !           respectively.
     !
     ! maxl      is the maximum Krylov dimension. This is an
     !           optional input to the CVSPTFQMR solver. Pass 0 to
     !           use the default value CVSPILS_MAXL=5.
     !
     ! The return value of FCVSptfqmr is one of:
     !    CVSPILS_SUCCESS   if successful
     !    CVSPILS_MEM_NULL  if the cvode memory was NULL
     !    CVSPILS_MEM_FAIL  if there was a memory allocation failure
     !    CVSPILS_ILL_INPUT if a required vector operation is missing
     ! The above constants are defined in cvode_spils.h
     !
     ! -----------------------------------------------------------------

     integer(c_int) function FCVSptfqmr(cvode_mem, pretype, maxl) &
          bind(C,name='CVSptfqmr')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: pretype
       integer(c_int), value :: maxl
     end function FCVSptfqmr

     ! =================================================================
     ! Interfaces from cvode_superlumt.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! Function : FCVSuperLUMT
     ! -----------------------------------------------------------------
     ! A call to the FCVSuperLUMT function links the main integrator
     ! with the CVSuperLUMT linear solver module.
     !
     ! cv_mem is the pointer to integrator memory returned by
     !     FCVCreate.
     !
     !
     ! FCVSuperLUMT returns:
     !     CVSLU_SUCCESS   = 0  if successful
     !     CVSLU_LMEM_FAIL = -1 if there was a memory allocation failure
     !     CVSLU_ILL_INPUT = -2 if NVECTOR found incompatible
     !
     ! NOTE: The CVSuperLUMT linear solver assumes a serial implementation
     !       of the NVECTOR package. Therefore, CVSuperLUMT will first
     !       test for a compatible N_Vector internal representation
     !       by checking that the functions N_VGetArrayPointer and
     !       N_VSetArrayPointer exist.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVSuperLUMT(cvode_mem, num_threads, n, nnz) &
          bind(C,name='CVSuperLUMT')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: num_threads
       integer(c_int), value :: n
       integer(c_int), value :: nnz
     end function FCVSuperLUMT

     ! -----------------------------------------------------------------
     ! Optional Input Specification Functions
     ! -----------------------------------------------------------------
     !
     ! FCVSuperLUMTSetOrdering sets the ordering used by CVSuperLUMT for
     ! reducing fill.
     ! Options are:
     ! 0 for natural ordering
     ! 1 for minimal degree ordering on A'*A
     ! 2 for minimal degree ordering on A'+A
     ! 3 for approximate minimal degree ordering for unsymmetric matrices
     ! The default used in SUNDIALS is 3 for COLAMD.
     ! -----------------------------------------------------------------

     integer(c_int) function FCVSuperLUMTSetOrdering(cvode_mem, ordering_choice) &
          bind(C,name='CVSuperLUMTSetOrdering')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: ordering_choice
     end function FCVSuperLUMTSetOrdering

  end interface

end module cvode_interface
