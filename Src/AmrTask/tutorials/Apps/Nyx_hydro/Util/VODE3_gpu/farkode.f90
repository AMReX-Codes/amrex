! ------------------------------------------------------------------
! Programmer(s): David J. Gardner @ LLNL
!                Daniel R. Reynolds @ SMU
!                modified Jean M. Sexton @ LBL
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
! This file contains a Fortran module for interfacing directly with
! ARKODE using the ISO_C_BINDING module.
! ------------------------------------------------------------------

module farkode_mod

  use, intrinsic :: iso_c_binding, only : c_int

  ! =================================================================
  !              A R K O D E     C O N S T A N T S
  ! =================================================================

  ! -----------------------------------------------------------------
  ! Enumerations for inputs to ARKodeCreate and ARKode.
  ! -----------------------------------------------------------------

  ! itask
  integer(c_int), parameter :: ARK_NORMAL   = 1
  integer(c_int), parameter :: ARK_ONE_STEP = 2

  ! -----------------------------------------------------------------
  ! ARKODE return flags
  ! -----------------------------------------------------------------

  integer(c_int), parameter :: ARK_SUCCESS           =   0
  integer(c_int), parameter :: ARK_TSTOP_RETURN      =   1
  integer(c_int), parameter :: ARK_ROOT_RETURN       =   2

  integer(c_int), parameter :: ARK_WARNING           =  99

  integer(c_int), parameter :: ARK_TOO_MUCH_WORK     =  -1
  integer(c_int), parameter :: ARK_TOO_MUCH_ACC      =  -2
  integer(c_int), parameter :: ARK_ERR_FAILURE       =  -3
  integer(c_int), parameter :: ARK_CONV_FAILURE      =  -4

  integer(c_int), parameter :: ARK_LINIT_FAIL        =  -5
  integer(c_int), parameter :: ARK_LSETUP_FAIL       =  -6
  integer(c_int), parameter :: ARK_LSOLVE_FAIL       =  -7
  integer(c_int), parameter :: ARK_RHSFUNC_FAIL      =  -8
  integer(c_int), parameter :: ARK_FIRST_RHSFUNC_ERR =  -9
  integer(c_int), parameter :: ARK_REPTD_RHSFUNC_ERR =  -10
  integer(c_int), parameter :: ARK_UNREC_RHSFUNC_ERR =  -11
  integer(c_int), parameter :: ARK_RTFUNC_FAIL       =  -12

  integer(c_int), parameter :: ARK_MEM_FAIL          =  -20
  integer(c_int), parameter :: ARK_MEM_NULL          =  -21
  integer(c_int), parameter :: ARK_ILL_INPUT         =  -22
  integer(c_int), parameter :: ARK_NO_MALLOC         =  -23
  integer(c_int), parameter :: ARK_BAD_K             =  -24
  integer(c_int), parameter :: ARK_BAD_T             =  -25
  integer(c_int), parameter :: ARK_BAD_DKY           =  -26
  integer(c_int), parameter :: ARK_TOO_CLOSE         =  -27

  ! -----------------------------------------------------------------
  ! ARKDLS return values
  ! -----------------------------------------------------------------

  integer(c_int), parameter :: ARKDLS_SUCCESS         =  0
  integer(c_int), parameter :: ARKDLS_MEM_NULL        = -1
  integer(c_int), parameter :: ARKDLS_LMEM_NULL       = -2
  integer(c_int), parameter :: ARKDLS_ILL_INPUT       = -3
  integer(c_int), parameter :: ARKDLS_MEM_FAIL        = -4

  ! Additional last_flag values
  integer(c_int), parameter :: ARKDLS_JACFUNC_UNREARKR = -5
  integer(c_int), parameter :: ARKDLS_JACFUNC_REARKR   = -6
  integer(c_int), parameter :: ARKDLS_SUNMAT_FAIL     = -7

  ! -----------------------------------------------------------------
  ! ARKSPILS return values
  ! -----------------------------------------------------------------

  integer(c_int), parameter :: ARKSPILS_SUCCESS    =  0
  integer(c_int), parameter :: ARKSPILS_MEM_NULL   = -1
  integer(c_int), parameter :: ARKSPILS_LMEM_NULL  = -2
  integer(c_int), parameter :: ARKSPILS_ILL_INPUT  = -3
  integer(c_int), parameter :: ARKSPILS_MEM_FAIL   = -4
  integer(c_int), parameter :: ARKSPILS_PMEM_NULL  = -5
  integer(c_int), parameter :: ARKSPILS_SUNLS_NULL = -6

  ! -----------------------------------------------------------------
  ! ARKDIAG return values
  ! -----------------------------------------------------------------

  integer(c_int), parameter :: ARKDIAG_SUCCESS         =  0
  integer(c_int), parameter :: ARKDIAG_MEM_NULL        = -1
  integer(c_int), parameter :: ARKDIAG_LMEM_NULL       = -2
  integer(c_int), parameter :: ARKDIAG_ILL_INPUT       = -3
  integer(c_int), parameter :: ARKDIAG_MEM_FAIL        = -4

  ! Additional last_flag values
  integer(c_int), parameter :: ARKDIAG_INV_FAIL        = -5
  integer(c_int), parameter :: ARKDIAG_RHSFUNC_UNREARKR = -6
  integer(c_int), parameter :: ARKDIAG_RHSFUNC_REARKR   = -7

  ! =================================================================
  !          U S E R - C A L L A B L E   R O U T I N E S
  ! =================================================================

  interface

     ! =================================================================
     ! Interfaces for arkode.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! ARKodeCreate
     ! -----------------------------------------------------------------

     type(c_ptr) function FARKodeCreate() &
          bind(C,name='ARKodeCreate')
       use, intrinsic :: iso_c_binding
       implicit none
     end function FARKodeCreate

     ! -----------------------------------------------------------------
     ! Integrator optional input specification functions
     ! -----------------------------------------------------------------

     integer(c_int) function FARKodeSetErrHandlerFn(arkode_mem, ehfun, eh_data) &
          bind(C,name='ARKodeSetErrHandlerFn')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       type(c_funptr), value :: ehfun
       type(c_ptr),    value :: eh_data
     end function FARKodeSetErrHandlerFn

     ! -----------------------------------------------------------------
     ! NOT INTERFACED: ARKodeSetErrFile
     ! -----------------------------------------------------------------

     integer(c_int) function FARKodeSetUserData(arkode_mem, user_data) &
          bind(C,name='ARKodeSetUserData')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       type(c_ptr), value :: user_data
     end function FARKodeSetUserData

     integer(c_int) function FARKodeSetMaxOrd(arkode_mem, maxord) &
          bind(C,name='ARKodeSetMaxOrd')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       integer(c_int), value :: maxord
     end function FARKodeSetMaxOrd

     integer(c_int) function FARKodeSetMaxNumSteps(arkode_mem, mxsteps) &
          bind(C,name='ARKodeSetMaxNumSteps')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: arkode_mem
       integer(c_long), value :: mxsteps
     end function FARKodeSetMaxNumSteps

     integer(c_int) function FARKodeSetMaxHnilWarns(arkode_mem, mxhnil) &
          bind(C,name='ARKodeSetMaxHnilWarns')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       integer(c_int), value :: mxhnil
     end function FARKodeSetMaxHnilWarns

     integer(c_int) function FARKodeSetStabLimDet(arkode_mem, stldet) &
          bind(C,name='ARKodeSetStabLimDet')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       integer(c_int), value :: stldet
     end function FARKodeSetStabLimDet

     integer(c_int) function FARKodeSetInitStep(arkode_mem, hin) &
          bind(C,name='ARKodeSetInitStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       real(c_double), value :: hin
     end function FARKodeSetInitStep

     integer(c_int) function FARKodeSetMinStep(arkode_mem, hmin) &
          bind(C,name='ARKodeSetMinStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       real(c_double), value :: hmin
     end function FARKodeSetMinStep

     integer(c_int) function FARKodeSetMaxStep(arkode_mem, hmax) &
          bind(C,name='ARKodeSetMaxStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       real(c_double), value :: hmax
     end function FARKodeSetMaxStep

     integer(c_int) function FARKodeSetAdaptivityMethod(arkode_mem, imethod, idefault, pq, adapt_params) &
	bind(C,name='ARKodeSetAdaptivityMethod')
       use, intrinsic :: iso_c_binding
       type(c_ptr),    value :: arkode_mem
       integer(c_int), value :: imethod
       integer(c_int), value :: idefault
       integer(c_int), value :: pq
       type(c_ptr),    value :: adapt_params
     end function FARKodeSetAdaptivityMethod

     integer(c_int) function FARKodeSetStopTime(arkode_mem, tstop) &
          bind(C,name='ARKodeSetStopTime')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       real(c_double), value :: tstop
     end function FARKodeSetStopTime

     integer(c_int) function FARKodeSetMaxErrTestFails(arkode_mem, maxnef) &
          bind(C,name='ARKodeSetMaxErrTestFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       integer(c_int), value :: maxnef
     end function FARKodeSetMaxErrTestFails

     integer(c_int) function FARKodeSetMaxNonlinIters(arkode_mem, maxcor) &
          bind(C,name='ARKodeSetMaxNonlinIters')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       integer(c_int), value :: maxcor
     end function FARKodeSetMaxNonlinIters

     integer(c_int) function FARKodeSetMaxConvFails(arkode_mem, maxncf) &
          bind(C,name='ARKodeSetMaxConvFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       integer(c_int), value :: maxncf
     end function FARKodeSetMaxConvFails

     integer(c_int) function FARKodeSetNonlinConvCoef(arkode_mem, nlscoef) &
          bind(C,name='ARKodeSetNonlinConvCoef')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       real(c_double), value :: nlscoef
     end function FARKodeSetNonlinConvCoef

     integer(c_int) function FARKodeSetIterType(arkode_mem, iter) &
          bind(C,name='ARKodeSetIterType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       real(c_double), value :: iter
     end function FARKodeSetIterType

     integer(c_int) function FARKodeSetRootDirection(arkode_mem, rootdir) &
          bind(C,name='ARKodeSetRootDirection')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       real(c_double)     :: rootdir
     end function FARKodeSetRootDirection

     integer(c_int) function FARKodeSetNoInactiveRootWarn(arkode_mem) &
       bind(C,name='ARKodeSetNoInactiveRootWarn')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
     end function FARKodeSetNoInactiveRootWarn

     ! -----------------------------------------------------------------
     ! ARKodeInit
     ! -----------------------------------------------------------------

     integer(c_int) function FARKodeInit(arkode_mem, fe, fi, t0, y0) &
          bind(C,name='ARKodeInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       type(c_funptr), value :: fe
       type(c_funptr), value :: fi
       real(c_double), value :: t0
       type(c_ptr),    value :: y0
     end function FARKodeInit

     ! -----------------------------------------------------------------
     ! ARKodeReInit
     ! -----------------------------------------------------------------

     integer(c_int) function FARKodeReInit(arkode_mem, fe, fi, t0, y0) &
          bind(C,name='ARKodeReInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       type(c_funptr), value :: fe
       type(c_funptr), value :: fi
       real(c_double), value :: t0
       type(c_ptr),    value :: y0
     end function FARKodeReInit


     ! -----------------------------------------------------------------
     ! ARKodeSStolerances
     ! ARKodeSVtolerances
     ! ARKodeWFtolerances
     ! -----------------------------------------------------------------

     integer(c_int) function FARKodeSStolerances(arkode_mem, reltol, abstol) &
          bind(C,name='ARKodeSStolerances')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       real(c_double), value :: reltol
       real(c_double), value :: abstol
     end function FARKodeSStolerances

     integer(c_int) function FARKodeSVtolerances(arkode_mem, reltol, abstol) &
          bind(C,name='ARKodeSVtolerances')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       real(c_double), value :: reltol
       type(c_ptr),    value :: abstol
     end function FARKodeSVtolerances

     integer(c_int) function FARKodeWFtolerances(arkode_mem, efun) &
          bind(C,name='ARKodeWFtolerances')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       type(c_funptr), value :: efun
     end function FARKodeWFtolerances

     ! -----------------------------------------------------------------
     ! ARKodeRootInit
     ! -----------------------------------------------------------------

     integer(c_int) function FARKodeRootInit(arkode_mem, nrtfn, g) &
          bind(C,name='ARKodeRootInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       integer(c_int), value :: nrtfn
       type(c_funptr), value :: g
     end function FARKodeRootInit

     ! -----------------------------------------------------------------
     ! ARKode
     ! -----------------------------------------------------------------

     integer(c_int) function FARKode(arkode_mem, tout, yout, tret, itask) &
          bind(C,name='ARKode')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       real(c_double), value :: tout
       type(c_ptr),    value :: yout
       real(c_double)        :: tret
       integer(c_int), value :: itask
     end function FARKode

     ! -----------------------------------------------------------------
     ! ARKodeGetDky
     ! -----------------------------------------------------------------

     integer(c_int) function FARKodeGetDky(arkode_mem, t, k, dky) &
          bind(C,name='ARKodeGetDky')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       real(c_double), value :: t
       integer(c_int), value :: k
       type(c_ptr),    value :: dky
     end function FARKodeGetDky

     ! -----------------------------------------------------------------
     ! Integrator optional output extraction functions
     ! -----------------------------------------------------------------

     integer(c_int) function FARKodeGetWorkSpace(arkode_mem, lenrw, leniw) &
          bind(C,name='ARKodeGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: lenrw
       integer(c_long)    :: leniw
     end function FARKodeGetWorkSpace

     integer(c_int) function FARKodeGetNumSteps(arkode_mem, nsteps) &
          bind(C,name='ARKodeGetNumSteps')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nsteps
     end function FARKodeGetNumSteps

     integer(c_int) function FARKodeGetNumStepAttempts(arkode_mem, nst_a) &
          bind(C,name='ARKodeGetNumStepAttempts')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nst_a
     end function FARKodeGetNumStepAttempts

     integer(c_int) function FARKodeGetNumRhsEvals(arkode_mem, nfe,nfi) &
          bind(C,name='ARKodeGetNumRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nfe
       integer(c_long)    :: nfi
     end function FARKodeGetNumRhsEvals

     integer(c_int) function FARKodeGetNumLinSolvSetups(arkode_mem, nlinsetups) &
          bind(C,name='ARKodeGetNumLinSolvSetups')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nlinsetups
     end function FARKodeGetNumLinSolvSetups

     integer(c_int) function FARKodeGetNumErrTestFails(arkode_mem, netfails) &
          bind(C,name='ARKodeGetNumErrTestFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: netfails
     end function FARKodeGetNumErrTestFails

     integer(c_int) function FARKodeGetNumStabLimOrderReds(arkode_mem, nslred) &
          bind(C,name='ARKodeGetNumStabLimOrderReds')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nslred
     end function FARKodeGetNumStabLimOrderReds

     integer(c_int) function FARKodeGetActualInitStep(arkode_mem, hinused) &
          bind(C,name='ARKodeGetActualInitStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       real(c_double)     :: hinused
     end function FARKodeGetActualInitStep

     integer(c_int) function FARKodeGetLastStep(arkode_mem, hlast) &
          bind(C,name='ARKodeGetLastStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       real(c_double)     :: hlast
     end function FARKodeGetLastStep

     integer(c_int) function FARKodeGetCurrentStep(arkode_mem, hcur) &
         bind(C,name='ARKodeGetCurrentStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       real(c_double)     :: hcur
     end function FARKodeGetCurrentStep

     integer(c_int) function FARKodeGetCurrentTime(arkode_mem, tcur) &
          bind(C,name='ARKodeGetCurrentTime')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       real(c_double)     :: tcur
     end function FARKodeGetCurrentTime

     integer(c_int) function FARKodeGetTolScaleFactor(arkode_mem, tolsfac) &
          bind(C,name='ARKodeGetTolScaleFactor')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       real(c_double)     :: tolsfac
     end function FARKodeGetTolScaleFactor

     integer(c_int) function FARKodeGetErrWeights(arkode_mem, eweight) &
          bind(C,name='ARKodeGetEffWeights')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       type(c_ptr), value :: eweight
     end function FARKodeGetErrWeights

     integer(c_int) function FARKodeGetEstLocalErrors(arkode_mem, ele) &
          bind(C,name='ARKodeGetEstLocalErrors')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       type(c_ptr), value :: ele
     end function FARKodeGetEstLocalErrors

     integer(c_int) function FARKodeGetNumGEvals(arkode_mem, ngevals) &
          bind(C,name='ARKodeGetNumGEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: ngevals
     end function FARKodeGetNumGEvals

     integer(c_int) function FARKodeGetRootInfo(arkode_mem, rootsfound) &
          bind(C,name='ARKodeGetRootInfo')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_int)     :: rootsfound
     end function FARKodeGetRootInfo

     integer(c_int) function FARKodeGetIntegratorStats(arkode_mem, nsteps, nfevals, &
          nlinsetups, netfails, qlast, qcur, hinused, hlast, hcur, tcur) &
          bind(C,name='ARKodeGetIntegratorStats')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
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
     end function FARKodeGetIntegratorStats

     ! -----------------------------------------------------------------
     ! Nonlinear solver optional output extraction functions
     ! -----------------------------------------------------------------

     integer(c_int) function FARKodeGetNumNonlinSolvIters(arkode_mem, nniters) &
          bind(C,name='ARKodeGetNumNonlinSolvIters')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nniters
     end function FARKodeGetNumNonlinSolvIters

     integer(c_int) function FARKodeGetNumNonlinSolvConvFails(arkode_mem, nncfails) &
          bind(C,name='ARKodeGetNumNonlinSolvConvFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nncfails
     end function FARKodeGetNumNonlinSolvConvFails

     integer(c_int) function FARKodeGetNonlinSolvStats(arkode_mem, nniters, nncfails) &
          bind(C,name='ARKodeGetNonlinSolvStats')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nniters
       integer(c_long)    :: nncfails
     end function FARKodeGetNonlinSolvStats

     ! -----------------------------------------------------------------
     ! NOT INTERFACED: ARKodeGetReturnFlagName
     ! -----------------------------------------------------------------

     ! -----------------------------------------------------------------
     ! ARKodeFree
     ! -----------------------------------------------------------------

     subroutine FARKodeFree(arkode_mem) &
          bind(C,name='ARKodeFree')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: arkode_mem ! DO NOT use value attribute input is void**
     end subroutine FARKodeFree

     ! =================================================================
     ! Interfaces from arkode_direct.h
     ! =================================================================

     integer(c_int) function FARKDlsSetLinearSolver(arkode_mem, LS, A) &
          bind(C,name='ARKDlsSetLinearSolver')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       type(c_ptr), value :: LS
       type(c_ptr), value :: A
     end function FARKDlsSetLinearSolver

     ! -----------------------------------------------------------------
     ! Optional inputs to the ARKDLS linear solver
     ! -----------------------------------------------------------------

     integer(c_int) function FARKDlsSetJacFn(arkode_mem, jac) &
          bind(C,name='ARKDlsSetJacFn')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       type(c_funptr), value :: jac
     end function FARKDlsSetJacFn

     ! -----------------------------------------------------------------
     ! Optional outputs from the ARKDLS linear solver
     ! -----------------------------------------------------------------

     integer(c_int) function FARKDlsGetWorkSpace(arkode_mem, lenrwLS, leniwLS) &
          bind(C,name='ARKDlsGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FARKDlsGetWorkSpace

     integer(c_int) function FARKDlsGetNumJacEvals(arkode_mem, njevals) &
          bind(C,name='ARKDlsGetNumJacEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: njevals
     end function FARKDlsGetNumJacEvals

     integer(c_int) function FARKDlsGetNumRhsEvals(arkode_mem, nfevalsLS) &
          bind(C,name='ARKDlsGetNumRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nfevalsLS
     end function FARKDlsGetNumRhsEvals

     integer(c_int) function FARKDlsGetLastFlag(arkode_mem, flag) &
          bind(C,name='ARKDlsGetLastFlag')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: flag
     end function FARKDlsGetLastFlag

     ! -----------------------------------------------------------------
     ! NOT INTERFACED: ARKDlsGetReturnFlagName
     ! -----------------------------------------------------------------


     ! =================================================================
     ! Interfaces from arkode_spils.h
     ! =================================================================

     integer(c_int) function FARKSpilsSetLinearSolver(arkode_mem, LS) &
          bind(C,name='ARKSpilsSetLinearSolver')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       type(c_ptr), value :: LS
     end function FARKSpilsSetLinearSolver

     ! -----------------------------------------------------------------
     ! Optional inputs to the ARKSPILS linear solver
     ! -----------------------------------------------------------------

     integer(c_int) function FARKSpilsSetEpsLin(arkode_mem, eplifac) &
          bind(C,name='ARKSpilsSetEpsLin')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       real(c_double), value :: eplifac
     end function FARKSpilsSetEpsLin

     integer(c_int) function FARKSpilsSetPreconditioner(arkode_mem, pset, psolve) &
          bind(C,name='ARKSpilsSetPreconditioner')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       type(c_funptr), value :: pset
       type(c_funptr), value :: psolve
     end function FARKSpilsSetPreconditioner

     integer(c_int) function FARKSpilsSetJacTimes(arkode_mem, jtsetup, jtimes) &
          bind(C,name='ARKSpilsSetJacTimesVecFn')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: arkode_mem
       type(c_funptr), value :: jtsetup
       type(c_funptr), value :: jtimes
     end function FARKSpilsSetJacTimes

     ! -----------------------------------------------------------------
     ! Optional outputs from the ARKSPILS linear solver
     ! -----------------------------------------------------------------

     integer(c_int) function FARKSpilsGetWorkSpace(arkode_mem, lenrwLS, leniwLS) &
          bind(C,name='ARKSpilsGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FARKSpilsGetWorkSpace

     integer(c_int) function FARKSpilsGetNumPrecEvals(arkode_mem, npevals) &
          bind(C,name='ARKSpilsGetNumPrecEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: npevals
     end function FARKSpilsGetNumPrecEvals

     integer(c_int) function FARKSpilsGetNumPrecSolves(arkode_mem, npsolves) &
          bind(C,name='ARKSpilsGetNumPrecSolves')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: npsolves
     end function FARKSpilsGetNumPrecSolves

     integer(c_int) function FARKSpilsGetNumLinIters(arkode_mem, nliters) &
          bind(C,name='ARKSpilsGetNumLinIters')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nliters
     end function FARKSpilsGetNumLinIters

     integer(c_int) function FARKSpilsGetNumConvFails(arkode_mem, nlcfails) &
          bind(C,name='ARKSpilsGetNumConvFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nlcfails
     end function FARKSpilsGetNumConvFails

     integer(c_int) function FARKSpilsGetNumJTSetupEvals(arkode_mem, njtsetups) &
          bind(C,name='ARKSpilsGetNumJTSetupEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: njtsetups
     end function FARKSpilsGetNumJTSetupEvals

     integer(c_int) function FARKSpilsGetNumJtimesEvals(arkode_mem, njvevals) &
          bind(C,name='ARKSpilsGetNumJtimesEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: njvevals
     end function FARKSpilsGetNumJtimesEvals

     integer(c_int) function FARKSpilsGetNumRhsEvals(arkode_mem, nfevalsLS) &
          bind(C,name='ARKSpilsGetNumRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nfevalsLS
     end function FARKSpilsGetNumRhsEvals

     integer(c_int) function FARKSpilsGetLastFlag(arkode_mem, flag) &
          bind(C,name='ARKSpilsGetLastFalg')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: flag
     end function FARKSpilsGetLastFlag

     ! -----------------------------------------------------------------
     ! NOT INTERFACED: ARKSpilsGetReturnFlagName
     ! -----------------------------------------------------------------

     ! =================================================================
     ! Interfaces for arkode_diag.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! ARKDiag
     ! -----------------------------------------------------------------

     integer(c_int) function FARKDiag(arkode_mem) &
          bind(C,name='ARKDiag')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
     end function FARKDiag

     ! -----------------------------------------------------------------
     ! Optional outputs from the ARKDIAG linear solver
     ! -----------------------------------------------------------------

     integer(c_int) function FARKDiagGetWorkSpace(arkode_mem, lenrwLS, leniwLS) &
          bind(C,name='ARKDiagGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FARKDiagGetWorkSpace

     integer(c_int) function FARKDiagGetNumRhsEvals(arkode_mem, nfevalsLS) &
          bind(C,name='ARKDiagGetNumRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nfevalsLS
     end function FARKDiagGetNumRhsEvals

     integer(c_int) function FARKDiagGetLastFlag(arkode_mem, flag) &
          bind(C,name='ARKDiagGetLastFlag')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: flag
     end function FARKDiagGetLastFlag

     ! -----------------------------------------------------------------
     ! NOT INTERFACED: ARKDiagGetReturnFlagName
     ! -----------------------------------------------------------------

     ! =================================================================
     ! Interfaces from arkode_bandpre.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! ARKBandPrecInit
     ! -----------------------------------------------------------------

     integer(c_int) function FARKBandPrecInit(arkode_mem, N, mu, ml) &
          bind(C,name='ARKBandPrecInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: arkode_mem
       integer(c_long), value :: N
       integer(c_long), value :: mu
       integer(c_long), value :: ml
     end function FARKBandPrecInit

     ! -----------------------------------------------------------------
     ! Optional output functions : ARKBandPrecGet*
     ! -----------------------------------------------------------------

     integer(c_int) function FARKBandPrecGetWorkSpace(arkode_mem, lenrwLS, leniwLS) &
          bind(C,name='ARKBandPrecGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FARKBandPrecGetWorkSpace

     integer(c_int) function FARKBandPrecGetNumRhsEvals(arkode_mem, nfevalsBP) &
          bind(C,name='ARKBandPrecGetNumRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: nfevalsBP
     end function FARKBandPrecGetNumRhsEvals

     ! =================================================================
     ! Interfaces for arkode_bbdpre.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! ARKBBDPrecInit
     ! -----------------------------------------------------------------

     integer(c_int) function FARKBBDPrecInit(arkode_mem, Nlocal, mudq, mldq, &
          mukeep, mlkeep, dqrely, gloc, cfn) &
          bind(C,name='ARKBBDPrecInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: arkode_mem
       integer(c_long), value :: Nlocal
       integer(c_long), value :: mudq
       integer(c_long), value :: mldq
       integer(c_long), value :: mukeep
       integer(c_long), value :: mlkeep
       real(c_double),  value :: dqrely
       type(c_funptr),  value :: gloc
       type(c_funptr),  value :: cfn
     end function FARKBBDPrecInit

     ! -----------------------------------------------------------------
     ! ARKBBDPrecReInit
     ! -----------------------------------------------------------------

     integer(c_int) function FARKBBDPrecReInit(arkode_mem, mudq, mldq, dqrely) &
          bind(C,name='ARKBBNPrecReInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: arkode_mem
       integer(c_long), value :: mudq
       integer(c_long), value :: mldq
       real(c_double),  value :: dqrely
     end function FARKBBDPrecReInit

     ! -----------------------------------------------------------------
     ! BBDPRE optional output extraction routines
     ! -----------------------------------------------------------------

     integer(c_int) function FARKBBDPrecGetWorkSpace(arkode_mem, lenrwLS, leniwLS) &
          bind(C,name='ARKBBDPrecGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FARKBBDPrecGetWorkSpace

     integer(c_int) function FARKBBDPrecGetNumGfnEvals(arkode_mem, ngevalsBBDP) &
          bind(C,name='ARKBBDPrecGetNumGfnEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: arkode_mem
       integer(c_long)    :: ngevalsBBDP
     end function FARKBBDPrecGetNumGfnEvals

  end interface

end module farkode_mod
