subroutine t_cc_multigrid()
  use BoxLib
  use f2kcli
  use stencil_module
  use coeffs_module
  use mg_module
  use list_box_module
  use ml_boxarray_module
  use itsol_module
  use sparse_solve_module
  use bl_mem_stat_module
  use bl_timer_module
  use box_util_module
  use bl_IO_module
  use fabio_module
  use mt19937_module
  use bl_prof_module

  implicit none

  logical, allocatable :: pmask(:)

  integer, allocatable :: domain_bc(:,:)

  type(boxarray) ba, pdv

  type(box) pd

  type(multifab), allocatable :: coeffs(:)

  type(multifab) :: uu, rh, ss
  type(imultifab) :: mm
  type(sparse) :: sparse_object
  type(ml_boxarray) :: mba
  type(layout) :: la
  type(mg_tower) :: mgt
  real(dp_t), allocatable :: dh(:)
  integer i, dm, ns

  integer :: narg, farg
  character(len=128) :: fname
  integer :: ba_maxsize, pd_xyz(MAX_SPACEDIM)
  logical :: pd_pmask(MAX_SPACEDIM)
  integer :: test
  real(dp_t) :: snrm(2)

  real(dp_t), allocatable :: qq(:)
  logical :: qq_history
  character(len=128) :: qq_fname
  integer :: qq_un, qq_num_iter

  ! MG solver defaults
  integer :: bottom_solver, bottom_max_iter
  real(dp_t) :: bottom_solver_eps
  real(dp_t) :: eps
  real(dp_t) :: h_finest
  integer :: max_iter
  integer :: min_width
  integer :: max_nlevel
  integer :: verbose
  integer :: nu1, nu2, nub, nuf, gamma, cycle, solver, smoother
  real(dp_t) :: omega
  integer :: ng, nc
  character(len=128) :: test_set
  integer :: test_lev
  logical :: test_set_mglib
  logical :: test_set_hgproj
  logical :: test_random_boxes
  integer :: stat
  logical :: neumann
  real(dp_t), allocatable :: tmval(:)

  integer :: random_min_size, random_max_size
  integer :: random_blocking_factor, random_num_boxes, random_iseed

  type(timer) :: tm(2)

  character(len=128) :: probin_env
  integer :: un, ierr

  logical :: mem_check
  logical :: ldummy

  logical :: lexist
  logical :: need_inputs
  logical :: fabio

  integer :: stencil_order
  logical, allocatable :: nodal(:)
  logical uu_rand

  real(dp_t), allocatable :: xa(:), xb(:), pxa(:), pxb(:)
  real(dp_t), parameter :: ONE = 1.0_dp_t
  real(dp_t), parameter :: ZERO = 0.0_dp_t

  character(len=128) defect_dirname
  logical :: defect_history

  logical :: make_solvable
  integer :: ba_coarsen
  type(bl_prof_timer), save :: bpt_setup

  namelist /probin/ test
  namelist /probin/ ba_maxsize
  namelist /probin/ pd_xyz
  namelist /probin/ pd_pmask
  namelist /probin/ dm
  namelist /probin/ test_set
  namelist /probin/ test_set_mglib
  namelist /probin/ test_set_hgproj
  namelist /probin/ test_lev
  namelist /probin/ fabio

  namelist /probin/ h_finest

  namelist /probin/ test_random_boxes
  namelist /probin/ random_blocking_factor, random_min_size
  namelist /probin/ random_max_size, random_num_boxes, random_iseed

  namelist /probin/ eps, max_iter
  namelist /probin/ nu1, nu2, nub, nuf, gamma, cycle
  namelist /probin/ bottom_solver, bottom_solver_eps, bottom_max_iter
  namelist /probin/ solver, smoother
  namelist /probin/ min_width, max_nlevel
  namelist /probin/ stencil_order

  namelist /probin/ mem_check

  namelist /probin/ verbose
  namelist /probin/ neumann

  namelist /probin/ qq_fname
  namelist /probin/ qq_history

  namelist /probin/ uu_rand

  namelist /probin/ defect_dirname
  namelist /probin/ ba_coarsen

  !! Defaults:

  mem_check = .true.

  neumann = .FALSE.
  need_inputs = .TRUE.

  uu_rand = .false.

  test           = 0
  fabio         = .false.

  ba_maxsize     = 16
  pd_xyz         = 32

  pd_pmask      = .false.

  dm             = 2

  h_finest       = 1.0_dp_t

  test_set_mglib = .FALSE.
  test_set_hgproj = .FALSE.
  test_lev       = 0
  test_set    = ''

  test_random_boxes = .False.
  random_blocking_factor = 8
  random_min_size =  8
  random_max_size = 32
  random_num_boxes = 1
  random_iseed = 1

  ba_coarsen = 1

  ng                = mgt%ng
  nc                = mgt%nc
  max_nlevel        = mgt%max_nlevel
  max_iter          = mgt%max_iter
  eps               = mgt%eps
  smoother          = mgt%smoother
  nu1               = mgt%nu1
  nu2               = mgt%nu2
  nub               = mgt%nub
  nuf               = mgt%nuf
  gamma             = mgt%gamma
  omega             = mgt%omega
!  cycle             = mgt%cycle
  bottom_solver     = mgt%bottom_solver
  bottom_solver_eps = mgt%bottom_solver_eps
  bottom_max_iter   = mgt%bottom_max_iter
  min_width         = mgt%min_width
  verbose           = mgt%verbose

  qq_fname = ''
  qq_history = .false.

  defect_history = .false.
  defect_dirname = ''

  stencil_order = 3

  narg = command_argument_count()

  call get_environment_variable('PROBIN', probin_env, status = ierr)
  if ( need_inputs .AND. ierr == 0 ) then
     un = unit_new()
     open(unit=un, file = probin_env, status = 'old', action = 'read')
     read(unit=un, nml = probin)
     close(unit=un)
     need_inputs = .FALSE.
  end if

  farg = 1
  if ( need_inputs .AND. narg >= 1 ) then
     call get_command_argument(farg, value = fname)
     inquire(file = fname, exist = lexist )
     if ( lexist ) then
        farg = farg + 1
        un = unit_new()
        open(unit=un, file = fname, status = 'old', action = 'read')
        read(unit=un, nml = probin)
        close(unit=un)
        need_inputs = .FALSE.
     end if
  end if

  inquire(file = 'inputs', exist = lexist)
  if ( need_inputs .AND. lexist ) then
     un = unit_new()
     open(unit=un, file = 'inputs', status = 'old', action = 'read')
     read(unit=un, nml = probin)
     close(unit=un)
     need_inputs = .FALSE.
  end if

  if ( .true. ) then
     do while ( farg <= narg )
        call get_command_argument(farg, value = fname)
        select case (fname)

        case ('--test')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname,*) test

        case ('--defect_dirname')
           farg = farg + 1
           call get_command_argument(farg, value = defect_dirname)
           defect_history = .true.

        case ('--qq_fname')
           farg = farg + 1
           call get_command_argument(farg, value = qq_fname)
        case ('--qq_history')
           qq_history = .true.

        case ('--uu_rand')
           uu_rand = .true.

        case ('--verbose')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) verbose

        case ('--dim')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) dm

        case ('--solver')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) solver
        case ('--gamma')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname,*) gamma
        case ('--smoother')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) smoother
        case ('--nu1')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) nu1
        case ('--nu2')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) nu2
        case ('--nub')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) nub
        case ('--nuf')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) nuf
!        case ('--cycle')
!           farg = farg + 1
!           call get_command_argument(farg, value = fname)
!           read(fname, *) cycle
        case ('--omega')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) omega

        case ('--h_finest')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) h_finest

        case ('--min_width')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) min_width
        case ('--max_nlevel')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) max_nlevel
        case ('--max_iter')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) max_iter
        case ('--eps')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) eps

        case ('--bottom_solver')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bottom_solver
        case ('--bottom_solver_eps')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bottom_solver_eps
        case ('--bottom_max_iter')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bottom_max_iter

        case ('--maxsize')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) ba_maxsize
        case ('--pd_x')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) pd_xyz(1)
        case ('--pd_y')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) pd_xyz(2)
        case ('--pd_z')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) pd_xyz(3)
        case ('--pd_xyz')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) i
           pd_xyz = i

        case ('--pd_pmask_x')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) pd_pmask(1)
        case ('--pd_pmask_y')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) pd_pmask(2)
        case ('--pd_pmask_z')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) pd_pmask(3)
        case ('--pd_pmask')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) ldummy
           pd_pmask = ldummy

        case ('--fabio')
           fabio = .True.

        case ('--neumann')
           neumann = .True.

        case ('--test_random_boxes')
           test_random_boxes = .True.
        case ('--random_blocking_factor')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname,*) random_blocking_factor
        case ('--random_min_size')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname,*) random_min_size
        case ('--random_max_size')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname,*) random_max_size
        case ('--random_num_boxes')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname,*) random_num_boxes
        case ('--random_iseed')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname,*) random_iseed

        case ('--test_set')
           farg = farg + 1
           call get_command_argument(farg, value = test_set)
        case ('--test_lev')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) test_lev
        case ('--test_set_mglib')
           test_set_mglib = .True.
        case ('--test_set_hgproj')
           test_set_hgproj = .True.

        case ('--stencil_order')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) stencil_order

        case ('--ba_coarsen')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) ba_coarsen

        case ('--')
           farg = farg + 1
           exit

        case default
           if ( .not. parallel_q() ) then
              write(*,*) 'UNKNOWN option = ', fname
              call bl_error("MAIN")
           end if
        end select

        farg = farg + 1
     end do
  end if

  if ( test_set /= '' ) then
     if ( test_set_mglib ) then
        call read_a_mglib_grid(mba, test_set)
     else if ( test_set_hgproj ) then
        call read_a_hgproj_grid(mba, test_set)
     else
        call ml_boxarray_read_boxes(mba, test_set)
     end if
     dm = mba%dim
     if ( test_lev == 0 ) test_lev = mba%nlevel
     pd = mba%pd(test_lev)
     call boxarray_build_copy(ba, mba%bas(test_lev))
     call ml_boxarray_destroy(mba)
  else
     call box_build_2(pd, (/(0,i=1,dm)/), pd_xyz(1:dm)-1)
     if ( test_random_boxes ) then
        call init_genrand(random_iseed)
        call build_random_boxarray(ba, pd, &
             random_num_boxes, random_min_size, &
             random_max_size, random_blocking_factor)
     else
        call boxarray_build_bx(ba, pd)
        call boxarray_maxsize(ba, ba_maxsize)
     end if
  end if

  if ( ba_coarsen > 1 ) then
     call boxarray_coarsen(ba, ba_coarsen)
  end if

  if ( test /= 0 .AND. max_iter == mgt%max_iter ) then
     max_iter = 1000
  end if

  allocate(pmask(dm))
  pmask = pd_pmask(1:dm)

  allocate(domain_bc(dm,2))
  domain_bc = BC_DIR
  if ( neumann ) domain_bc = BC_NEU
  where(pmask)
     domain_bc(:,1) = BC_PER
     domain_bc(:,2) = BC_PER
  end where

  make_solvable = .true.
  if ( any(domain_bc  == BC_DIR) ) make_solvable = .false.

  if ( parallel_IOProcessor() ) then
     call box_print(pd, 'problem domain')
     print *, 'Problem Domain Volume(ba) ', volume(pd)
     call boxarray_print(ba, "Problem Boxes")
     print *, 'Problem Num Boxes ', nboxes(ba)
     print *, 'Problem Boxes Volume(ba) ', volume(ba)
     print *, 'PMASK: ', pmask
     !$OMP PARALLEL
     !$OMP MASTER
     print *, 'MAX_THREADS = ', omp_get_max_threads()
     !$OMP END MASTER
     !$OMP END PARALLEL
     print *, 'STENCIL_ORDER = ', stencil_order
     print *, 'BA_COARSEN = ', ba_coarsen
  end if

  do i = 1, dm
     if ( pmask(i) ) domain_bc(i,:) = BC_PER
  end do

  call build(la, ba, pd = pd, pmask = pmask)

  allocate(dh(dm))
  dh = h_finest

  call multifab_build(uu, la, 1, 1)
  call multifab_build(rh, la, 1, 0)
  if ( uu_rand ) then
     call mf_rand(uu)
  else
     call setval(uu, ZERO, all=.true.)
  end if
  call setval(rh, val = ZERO, all=.true.)
  call mf_init(rh)
  if ( fabio ) then
     call fabio_multifab_write_d(rh, "tdir1", "rh")
  end if

  ns = 1 + dm*3
  allocate(nodal(dm))
  nodal = .False.

  if ( parallel_IOProcessor() ) then
     print *, 'dh = ', dh
  end if

  call mg_tower_build(mgt, la, pd, domain_bc, &
       dh = dh, &
       ns = ns, &
       smoother = smoother, &
       nu1 = nu1, &
       nu2 = nu2, &
       nuf = nuf, &
       nub = nub, &
       gamma = gamma, &
!       cycle = cycle, &
       omega = omega, &
       bottom_solver = bottom_solver, &
       bottom_max_iter = bottom_max_iter, &
       bottom_solver_eps = bottom_solver_eps, &
       max_iter = max_iter, &
       max_nlevel = max_nlevel, &
       min_width = min_width, &
       eps = eps, &
       verbose = verbose, &
       nodal = nodal)

  allocate(coeffs(mgt%nlevels))
  call multifab_build(coeffs(mgt%nlevels), la, 1+dm, 1)
  call setval(coeffs(mgt%nlevels), ZERO, 1, all=.true.)
  ! call setval(coeffs(mgt%nlevels), ONE, 1, all=.true.)
  call setval(coeffs(mgt%nlevels), ONE, 2, dm, all=.true.)

  allocate(xa(dm), xb(dm), pxa(dm), pxb(dm))
  pxa = ZERO
  pxb = ZERO
  xa = ZERO
  xb = ZERO
  select case (test)
  case (0)
     call timer_start(tm(1))
     do i = mgt%nlevels-1, 1, -1
        call multifab_build(coeffs(i), mgt%ss(i)%la, 1+dm, 1)
        call setval(coeffs(i), ZERO, 1, dm+1, all=.true.)
        call coarsen_coeffs(coeffs(i+1), coeffs(i))
     end do
     do i = mgt%nlevels, 1, -1
        pdv = layout_boxarray(mgt%ss(i)%la)
        call stencil_fill_cc(mgt%ss(i), coeffs(i), mgt%dh(:,i), pdv, &
             mgt%mm(i), xa, xb, pxa, pxb, pd, stencil_order, domain_bc)
call fill_boundary(coeffs(i))
     end do
     if ( bottom_solver == 3 ) then
        call copy(mgt%ss1, mgt%ss(1))
        call copy(mgt%mm1, mgt%mm(1))
        if ( parallel_IOProcessor() ) then
           call sparse_build(mgt%sparse_object, mgt%ss1, mgt%mm1, &
                mgt%ss1%la, stencil_order, verbose)
        end if
     end if
     call timer_stop(tm(1))
     if ( qq_history ) then
        allocate(qq(0:max_iter))
     else
        allocate(qq(0))
     end if
     call timer_start(tm(2))
     call mg_tower_solve(mgt, uu, rh, &
          stat = stat, &
          qq = qq, &
          num_iter = qq_num_iter, &
          defect_history = defect_history, &
          defect_dirname = defect_dirname &
          )
     call timer_stop(tm(2))
     if ( stat /= 0 ) then
        call bl_warn("MG_TOWER_SOLVE: failed : ", stat)
     end if
     do i = mgt%nlevels-1, 1, -1
        call multifab_destroy(coeffs(i))
     end do
  case (1)
     la = mgt%dd(Mgt%nlevels)%la
     pdv = layout_boxarray(la)
     call  multifab_build(ss, la, ns, 0)
     call imultifab_build(mm, la,  1, 0)
     call timer_start(tm(1))
     call build(bpt_setup, "BICG_SETUP")
     call stencil_fill_cc(ss, coeffs(mgt%nlevels), mgt%dh(:,mgt%nlevels), &
          pdv, mm, xa, xb, pxa, pxb, pd, stencil_order, domain_bc)
call fill_boundary(coeffs(mgt%nlevels))
     call timer_stop(tm(1))
     call destroy(bpt_setup)
     call timer_start(tm(2))
     call itsol_BiCGStab_solve(ss, uu, rh, mm, eps, max_iter, verbose, stat)
     call timer_stop(tm(2))
     if ( stat /= 0 ) then
        call bl_warn("itsol_BiCGStab_solve: failed : ", stat)
     end if
  case (2)
     la = mgt%dd(Mgt%nlevels)%la
     pdv = layout_boxarray(la)
     call  multifab_build(ss, la, ns, 0)
     call imultifab_build(mm, la,  1, 0)
     call build(bpt_setup, "CG_SETUP")
     call timer_start(tm(1))
     call stencil_fill_cc(ss, coeffs(mgt%nlevels), mgt%dh(:,mgt%nlevels), &
          pdv, mm, xa, xb, pxa, pxb, pd, stencil_order, domain_bc)
call fill_boundary(coeffs(mgt%nlevels))
     call timer_stop(tm(1))
     call timer_start(tm(2))
     call itsol_CG_solve(ss, uu, rh, mm, eps, max_iter, verbose, stat)
     call timer_stop(tm(2))
     call destroy(bpt_setup)
     if ( stat /= 0 ) then
        call bl_warn("itsol_CG_SOLVE: failed : ", stat)
     end if
  case (3)
     if ( parallel_nprocs() > 1 ) then
        call bl_error("cc_single: bottom solver test can't be called in parallel")
     end if
     la = mgt%dd(Mgt%nlevels)%la
     pdv = layout_boxarray(la)
     call multifab_build(ss, la, ns, 0)
     call imultifab_build(mm, la,  1, 0)
     call build(bpt_setup, "SS_SETUP")
     call timer_start(tm(1))
     call stencil_fill_cc(ss, coeffs(mgt%nlevels), mgt%dh(:,mgt%nlevels), &
          pdv, mm, xa, xb, pxa, pxb, pd, stencil_order, domain_bc)
call fill_boundary(coeffs(mgt%nlevels))
     call timer_stop(tm(1))
     call destroy(bpt_setup)
     call sparse_build(sparse_object, ss, mm, la, stencil_order, verbose)
     call timer_start(tm(2))
     call sparse_solve(sparse_object, uu, rh, eps, max_iter, verbose, stat)
     if ( stat /= 0 ) then
        call bl_warn("SPARSE SOLVE FAILED DUE TO BREAKDOWN")
     end if
     call timer_stop(tm(2))
  end select

  call multifab_fill_boundary(uu)
  ! call print(uu, "SOLN")

  call multifab_destroy(coeffs(mgt%nlevels))
  deallocate(coeffs)

  if ( parallel_IOProcessor() ) print *, 'TIMER RESOLUTION: ', timer_tick()
  call timer_print(tm(1), "setup  time")
  call timer_print(tm(2), "solver time")
  if ( test == 0 ) then
     allocate(tmval(mgt%nlevels))
     tmval(1) = timer_value(mgt%tm(1), TOTAL = .true.)
     do i = 2, mgt%nlevels
        tmval(i) = timer_value(mgt%tm(i), TOTAL = .true. ) &
             - timer_value(mgt%tm(i-1), TOTAL = .true.)
     end do
     if ( parallel_IOProcessor() ) print *, 'MG type per level'
     do i = 1, mgt%nlevels
        call timer_print(mgt%tm(i), TOTAL = .TRUE., ADVANCE = 'no')
        if ( parallel_IOProcessor() ) &
             write(*, '(": ",f15.4)') tmval(i)
     end do
  end if

  if ( mem_check ) then
     if ( parallel_IOProcessor() ) print *, 'MEMORY STATS'
     call print(multifab_mem_stats(),  " multifab before")
     call print(imultifab_mem_stats(), "imultifab before")
     call print(fab_mem_stats(),       "      fab before")
     call print(ifab_mem_stats(),      "     ifab before")
     call print(boxarray_mem_stats(),  " boxarray before")
     call print(boxassoc_mem_stats(),  " boxassoc before")
     call print(layout_mem_stats(),    "   layout before")
  end if

  snrm(1) = norm_l2(uu)
  snrm(2) = norm_inf(uu)
  if ( parallel_IOProcessor() ) &
       print *, 'solution norm = ', snrm(1), "/", snrm(2)

  if ( fabio ) then
     call fabio_multifab_write_d(mgt%ss(mgt%nlevels), "tdir1", "ss")
     call fabio_multifab_write_d(uu, "tdir1", "uu")
  end if

  if ( test == 3 ) then
     call sparse_destroy(sparse_object)
  end if
  if ( test > 0 ) then
     call destroy(ss)
     call destroy(mm)
  end if
  call mg_tower_destroy(mgt)

  call multifab_destroy(uu)
  call multifab_destroy(rh)
  call layout_destroy(la)
  call boxarray_destroy(ba)

  if ( mem_check ) then
     call print(multifab_mem_stats(),  " multifab at end")
     call print(imultifab_mem_stats(), "imultifab at end")
     call print(fab_mem_stats(),       "      fab at end")
     call print(ifab_mem_stats(),      "     ifab at end")
     call print(boxarray_mem_stats(),  " boxarray at end")
     call print(boxassoc_mem_stats(),  " boxassoc at end")
     call print(layout_mem_stats(),    "   layout at end")
  end if

  if ( qq_history ) then
     if ( qq_fname /= '' ) then
        qq_un = unit_new()
        open(unit = qq_un,         &
             file=qq_fname,        &
             status = 'replace',   &
             action = 'write',     &
             form = 'formatted',   &
             access = 'sequential' &
             )
     else
        qq_un = unit_stdout()
     end if
     do i = lbound(qq,1), qq_num_iter
        write(unit=qq_un,fmt='(i3,1x,es30.15)')  i, qq(i)
     end do
  end if

contains

  subroutine mf_init(mf)
    type(multifab), intent(inout) :: mf
    integer i
    type(box) bx
    do i = 1, mf%nboxes; if ( remote(mf,i) ) cycle
       bx = get_box(mf,i)
       bx%lo(1:bx%dim) = (bx%hi(1:bx%dim) + bx%lo(1:bx%dim))/2
       bx%hi(1:bx%dim) = bx%lo(1:bx%dim)
       call setval(mf%fbs(i), ONE, bx)
       if ( make_solvable ) then
          bx%lo(1:bx%dim) = (bx%hi(1:bx%dim) + bx%lo(1:bx%dim))/2 + 1
          bx%hi(1:bx%dim) = bx%lo(1:bx%dim)
          call setval(mf%fbs(i), -ONE, bx)
       end if
    end do
  end subroutine mf_init

  subroutine mf_rand(mf)
    type(multifab), intent(inout) :: mf
    integer :: i
    real(kind=dp_t), pointer :: mp(:,:,:,:)
    do i = 1, mf%nboxes; if ( remote(mf,i) ) cycle
       mp => dataptr(mf, i)
       call mt_random_number(mp)
    end do
  end subroutine mf_rand

  subroutine mf_init_1(mf)
    type(multifab), intent(inout) :: mf
    integer i, ng
    type(box) bx
    real(kind=dp_t), pointer :: mp(:,:,:,:)

    ng = mf%ng
    do i = 1, mf%nboxes; if ( remote(mf,i) ) cycle
       mp => dataptr(mf, i)
       bx = get_ibox(mf, i)
       select case ( mf%dim ) 
       case (2)
          call pot_init_2d(mp(:,:,1,1), bx, ng)
       end select
    end do
  end subroutine mf_init_1

  subroutine pot_init_2d(po, bx, ng)

    type(box), intent(in) :: bx
    integer, intent(in)   :: ng
    real(dp_t), intent(inout) :: po(bx%lo(1)-ng:,bx%lo(2)-ng:)
    integer :: i, j

    po = 0
    do j = bx%lo(2), bx%hi(2)
       do i = bx%lo(1), bx%hi(1)
          po(i,j) = po(i,j) - gau_2d(i, j)
       end do
    end do

  end subroutine pot_init_2d
  function gau_2d(i, j) result(r)
    real(dp_t) :: r
    integer, intent(in) :: i, j
    integer :: n
    real(dp_t) :: atom_str(1)
    real(dp_t) :: atom_wid(1)
    real(dp_t) :: atom_loc(2,1)
    n = 1
    atom_str = 1
    atom_wid = 20
    atom_loc(1,1) = (lwb(pd,1)+upb(pd,1))/2
    atom_loc(2,1) = (lwb(pd,2)+upb(pd,2))/2
    r = atom_str(n)*exp(-((i*dh(1) - atom_loc(1,n))**2 + (j*dh(2) - atom_loc(2,n))**2)/atom_wid(n)**2)
  end function gau_2d

end subroutine t_cc_multigrid
