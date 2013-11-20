subroutine wrapper()

  use BoxLib
  use f2kcli
  use bl_IO_module
  use ml_layout_module
  use mg_module
  use box_util_module
  use mt19937_module
  use bl_timer_module

  use    cc_rhs_module
  use nodal_rhs_module

  implicit none

  interface
     subroutine t_cc_ml_multigrid(mla, mgt, rh, coeffs_type, domain_bc, &
                                  do_diagnostics, eps, stencil_order, fabio)
       use mg_module    
       use ml_boxarray_module    
       use ml_layout_module    
       type(ml_layout  ), intent(inout) :: mla
       type(mg_tower) , intent(inout) :: mgt(:)
       type(multifab) , intent(inout) :: rh(:)
       integer        , intent(in   ) :: coeffs_type
       integer        , intent(in   ) :: domain_bc(:,:)
       integer        , intent(in   ) :: do_diagnostics
       real(dp_t)     , intent(in   ) :: eps
       integer        , intent(in   ) :: stencil_order
       logical        , intent(in   ) :: fabio
     end subroutine t_cc_ml_multigrid

     subroutine t_nodal_ml_multigrid(mla, mgt, rh, coeffs_type, domain_bc, &
          do_diagnostics, eps, fabio, stencil_type)
       use mg_module    
       use ml_boxarray_module    
       use ml_layout_module    
       type(ml_layout  ), intent(inout) :: mla
       type(mg_tower) , intent(inout) :: mgt(:)
       type(multifab) , intent(inout) :: rh(:)
       integer        , intent(in   ) :: coeffs_type
       integer        , intent(in   ) :: domain_bc(:,:)
       integer        , intent(in   ) :: do_diagnostics
       real(dp_t)     , intent(in   ) :: eps
       logical        , intent(in   ) :: fabio
       integer        , intent(in   ) :: stencil_type
     end subroutine t_nodal_ml_multigrid
  end interface

  logical :: need_inputs
  logical :: tl
  integer :: narg, farg
  real(dp_t) :: h_finest, wce, wcb

  type(mg_tower   ) :: mgt_default
  type(ml_boxarray) :: mba
  type(ml_layout  ) :: mla

  type(mg_tower), allocatable :: mgt(:)
  type(boxarray)  ::  ba
  integer, allocatable :: ref_ratio(:)
  integer, allocatable :: domain_bc(:,:)
  real(dp_t), allocatable :: dh(:,:)

  character(len=128) :: fname
  character(len=128) :: probin_env

  integer :: stencil_order
  integer :: do_diagnostics
  integer :: un, ierr
  logical :: lexist

  logical, allocatable :: nodal(:)
  logical, allocatable :: pmask(:)
  logical :: nodal_in, dense_in, lininterp

  integer :: i,n,nlevs,ns,dm

  type(box) :: pd

  ! MG solver defaults
  integer :: bottom_solver, bottom_max_iter
  integer :: bottom_solver_in, bottom_max_iter_in
  real(dp_t) :: bottom_solver_eps
  real(dp_t) :: eps
  integer :: max_iter
  integer :: min_width
  integer :: max_nlevel, max_nlevel_in
  integer :: max_lev_of_mba
  integer :: verbose, cg_verbose, ptype
  integer :: nu1, nu2, nub, nuf, solver, smoother
  integer :: ng, nc
  character(len=128) :: test_set
  integer :: test, test_lev, cycle_type, rhs_type, coeffs_type
  logical :: test_set_mglib
  logical :: test_set_hgproj
  logical :: test_random_boxes

  integer :: random_min_size, random_max_size
  integer :: random_blocking_factor, random_num_boxes, random_iseed
  integer :: random_rr
  
  integer :: pd_xyz(MAX_SPACEDIM)
  logical :: pd_pmask(MAX_SPACEDIM)
  integer :: ba_maxsize

  integer :: bcx_lo,bcx_hi,bcy_lo,bcy_hi,bcz_lo,bcz_hi

  integer :: stencil_type

  integer :: bc_all

  logical :: fabio

  type(multifab), allocatable :: rh(:)

  namelist /probin/ cycle_type
  namelist /probin/ rhs_type
  namelist /probin/ coeffs_type
  namelist /probin/ test
  namelist /probin/ nodal_in
  namelist /probin/ dense_in
  namelist /probin/ ba_maxsize
  namelist /probin/ pd_xyz
  namelist /probin/ pd_pmask
  namelist /probin/ h_finest
  namelist /probin/ bcx_lo, bcx_hi
  namelist /probin/ bcy_lo, bcy_hi
  namelist /probin/ bcz_lo, bcz_hi
  namelist /probin/ bc_all
  namelist /probin/ dm
  namelist /probin/ test_set
  namelist /probin/ test_set_mglib
  namelist /probin/ test_set_hgproj
  namelist /probin/ test_lev
  namelist /probin/ fabio

  namelist /probin/ test_random_boxes
  namelist /probin/ random_blocking_factor, random_min_size
  namelist /probin/ random_max_size, random_num_boxes, random_iseed
  namelist /probin/ random_rr

  namelist /probin/ eps, max_iter
  namelist /probin/ nu1, nu2, nub, nuf
  namelist /probin/ bottom_solver, bottom_solver_eps, bottom_max_iter
  namelist /probin/ solver, smoother
  namelist /probin/ min_width, max_nlevel
  namelist /probin/ stencil_order
  namelist /probin/ max_lev_of_mba

  namelist /probin/ verbose, cg_verbose, ptype
  namelist /probin/ do_diagnostics

  !! Defaults:

  need_inputs = .TRUE.
  fabio       = .FALSE.

  max_lev_of_mba = Huge(max_lev_of_mba)

  test           = 0
  cycle_type     = 3 ! Default to V-cycle 
  rhs_type       = 4 ! Default to sum of sin's
  coeffs_type    = 0 ! Default to constant coefficients = 1
  nodal_in       = .false.
  dense_in       = .false.
  lininterp      = .true.

  ba_maxsize     = 32
  pd_xyz         = 32
  pd_pmask       = .FALSE.
  h_finest       = 1.0_dp_t

  bc_all = BC_UNDEF
! bcx_lo         = BC_NEU
! bcy_lo         = BC_NEU
! bcz_lo         = BC_NEU
! bcx_hi         = BC_NEU
! bcy_hi         = BC_NEU
! bcz_hi         = BC_NEU

  bcx_lo         = BC_DIR
  bcy_lo         = BC_DIR
  bcz_lo         = BC_DIR
  bcx_hi         = BC_DIR
  bcy_hi         = BC_DIR
  bcz_hi         = BC_DIR

! bcx_lo         = BC_NEU
! bcx_hi         = BC_NEU
! bcy_lo         = BC_DIR
! bcy_hi         = BC_DIR

  dm             = 2

  test_set_mglib = .FALSE.
  test_set_hgproj = .FALSE.
  test_lev       = 0
  test_set    = ''

  ng                = mgt_default%ng
  nc                = mgt_default%nc
  smoother          = mgt_default%smoother
  nu1               = mgt_default%nu1
  nu2               = mgt_default%nu2
  nub               = mgt_default%nub
  nuf               = mgt_default%nuf
  bottom_solver     = mgt_default%bottom_solver
  bottom_max_iter   = mgt_default%bottom_max_iter
  bottom_solver_eps = mgt_default%bottom_solver_eps
  max_iter          = mgt_default%max_iter
  max_nlevel        = mgt_default%max_nlevel
  min_width         = mgt_default%min_width
  eps               = mgt_default%eps
  verbose           = mgt_default%verbose
  ptype             = mgt_default%ptype
  cg_verbose        = mgt_default%cg_verbose

  do_diagnostics    = 0

  stencil_order = 3

  test_random_boxes = .False.
  random_blocking_factor = 8
  random_min_size =  8
  random_max_size = 32
  random_num_boxes = 1
  random_iseed = 1
  random_rr = 2

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

        case ('--cycle_type')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname,*) cycle_type

        case ('--rhs_type')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname,*) rhs_type

        case ('--coeffs_type')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname,*) coeffs_type

        case ('--ptype')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname,*) ptype

        case ('--verbose')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) verbose

        case ('--cg_verbose')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) cg_verbose

        case ('--dim')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) dm

        case ('--solver')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) solver
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
           read(fname, *) tl
           pd_pmask = tl

        case ('--nodal')
           nodal_in = .true.
        case ('--cc')
           nodal_in = .false.
        case ('--lininterp')
           lininterp = .true.

        case ('--pc')
           !
           ! An easy way to turn off lininterp.
           !
           lininterp = .false.

        case ('--h_finest')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) h_finest

        case ('--dense')
           dense_in = .true.

        case ('--bcx_lo')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcx_lo
        case ('--bcy_lo')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcy_lo
        case ('--bcz_lo')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcz_lo
        case ('--bcx_hi')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcx_hi
        case ('--bcy_hi')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcy_hi
        case ('--bcz_hi')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcz_hi

        case ('--dirichlet')
           bc_all = BC_DIR
        case ('--neumann')
           bc_all = BC_NEU
        case ('--periodic')
           bc_all = BC_PER
           pd_pmask  = .TRUE.

        case ('--fabio')
           fabio = .true.

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

        case ('--do_diagnostics')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) do_diagnostics

        case ('--max_lev_of_mba')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) max_lev_of_mba

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
        case ('--random_rr')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname,*) random_rr

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

     if ( test_set_hgproj ) then
        call read_a_hgproj_grid(mba, test_set, max_lev_of_mba)
     else if ( test_set_mglib ) then
        call read_a_mglib_grid(mba, test_set)
     else
        if ( parallel_ioprocessor() ) &
            print *,'READING ',test_set
        call ml_boxarray_read_boxes(mba, test_set)
     end if

  else

     pd = make_box((/(0,i=1,dm)/), pd_xyz(1:dm)-1)
     if ( test_random_boxes ) then
        call init_genrand(random_iseed)
        if ( max_lev_of_mba == 2 ) then
           call build(mba, 2, dm)
           mba%rr(1,:) = random_rr
           mba%pd(1) = pd
           call boxarray_build_bx(mba%bas(1), pd)
           call boxarray_maxsize(mba%bas(1), ba_maxsize)
           pd = refine(pd, random_rr)
           mba%pd(2) = pd
           call build_random_boxarray(mba%bas(2), pd, &
                random_num_boxes, random_min_size, &
                random_max_size, random_blocking_factor)
        else 
           call build_random_boxarray(ba, pd, &
                random_num_boxes, random_min_size, &
                random_max_size, random_blocking_factor)
           call ml_boxarray_build_ba(mba, ba, pd)
           call destroy(ba)
        end if
     else
        call build(ba, pd)
        call boxarray_maxsize(ba, ba_maxsize)
        call build(mba, ba, pd)
        call destroy(ba)
     end if

     do i = 1, mba%nlevel
        call boxarray_simplify(mba%bas(i))
        call boxarray_maxsize(mba%bas(i), ba_maxsize)
     end do
  end if

  ! For sanity make sure the mba is clean'
  if ( .not. ml_boxarray_clean(mba) ) call bl_error("MBOXARRAY is not 'clean'")

  if ( .not. ml_boxarray_properly_nested(mba) ) call bl_error("MBOXARRAY is not 'properly nested'")

  dm = mba%dim
  
  if ( dm <= 0 ) then
      call bl_error("DM <= 0", dm)
  end if

  allocate(nodal(dm))
  allocate(pmask(dm))

  nlevs = min(max_lev_of_mba, mba%nlevel)

  if ( nodal_in ) then
     nodal = .true.
  else
     nodal = .false.
  end if

  if ( nodal_in ) then
     if ( dense_in ) then
        stencil_type = ND_DENSE_STENCIL
     else
        stencil_type = ND_CROSS_STENCIL
     end if
  else
     stencil_type = CC_CROSS_STENCIL
  end if

  pmask = pd_pmask(1:dm)

  ! Put the bc values from the inputs file into domain_bc
  allocate(domain_bc(dm,2))

  if ( bc_all /= BC_UNDEF ) then
     domain_bc = bc_all
  else
     domain_bc(1,1) = bcx_lo
     domain_bc(1,2) = bcx_hi
     if (dm > 1) then
        domain_bc(2,1) = bcy_lo
        domain_bc(2,2) = bcy_hi
     end if
     if (dm > 2) then
        domain_bc(3,1) = bcz_lo
        domain_bc(3,2) = bcz_hi
     end if
  end if
  where(pmask) 
     domain_bc(:,1) = BC_PER
     domain_bc(:,2) = BC_PER
  end where

  allocate(mgt(nlevs))
  allocate(ref_ratio(dm))

  allocate(dh(nlevs,dm))

  ! Define the grid spacing on the finest level to be 1
  dh(nlevs,:) = h_finest

  if ( nodal_in ) then
     if (dense_in) then
       if (dm .eq. 3) then
         if ( (dh(nlevs,1).eq.dh(nlevs,2)) .and. (dh(nlevs,1) .eq. dh(nlevs,3)) ) then
           ns = 21
         else
           ns = 27
         end if
       else 
         ns = 9
       end if
     else 
       ns = 1 + 2*dm
     end if
  else
     ns = 1 + dm*3
  end if

  call ml_layout_build(mla, mba, pmask = pmask)

  do n = nlevs, 1, -1

     if (n == 1) then

        max_nlevel_in = max_nlevel
        bottom_solver_in = bottom_solver
        bottom_max_iter_in = bottom_max_iter

     else

        ref_ratio = mba%rr(n-1,:)
        if ( parallel_IOProcessor() ) then
            print *, "REF_RATIO = ", ref_ratio
        end if

        if ( all(ref_ratio == 2) ) then
           max_nlevel_in = 1
        else if ( all(ref_ratio == 4) ) then
           max_nlevel_in = 2
        else
           call bl_error("WRAPPER: confused about ref_ratio")
        end if
        bottom_solver_in = 0
        bottom_max_iter_in = nu1

     end if

     if ( nodal_in ) min_width = 2

     ! Define the grid spacing on all lower levels in terms of the finest level.
     if (n < nlevs) then
        dh(n,:) = dh(n+1,:) * mba%rr(n,:)
     end if

     call mg_tower_build(mgt(n), mla%la(n), mba%pd(n), domain_bc, stencil_type,&
          dh = dh(n,:), &
          ns = ns, &
          smoother = smoother, &
          nu1 = nu1, &
          nu2 = nu2, &
          nub = nub, &
          nuf = nuf, &
          cycle_type = cycle_type, &
          bottom_solver = bottom_solver_in, &
          bottom_max_iter = bottom_max_iter_in, &
          bottom_solver_eps = bottom_solver_eps, &
          max_iter = max_iter, &
          max_nlevel = max_nlevel_in, &
          min_width = min_width, &
          eps = eps, &
          verbose = verbose, &
          cg_verbose = cg_verbose, &
          ptype = ptype, &
          use_lininterp = lininterp, &
          nodal = nodal)

  end do

  ! Allocate space for the RHS for the solve.
  allocate(rh(nlevs))

  call wall_second(wcb)
  if ( nodal_in ) then
     do n = nlevs, 1, -1
        call multifab_build( rh(n), mla%la(n), 1, 1, nodal)
     end do
     call nodal_rhs(mla, rh)
     call t_nodal_ml_multigrid(mla, mgt, rh, coeffs_type, domain_bc, do_diagnostics, eps, &
                               fabio, stencil_type)
  else
     do n = nlevs, 1, -1
        call multifab_build( rh(n), mla%la(n), 1, 0, nodal)
     end do
     call cc_rhs(mla, pd, rh, rhs_type)
     call t_cc_ml_multigrid(mla, mgt, rh, coeffs_type, domain_bc, do_diagnostics, eps, stencil_order, fabio)
  end if

  call wall_second(wce)

  wce = wce - wcb

  if ( parallel_ioprocessor() ) then
     print*, 'Wall clock for solve: ', wce
  end if

  do n = nlevs, 1, -1
     call multifab_destroy(rh(n))
  end do
  deallocate(rh)

  do n = 1, nlevs
     call mg_tower_destroy(mgt(n))
  end do

  call destroy(mba)
  call destroy(mla)

end subroutine wrapper
