module mg_module

  use multifab_module
  use stencil_module
  use sparse_solve_module
  use mg_tower_module
  use bl_timer_module

  implicit none

  real(kind=dp_t), parameter, private :: zero = 0.0_dp_t

  interface destroy
     module procedure mg_tower_destroy
  end interface

contains

  recursive subroutine mg_tower_build(mgt, la, pd, domain_bc, &
                            nu1, nu2, nuf, nub, gamma, cycle_type, &
                            smoother, omega, &
                            dh, &
                            ns, &
                            nc, ng, &
                            max_nlevel, max_bottom_nlevel, min_width, &
                            max_iter, eps, abs_eps, &
                            bottom_solver, bottom_max_iter, bottom_solver_eps, &
                            st_type, &
                            verbose, cg_verbose, nodal)

    use bl_IO_module
    use bl_prof_module
    type(mg_tower), intent(inout) :: mgt
    type(layout), intent(in   ) :: la
    type(box), intent(in) :: pd
    integer, intent(in) :: domain_bc(:,:)

    integer, intent(in), optional :: ns
    integer, intent(in), optional :: nc
    integer, intent(in), optional :: ng
    integer, intent(in), optional :: nu1, nu2, nuf, nub, gamma, cycle_type
    integer, intent(in), optional :: smoother
    logical, intent(in), optional :: nodal(:)
    real(dp_t), intent(in), optional :: omega
    real(kind=dp_t), intent(in), optional :: dh(:)
    real(kind=dp_t), intent(in), optional :: eps
    real(kind=dp_t), intent(in), optional :: abs_eps
    real(kind=dp_t), intent(in), optional :: bottom_solver_eps
    integer, intent(in), optional :: max_nlevel
    integer, intent(in), optional :: max_bottom_nlevel
    integer, intent(in), optional :: min_width
    integer, intent(in), optional :: max_iter
    integer, intent(in), optional :: bottom_solver
    integer, intent(in), optional :: bottom_max_iter
    integer, intent(in), optional :: verbose
    integer, intent(in), optional :: cg_verbose
    integer, intent(in), optional :: st_type

    integer :: lo_grid,hi_grid,lo_dom,hi_dom
    integer :: ng_for_res
    integer :: n, i, id
    type(layout) :: la1, la2
    type(boxarray) :: ba
    logical :: nodal_flag
    real(kind=dp_t) :: vol
    type(bl_prof_timer), save :: bpt

    ! These are added to help build the bottom_mgt
    type(  layout)  :: old_coarse_la, new_coarse_la
    type(     box)  :: coarse_pd,bxs
    type(boxarray)  :: new_coarse_ba
    integer         :: bottom_box_size
    real(kind=dp_t) :: coarse_dx(pd%dim)

    call build(bpt, "mgt_build")

    ! Need to set this here because so many things below depend on it.
    mgt%dim = get_dim(la)

    ! Paste in optional arguments
    if ( present(st_type)           ) mgt%st_type           = st_type
    if ( present(ng)                ) mgt%ng                = ng
    if ( present(nc)                ) mgt%nc                = nc
    if ( present(max_nlevel)        ) mgt%max_nlevel        = max_nlevel
    if ( present(max_bottom_nlevel) ) mgt%max_bottom_nlevel = max_bottom_nlevel
    if ( present(max_iter)          ) mgt%max_iter          = max_iter
    if ( present(eps)               ) mgt%eps               = eps
    if ( present(abs_eps)           ) mgt%abs_eps           = abs_eps
    if ( present(smoother)          ) mgt%smoother          = smoother
    if ( present(nu1)               ) mgt%nu1               = nu1
    if ( present(nu2)               ) mgt%nu2               = nu2
    if ( present(nuf)               ) mgt%nuf               = nuf
    if ( present(nub)               ) mgt%nub               = nub
    if ( present(gamma)             ) mgt%gamma             = gamma
    if ( present(omega)             ) mgt%omega             = omega
    if ( present(cycle_type)        ) mgt%cycle_type        = cycle_type
    if ( present(bottom_solver)     ) mgt%bottom_solver     = bottom_solver
    if ( present(bottom_solver_eps) ) mgt%bottom_solver_eps = bottom_solver_eps
    if ( present(bottom_max_iter)   ) mgt%bottom_max_iter   = bottom_max_iter
    if ( present(min_width)         ) mgt%min_width         = min_width
    if ( present(verbose)           ) mgt%verbose           = verbose
    if ( present(cg_verbose)        ) mgt%cg_verbose        = cg_verbose

    nodal_flag = .false.
    if ( present(nodal) ) then
       if ( all(nodal) ) then
          nodal_flag = .true.
       else if ( any(nodal) ) then
          call bl_error("mixed nodal/not nodal forbidden")
       end if
    end if

    if ( present(ns) ) then
       mgt%ns = ns
    else
       if ( nodal_flag ) then
          mgt%ns = 3**mgt%dim
       else
          mgt%ns = 1 + 3*mgt%dim
       end if
    end if

    if ( .not. nodal_flag ) then
       if ( .not. present(omega) ) then
          select case ( mgt%dim )
          case (1)
             mgt%omega = 0.6_dp_t
          case (2)
             select case (mgt%smoother)
             case ( MG_SMOOTHER_JACOBI )
                mgt%omega = 4.0_dp_t/5.0_dp_t
             end select
          case (3)
             select case (mgt%smoother)
             case ( MG_SMOOTHER_JACOBI )
                mgt%omega = 6.0_dp_t/7.0_dp_t
             case ( MG_SMOOTHER_GS_RB )
                mgt%omega = 1.15_dp_t
             end select
          end select
       end if
    end if
    ng_for_res = 0; if ( nodal_flag ) ng_for_res = 1

    n = max_mg_levels(get_boxarray(la), mgt%min_width)
    mgt%nlevels = min(n,mgt%max_nlevel) 

    n = mgt%nlevels
    mgt%nboxes = nboxes(la)

    allocate(mgt%face_type(mgt%nboxes,mgt%dim,2))
    allocate(mgt%skewed(mgt%nlevels,mgt%nboxes))
    allocate(mgt%skewed_not_set(mgt%nlevels))
    mgt%skewed_not_set = .true.

    allocate(mgt%cc(n), mgt%ff(n), mgt%dd(n), mgt%uu(n-1), mgt%ss(n), mgt%mm(n))
    allocate(mgt%pd(n),mgt%dh(mgt%dim,n))
    allocate(mgt%tm(n))
    if ( n == 1 ) then
       mgt%tm(n)%name = "SINGLE LEVEL"
    else
       mgt%tm(n)%name = "FINEST LEVEL"
       mgt%tm(1)%name = "COARSEST LEVEL"
    end if

    la1 = la

    do i = n, 1, -1
       call  multifab_build(mgt%cc(i), la1, mgt%nc, ng_for_res, nodal)
       call  multifab_build(mgt%ff(i), la1, mgt%nc, ng_for_res, nodal)
       call  multifab_build(mgt%dd(i), la1, mgt%nc, ng_for_res, nodal)
       call  multifab_build(mgt%ss(i), la1, ns, 0, nodal)

       call setval(mgt%cc(i), zero, all = .TRUE.)
       call setval(mgt%ff(i), zero, all = .TRUE.)
       call setval(mgt%dd(i), zero, all = .TRUE.)
       call setval(mgt%ss(i), zero, all = .TRUE.)

       call imultifab_build(mgt%mm(i), la1, 1, 0, nodal)
       if ( i /= n ) &
          call multifab_build(mgt%uu(i), la1, mgt%nc, mgt%ng, nodal)
       
       mgt%pd(i) = coarsen(pd, 2**(n-i))
       if ( i > 1 ) call layout_build_coarse(la2, la1, (/(2,i=1,mgt%dim)/))
       la1 = la2
    end do

    if ( n > 1 .and. mgt%bottom_solver == 3 ) then
       la2 = mgt%cc(1)%la
       call layout_build_derived(la1, la2)
       call build(mgt%rh1, la1, mgt%nc, 0, nodal)
       call build(mgt%uu1, la1, mgt%nc, 0, nodal)
       call build(mgt%ss1, la1, ns,     0, nodal)
       call build(mgt%mm1, la1, 1,      0, nodal)
    else if ( nodal_flag .and. (mgt%bottom_solver == 1 .or. mgt%bottom_solver == 2) ) then
       la2 = mgt%cc(1)%la
       call layout_build_derived(la1, la2)
       call build_nodal_dot_mask(mgt%nodal_mask,mgt%ss(1))
    end if

    mgt%uniform_dh = .true.
    if ( present(dh) ) then
       mgt%dh(:,mgt%nlevels) = dh(:)
       select case (mgt%dim)
       case (2)
          mgt%uniform_dh = (dh(1) == dh(2))
       case (3)
          mgt%uniform_dh = ((dh(1) == dh(2)) .and. (dh(1) == dh(3)))
       end select
    else
       mgt%dh(:,mgt%nlevels) = 1.0_dp_t
    end if

    do i = mgt%nlevels-1, 1, -1
       mgt%dh(:,i) = mgt%dh(:,i+1)*2.0_dp_t
    end do

    !   Set the face_type array to be BC_DIR or BC_NEU depending on domain_bc
    mgt%face_type = BC_INT
    mgt%bottom_singular = .true.
    do id = 1,mgt%dim
       lo_dom = lwb(pd,id)
       hi_dom = upb(pd,id)
       do i = 1,mgt%nboxes
          lo_grid =  lwb(get_box(mgt%ss(mgt%nlevels), i),id)
          if (lo_grid == lo_dom) mgt%face_type(i,id,1) = domain_bc(id,1)

          hi_grid = upb(get_box(mgt%ss(mgt%nlevels), i),id)
          if (hi_grid == hi_dom) mgt%face_type(i,id,2) = domain_bc(id,2)
       end do
    end do
    do id = 1,mgt%dim
       do i = 1,mgt%nboxes
          if (mgt%face_type(i,id,1) == BC_INT .or. &
               mgt%face_type(i,id,1) == BC_DIR) mgt%bottom_singular = .false.
          if (mgt%face_type(i,id,2) == BC_INT .or. &
               mgt%face_type(i,id,2) == BC_DIR) mgt%bottom_singular = .false.
       end do
    end do

    ! if ( mgt%bottom_solver == 0 .and. .not. present(bottom_max_iter) ) mgt%bottom_max_iter = 20

    if ( mgt%cycle_type == MG_WCycle ) mgt%gamma = 2

    ! if only the bottom solver is 'solving' make sure that its eps is
    ! in effect
    if ( mgt%nlevels == 1 ) then
       ba = get_boxarray(mgt%cc(1)%la)
       vol = boxarray_volume(ba)
       if (vol > 4**mgt%dim) &
         mgt%bottom_solver_eps = mgt%eps
    end if

    call destroy(bpt)

    ! If we're at a higher AMR level that coarsens to another mg_tower,
    !   or if you're here as part of a mac_applyop instead of a solve
    !   instead of coarsening within this one then don't bother
    !   creating the special bottom solver stuff
    if (mgt%nlevels == 1) mgt%bottom_solver = 1

    if (mgt%bottom_solver == 4) then

       allocate(mgt%bottom_mgt)

       ! Get the old/new coarse problem domain
       old_coarse_la = mgt%ss(1)%la
       coarse_pd = layout_get_pd(old_coarse_la)

       ! Get the new coarse boxarray and layout
       call box_build_2(bxs,coarse_pd%lo(1:mgt%dim),coarse_pd%hi(1:mgt%dim))
       call boxarray_build_bx(new_coarse_ba,bxs)

       ! This is how many levels could be built if we made just one grid
       n = max_mg_levels_bottom(new_coarse_ba,min_width)

       ! This is the user-imposed limit
       n = min(n,mgt%max_bottom_nlevel)

       if ( n .eq. 1) then
          call bl_error("DONT USE MG_BOTTOM_SOLVER == 4: BOTTOM GRID NOT PROPERLY DIVISIBLE")
       end if

       bottom_box_size = 2**n
       if (parallel_IOProcessor() .and. verbose .ge. 1) &
          print *,'F90mg: New levels in bottom_solver = 4: ',n

       call boxarray_maxsize(new_coarse_ba,bottom_box_size)
       call layout_build_ba(new_coarse_la,new_coarse_ba,coarse_pd, &
                            pmask=old_coarse_la%lap%pmask)

       call boxarray_destroy(new_coarse_ba)

       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          print *,'F90mg: Coarse problem domain for bottom_solver = 4: '
          call print(layout_get_pd(old_coarse_la))
          print *,'   ... Original boxes ',old_coarse_la%lap%nboxes
          print *,'   ... New      boxes ',new_coarse_la%lap%nboxes
       end if

       coarse_dx(:) = mgt%dh(:,1)

       call mg_tower_build(mgt%bottom_mgt, new_coarse_la, coarse_pd, &
                           domain_bc, &
                           dh = coarse_dx, &
                           ns = ns, &
                           nc = mgt%nc, &
                           ng = mgt%ng, &
                           smoother = smoother, &
                           nu1 = nu1, &
                           nu2 = nu2, &
                           gamma = gamma, &
                           cycle_type = cycle_type, &
                           omega = omega, &
                           bottom_solver = 1, &
                           bottom_max_iter = bottom_max_iter, &
                           bottom_solver_eps = bottom_solver_eps, &
                           max_iter = max_iter, &
                           max_nlevel = max_nlevel, &
                           min_width = min_width, &
                           eps = eps, &
                           abs_eps = abs_eps, &
                           verbose = verbose, &
                           cg_verbose = cg_verbose, &
                           nodal = nodal)
    end if

  end subroutine mg_tower_build

  subroutine mg_tower_print(mgt, unit, skip)

    use bl_IO_module

    type(mg_tower), intent(in) :: mgt
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: skip
    integer :: un, i, ii
    type(box) :: bb
    un = unit_stdout(unit)
    call unit_skip(un, skip)

    write(unit=un, fmt = '("F90MG settings...")')
    write(unit=un, fmt=*) 'nu1               = ', mgt%nu1
    call unit_skip(un, skip)
    write(unit=un, fmt=*) 'nu2               = ', mgt%nu2
    call unit_skip(un, skip)
    write(unit=un, fmt=*) 'nuf               = ', mgt%nuf
    call unit_skip(un, skip)
    write(unit=un, fmt=*) 'nub               = ', mgt%nub
    call unit_skip(un, skip)
    !   write(unit=un, fmt=*) 'ng                = ', mgt%ng
    !   call unit_skip(un, skip)
    !   write(unit=un, fmt=*) 'nc                = ', mgt%nc
    !   call unit_skip(un, skip)
    write(unit=un, fmt=*) 'min_width         = ', mgt%min_width
    call unit_skip(un, skip)
    write(unit=un, fmt=*) 'max_nlevel        = ', mgt%max_nlevel
    call unit_skip(un, skip)
    write(unit=un, fmt=*) 'max_iter          = ', mgt%max_iter
    call unit_skip(un, skip)
    !   write(unit=un, fmt=*) 'eps               = ', mgt%eps
    !   call unit_skip(un, skip)
    !   write(unit=un, fmt=*) 'abs_eps           = ', mgt%abs_eps
    !   call unit_skip(un, skip)
    !   write(unit=un, fmt=*) 'smoother          = ', mgt%smoother
    !   call unit_skip(un, skip)
    !   write(unit=un, fmt=*) 'gamma             = ', mgt%gamma
    !   call unit_skip(un, skip)
    !   write(unit=un, fmt=*) 'omega             = ', mgt%omega
    !   call unit_skip(un, skip)
    !   write(unit=un, fmt=*) 'cycle_type        = ', mgt%cycle_type
    !   call unit_skip(un, skip)
    write(unit=un, fmt=*) 'bottom_solver     = ', mgt%bottom_solver
    call unit_skip(un, skip)
    write(unit=un, fmt=*) 'bottom_solver_eps = ', mgt%bottom_solver_eps
    call unit_skip(un, skip)
    write(unit=un, fmt=*) 'bottom_max_iter   = ', mgt%bottom_max_iter
    call unit_skip(un, skip)

    if (mgt%verbose > 2) then
       print *,'F90MG: ',mgt%nlevels,' levels created for this solve'
       do i = mgt%nlevels,1,-1
          write(unit=un,fmt= '(" Level",i2)') i
          do ii = 1,nboxes(mgt%cc(i)%la)
             bb = get_box(mgt%cc(i)%la,ii)
             if (mgt%dim == 2) then
                write(unit=un,fmt= '("  [",i4,"]: (",i4,",",i4,") (",i4,",",i4,")",i4,i4 )') &
                     ii,bb%lo(1),bb%lo(2),bb%hi(1),bb%hi(2), &
                     bb%hi(1)-bb%lo(1)+1,bb%hi(2)-bb%lo(2)+1
             else if (mgt%dim == 3) then
                print *,'[',ii,']: ',bb%lo(1),bb%lo(2),bb%lo(3),bb%hi(1),bb%hi(2),bb%hi(3)
             end if
          end do
       end do
    else
       write(unit=un, fmt=*) 'F90MG: numlevels =  ',mgt%nlevels
       call unit_skip(un, skip)
    end if

  end subroutine mg_tower_print

  recursive subroutine mg_tower_destroy(mgt,destroy_la)

    type(mg_tower), intent(inout) :: mgt
    logical, intent(in), optional :: destroy_la

    logical :: ldestroy_la

    type(layout) :: la

    integer :: i

    ldestroy_la = .false.; if (present(destroy_la)) ldestroy_la = destroy_la

    la = mgt%cc(mgt%nlevels)%la

    do i = 1, mgt%nlevels
       call destroy(mgt%cc(i))
       call destroy(mgt%ff(i))
       call destroy(mgt%dd(i))
       call destroy(mgt%ss(i))
       call destroy(mgt%mm(i))
       if ( i /= mgt%nlevels ) then
          call destroy(mgt%uu(i))
       end if
    end do

    deallocate(mgt%cc, mgt%ff, mgt%dd, mgt%uu, mgt%mm, mgt%ss)
    deallocate(mgt%dh, mgt%pd)
    deallocate(mgt%tm)
    deallocate(mgt%face_type)
    deallocate(mgt%skewed)
    deallocate(mgt%skewed_not_set)
    if ( built_q(mgt%sparse_object) ) call destroy(mgt%sparse_object)
    if ( built_q(mgt%rh1)           ) call destroy(mgt%rh1)
    if ( built_q(mgt%uu1)           ) call destroy(mgt%uu1)
    if ( built_q(mgt%ss1)           ) call destroy(mgt%ss1)
    if ( built_q(mgt%mm1)           ) call destroy(mgt%mm1)

    if ( built_q(mgt%nodal_mask)    ) call destroy(mgt%nodal_mask)

    if (ldestroy_la) call layout_destroy(la)

    if ( associated(mgt%bottom_mgt) ) then
       call mg_tower_destroy(mgt%bottom_mgt, .true.)
       deallocate(mgt%bottom_mgt)
    end if

  end subroutine mg_tower_destroy

  function max_mg_levels(ba, min_size) result(r)

    type(boxarray), intent(in)           :: ba
    integer       , intent(in), optional :: min_size
    integer                              :: r
    integer, parameter :: rrr = 2
    type(box)          :: bx, bx1
    type(boxarray)     :: ba1
    real(kind=dp_t)    :: vol
    integer            :: i, rr, lmn, dm

    lmn = 1; if ( present(min_size) ) lmn = min_size
    r = 1
    rr = rrr
    dm = get_dim(ba)
    call copy(ba1,ba)
    do
       call boxarray_coarsen(ba1,rrr)
       vol = boxarray_volume(ba1)
       do i = 1, nboxes(ba)
          bx = get_box(ba,i)
          bx1 = coarsen(bx, rr)
          if ( any(extent(bx1) < lmn) ) then
             call destroy(ba1)
             return
          end if
          if ( bx /= refine(bx1, rr)  ) then
             call destroy(ba1)
             return
          end if

          ! We introduce the volume test to limit the case where we have a 
          !  single grid getting too small -- don't want to limit each grid in
          !  the case where there are many grids, so we need a test over the
          !  whole boxarray volume, not just the size of each grid.
          if ( vol <= 2**dm ) then
             call destroy(ba1)
             return
          end if

       end do
       rr = rr*rrr
       r  = r + 1
    end do

    call destroy(ba1)

  end function max_mg_levels

  function max_mg_levels_bottom(ba, min_size) result(r)

    type(boxarray), intent(in)           :: ba
    integer       , intent(in), optional :: min_size
    integer                              :: r
    integer, parameter :: rrr = 2
    type(box)          :: bx, bx1
    type(boxarray)     :: ba1
    real(kind=dp_t)    :: vol
    integer            :: i, rr, lmn, dm

    lmn = 1; if ( present(min_size) ) lmn = min_size
    r = 1
    rr = rrr
    dm = get_dim(ba)
    call copy(ba1,ba)
    do
       call boxarray_coarsen(ba1,rrr)
       vol = boxarray_volume(ba1)
       do i = 1, nboxes(ba)
          bx = get_box(ba,i)
          bx1 = coarsen(bx, rr)
          if ( any(extent(bx1) < lmn) .or. any(mod(extent(bx1),2) .eq. 1)) then
             call destroy(ba1)
             return
          end if
          if ( bx /= refine(bx1, rr)  ) then
             call destroy(ba1)
             return
          end if
       end do
       rr = rr*rrr
       r  = r + 1
    end do

    call destroy(ba1)

  end function max_mg_levels_bottom

  subroutine mg_tower_v_cycle(mgt, uu, rh)

    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in)    :: rh

    call mg_tower_cycle(mgt, MG_VCYCLE, mgt%nlevels, mgt%ss(mgt%nlevels), &
                        uu, rh, mgt%mm(mgt%nlevels), mgt%nu1, mgt%nu2, mgt%gamma)

  end subroutine mg_tower_v_cycle

  subroutine do_bottom_mgt(mgt, uu, rh)

    use bl_prof_module

    type( mg_tower), intent(inout) :: mgt
    type( multifab), intent(inout) :: uu
    type( multifab), intent(in   ) :: rh

    type(bl_prof_timer), save :: bpt

    type( multifab) :: bottom_uu
    type( multifab) :: bottom_rh
    integer         :: mglev
    logical         :: do_diag

    do_diag = .false.; if ( mgt%verbose >= 4 ) do_diag = .true.

    call build(bpt, "do_bottom_mgt")

    if (mgt%bottom_solver .ne. 4) &
       call bl_error("MG_TOWER_BOTTOM_SOLVE: must have bottom_solver == 4")

    mglev = mgt%bottom_mgt%nlevels

    call multifab_build(bottom_uu,mgt%bottom_mgt%ss(mglev)%la,1,uu%ng,uu%nodal)

    call setval(bottom_uu,0.d0,all=.true.)

    if (nodal_q(rh)) then
       call multifab_build(bottom_rh,mgt%bottom_mgt%ss(mglev)%la,1,1,rh%nodal)
       call setval(bottom_rh,ZERO,all=.true.)
    else
       call multifab_build(bottom_rh,mgt%bottom_mgt%ss(mglev)%la,1,0,rh%nodal)
    end if

    call multifab_copy_c(bottom_rh,1,rh,1,1,ng=0)

    call mg_tower_cycle(mgt%bottom_mgt, mgt%bottom_mgt%cycle_type, mglev, &
                        mgt%bottom_mgt%ss(mglev), bottom_uu, bottom_rh, &
                        mgt%bottom_mgt%mm(mglev), mgt%bottom_mgt%nu1, mgt%bottom_mgt%nu2, &
                        mgt%bottom_mgt%gamma)

    call multifab_copy_c(uu,1,bottom_uu,1,1,ng=0)
    call multifab_fill_boundary(uu)

    call destroy(bottom_uu)
    call destroy(bottom_rh)

    call destroy(bpt)

  end subroutine do_bottom_mgt

  subroutine mg_tower_bottom_solve(mgt, lev, ss, uu, rh, mm)

    use bl_prof_module
    use itsol_module, only: itsol_bicgstab_solve, itsol_cg_solve

    type( mg_tower), intent(inout) :: mgt
    type( multifab), intent(inout) :: uu
    type( multifab), intent(in) :: rh
    type( multifab), intent(in) :: ss
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: lev
    integer stat
    integer i
    type(bl_prof_timer), save :: bpt

    call build(bpt, "mgt_bottom_solve")

    stat = 0
    select case ( mgt%bottom_solver )
    case (0)
       do i = 1, mgt%nuf
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do
    case (1)
       if (nodal_q(rh)) then
          call itsol_bicgstab_solve(ss, uu, rh, mm, &
                                    mgt%bottom_solver_eps, mgt%bottom_max_iter, &
                                    mgt%cg_verbose, stat = stat, &
                                    singular_in = mgt%bottom_singular, &
                                    uniform_dh = mgt%uniform_dh,&
                                    nodal_mask = mgt%nodal_mask)
       else
          call itsol_bicgstab_solve(ss, uu, rh, mm, &
                                    mgt%bottom_solver_eps, mgt%bottom_max_iter, &
                                    mgt%cg_verbose, stat = stat, &
                                    singular_in = mgt%bottom_singular, &
                                    uniform_dh = mgt%uniform_dh)
       end if
       do i = 1, mgt%nub
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do
    case (2)
       if (nodal_q(rh)) then
          call itsol_cg_solve(ss, uu, rh, mm, &
                              mgt%bottom_solver_eps, mgt%bottom_max_iter, mgt%cg_verbose, &
                              stat = stat, singular_in = mgt%bottom_singular, &
                              uniform_dh = mgt%uniform_dh, nodal_mask=mgt%nodal_mask)
       else
          call itsol_cg_solve(ss, uu, rh, mm, &
                              mgt%bottom_solver_eps, mgt%bottom_max_iter, mgt%cg_verbose, &
                              stat = stat, singular_in = mgt%bottom_singular, &
                              uniform_dh = mgt%uniform_dh)
       end if
       do i = 1, mgt%nub
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do
    case (3)
       call copy(mgt%rh1, rh)
       if ( parallel_IOProcessor() ) then
          if (nodal_q(rh)) then
             call setval(mgt%uu1, zero)
             call sparse_nodal_solve(mgt%sparse_object, mgt%uu1, mgt%rh1, &
                                     mgt%bottom_solver_eps, mgt%bottom_max_iter, &
                                     mgt%verbose, stat)
          else
             call sparse_solve(mgt%sparse_object, mgt%uu1, mgt%rh1, &
                               mgt%bottom_solver_eps, mgt%bottom_max_iter, mgt%verbose, &
                               stat)
          end if
       end if
       call copy(uu, mgt%uu1)
       call parallel_bcast(stat, 1)

    case default
       call bl_error("MG_TOWER_BOTTOM_SOLVE: no such solver: ", mgt%bottom_solver)
    end select
    if ( stat /= 0 ) then
       if ( parallel_IOProcessor() ) then
          call bl_warn("BREAKDOWN in bottom_solver: trying smoother")
       end if
       do i = 1, 20
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do
    end if

    call destroy(bpt)

  end subroutine mg_tower_bottom_solve

  ! computes dd = ff - ss * uu
  subroutine mg_defect(ss, dd, ff, uu, mm, uniform_dh)

    use bl_prof_module
    use itsol_module, only: itsol_stencil_apply

    type(multifab), intent(in)    :: ff, ss
    type(multifab), intent(inout) :: dd, uu
    type(imultifab), intent(in)   :: mm
    logical, intent(in), optional :: uniform_dh
    type(bl_prof_timer), save     :: bpt

    call build(bpt, "mg_defect")
    call itsol_stencil_apply(ss, dd, uu, mm, uniform_dh)
    call saxpy(dd, ff, -1.0_dp_t, dd)
    call destroy(bpt)
  end subroutine mg_defect

  subroutine grid_res(mgt, ss, dd, ff, uu, mm, face_type, uniform_dh)

    use bl_prof_module
    use mg_defect_module

    type(multifab), intent(in)    :: ff, ss
    type(multifab), intent(inout) :: dd, uu
    type(imultifab), intent(in)   :: mm
    type(mg_tower), intent(inout) :: mgt
    integer, intent(in) :: face_type(:,:,:)
    logical, intent(in) :: uniform_dh
    integer :: i, n
    real(kind=dp_t), pointer :: dp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: dm,nodal_ng
    logical :: lcross
    type(bl_prof_timer), save :: bpt

    call build(bpt, "grid_res")

    nodal_ng = 0; if ( nodal_q(uu) ) nodal_ng = 1

    dm = mgt%dim
    lcross = (ncomp(ss) == (2*dm+1))

    call multifab_fill_boundary(uu, cross = lcross)

    do i = 1, mgt%nboxes
       if ( multifab_remote(dd, i) ) cycle
       dp => dataptr(dd, i)
       fp => dataptr(ff, i)
       up => dataptr(uu, i)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       do n = 1, mgt%nc
          select case(dm)
          case (1)
             call grid_laplace_1d(sp(:,1,1,:), dp(:,1,1,n), fp(:,1,1,n), up(:,1,1,n), &
                                  mgt%ng, nodal_ng)
          case (2)
             call grid_laplace_2d(sp(:,:,1,:), dp(:,:,1,n), fp(:,:,1,n), up(:,:,1,n), &
                                  mp(:,:,1,1), mgt%ng, face_type(i,:,:))
          case (3)
             call grid_laplace_3d(sp(:,:,:,:), dp(:,:,:,n), fp(:,:,:,n), up(:,:,:,n), &
                                  mp(:,:,:,1), mgt%ng, face_type(i,:,:), uniform_dh)
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine grid_res

  subroutine mg_tower_restriction(mgt, crse, fine, mm_fine, mm_crse)

    use bl_prof_module
    use mg_restriction_module, only: cc_restriction_1d, nodal_restriction_1d, &
         cc_restriction_2d, nodal_restriction_2d, cc_restriction_3d, nodal_restriction_3d

    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    type(imultifab), intent(in)   :: mm_fine
    type(imultifab), intent(in)   :: mm_crse
    type(mg_tower), intent(inout) :: mgt
    integer :: loc(mgt%dim), lof(mgt%dim)
    integer :: lom_fine(mgt%dim)
    integer :: lom_crse(mgt%dim)
    integer :: lo(mgt%dim), hi(mgt%dim)
    integer :: i, n, ir(mgt%dim)
    logical :: nodal_flag
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    integer        , pointer :: mp_fine(:,:,:,:)
    integer        , pointer :: mp_crse(:,:,:,:)

    integer :: mg_restriction_mode
    type(bl_prof_timer), save :: bpt

    call build(bpt, "mgt_restriction")

    ir = 2

    nodal_flag = nodal_q(crse)

    if ( nodal_flag ) then
       call multifab_fill_boundary(fine)
       call setval(crse, 0.0_dp_t, all=.true.)
       mg_restriction_mode = 1
    end if

    do i = 1, mgt%nboxes
       if ( multifab_remote(crse, i) ) cycle

       cp       => dataptr(crse, i)
       fp       => dataptr(fine, i)
       mp_fine  => dataptr(mm_fine, i)
       mp_crse  => dataptr(mm_crse, i)

       loc      = lwb(get_pbox(crse,i))
       lof      = lwb(get_pbox(fine,i))
       lom_fine = lwb(get_pbox(mm_fine,i))
       lom_crse = lwb(get_pbox(mm_crse,i))
       lo       = lwb(get_ibox(crse,i))
       hi       = upb(get_ibox(crse,i))

       do n = 1, mgt%nc
          select case ( mgt%dim)
          case (1)
             if ( .not. nodal_flag ) then
                call cc_restriction_1d(cp(:,1,1,n), loc, fp(:,1,1,n), lof, lo, hi, ir)
             else
                call nodal_restriction_1d(cp(:,1,1,n), loc, fp(:,1,1,n), lof, &
                                          mp_fine(:,1,1,1), lom_fine, &
                                          mp_crse(:,1,1,1), lom_crse, lo, hi, ir, .false., &
                                          mg_restriction_mode)
             end if
          case (2)
             if ( .not. nodal_flag ) then
                call cc_restriction_2d(cp(:,:,1,n), loc, fp(:,:,1,n), lof, lo, hi, ir)
             else
                call nodal_restriction_2d(cp(:,:,1,n), loc, fp(:,:,1,n), lof, &
                                          mp_fine(:,:,1,1), lom_fine, &
                                          mp_crse(:,:,1,1), lom_crse, lo, hi, ir, .false., &
                                          mg_restriction_mode)
             end if
          case (3)
             if ( .not. nodal_flag ) then
                call cc_restriction_3d(cp(:,:,:,n), loc, fp(:,:,:,n), lof, lo, hi, ir)
             else
                call nodal_restriction_3d(cp(:,:,:,n), loc, fp(:,:,:,n), lof, &
                                          mp_fine(:,:,:,1), lom_fine, &
                                          mp_crse(:,:,:,1), lom_crse, lo, hi, ir, .false., &
                                          mg_restriction_mode)
             end if
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine mg_tower_restriction

  subroutine mg_tower_smoother(mgt, lev, ss, uu, ff, mm)

    use bl_prof_module
    use mg_smoother_module, only: gs_line_solve_1d, gs_rb_smoother_1d, jac_smoother_2d, &
         jac_smoother_3d, gs_lex_smoother_2d, gs_lex_smoother_3d, nodal_smoother_1d, &
         nodal_smoother_2d, nodal_smoother_3d, gs_rb_smoother_2d,  gs_rb_smoother_3d, &
         minion_smoother_2d, minion_smoother_3d, nodal_line_solve_1d
    use itsol_module, only: itsol_bicgstab_solve, itsol_cg_solve

    integer        , intent(in   ) :: lev
    type( mg_tower), intent(inout) :: mgt
    type( multifab), intent(inout) :: uu
    type( multifab), intent(in   ) :: ff
    type( multifab), intent(in   ) :: ss
    type(imultifab), intent(in   ) :: mm

    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: i, k, n, nn, stat, npts
    integer :: lo(mgt%dim)
    type(bl_prof_timer), save :: bpt
    logical :: lcross, pmask(mgt%dim)
    real(kind=dp_t) :: local_eps

    pmask = layout_get_pmask(uu%la)

    call build(bpt, "mgt_smoother")

    lcross = (ncomp(ss) == 2*mgt%dim+1) 

    if (mgt%skewed_not_set(lev)) then 
       do i = 1, mgt%nboxes
          if ( remote(mm,i) ) cycle
          mp => dataptr(mm, i)
          mgt%skewed(lev,i) = skewed_q(mp)
       end do
       mgt%skewed_not_set(lev) = .false.
    end if

    if ( cell_centered_q(uu) ) then

       if (mgt%dim .eq. 1 .and. mgt%nboxes .eq. 1) then

          call multifab_fill_boundary(uu, cross = lcross)

          ! We do these line solves as a preconditioner
          do i = 1, mgt%nboxes
             if ( remote(ff, i) ) cycle
             up => dataptr(uu, i)
             fp => dataptr(ff, i)
             sp => dataptr(ss, i)
             mp => dataptr(mm, i)
             lo =  lwb(get_box(ss, i))
             do n = 1, mgt%nc
                call gs_line_solve_1d(sp(:,1,1,:), up(:,1,1,n), fp(:,1,1,n), &
                                      mp(:,1,1,1), lo, mgt%ng)
             end do
          end do

          call multifab_fill_boundary(uu, cross = lcross)

          !         local_eps = 0.001 
          !         call itsol_cg_solve(&
          !              ss, uu, ff, mm, &
          !              local_eps, 1050, mgt%cg_verbose, stat = stat, &
          !              local_eps, 1050, 1, stat = stat, &
          !              singular_in = mgt%bottom_singular, uniform_dh = mgt%uniform_dh)

          local_eps = 0.001 
          npts = multifab_volume(ff)
          call itsol_bicgstab_solve(ss, uu, ff, mm, &
                                    local_eps, npts, mgt%cg_verbose, stat = stat, &
                                    singular_in = mgt%bottom_singular, &
                                    uniform_dh = mgt%uniform_dh)

          call multifab_fill_boundary(uu, cross = lcross)

          if ( stat /= 0 ) then
             if ( parallel_IOProcessor() ) call bl_error("BREAKDOWN in 1d CG solve")
          end if

       else

          select case ( mgt%smoother )

          case ( MG_SMOOTHER_GS_RB )

             do nn = 0, 1
                call multifab_fill_boundary(uu, cross = lcross)

                do i = 1, mgt%nboxes
                   if ( remote(ff, i) ) cycle
                   up => dataptr(uu, i)
                   fp => dataptr(ff, i)
                   sp => dataptr(ss, i)
                   mp => dataptr(mm, i)
                   lo =  lwb(get_box(ss, i))
                   do n = 1, mgt%nc
                      select case ( mgt%dim)
                      case (1)
                         call gs_rb_smoother_1d(mgt%omega, sp(:,1,1,:), up(:,1,1,n), &
                                                fp(:,1,1,n), mp(:,1,1,1), lo, mgt%ng, nn, &
                                                mgt%skewed(lev,i))
                      case (2)
                         call gs_rb_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,n), &
                                                fp(:,:,1,n), mp(:,:,1,1), lo, mgt%ng, nn, &
                                                mgt%skewed(lev,i))
                      case (3)
                         call gs_rb_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,n), &
                                                fp(:,:,:,n), mp(:,:,:,1), lo, mgt%ng, nn, &
                                                mgt%skewed(lev,i))
                      end select
                   end do
                end do

             end do

          case ( MG_SMOOTHER_EFF_RB )

             call multifab_fill_boundary(uu, cross = lcross)
             do i = 1, mgt%nboxes
                if ( remote(ff, i) ) cycle
                up => dataptr(uu, i)
                fp => dataptr(ff, i)
                sp => dataptr(ss, i)
                mp => dataptr(mm, i)
                lo =  lwb(get_box(ss, i))
                do n = 1, mgt%nc
                   select case ( mgt%dim)
                   case (1)
                      call gs_rb_smoother_1d(mgt%omega, sp(:,1,1,:), up(:,1,1,n), &
                                             fp(:,1,1,n), mp(:,1,1,1), lo, mgt%ng, 0, &
                                             mgt%skewed(lev,i))
                      call gs_rb_smoother_1d(mgt%omega, sp(:,1,1,:), up(:,1,1,n), &
                                             fp(:,1,1,n), mp(:,1,1,1), lo, mgt%ng, 1, &
                                             mgt%skewed(lev,i))
                   case (2)
                      call gs_rb_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,n), &
                                             fp(:,:,1,n), mp(:,:,1,1), lo, mgt%ng, 0, &
                                             mgt%skewed(lev,i))
                      call gs_rb_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,n), &
                                             fp(:,:,1,n), mp(:,:,1,1), lo, mgt%ng, 1, &
                                             mgt%skewed(lev,i))
                   case (3)
                      call gs_rb_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,n), &
                                             fp(:,:,:,n), mp(:,:,:,1), lo, mgt%ng, 0, &
                                             mgt%skewed(lev,i))
                      call gs_rb_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,n), &
                                             fp(:,:,:,n), mp(:,:,:,1), lo, mgt%ng, 1, &
                                             mgt%skewed(lev,i))
                   end select
                end do
             end do

          case ( MG_SMOOTHER_MINION_CROSS )

             call multifab_fill_boundary(uu, cross = lcross)
             do i = 1, mgt%nboxes
                if ( remote(ff, i) ) cycle
                up => dataptr(uu, i)
                fp => dataptr(ff, i)
                sp => dataptr(ss, i)
                mp => dataptr(mm, i)
                lo =  lwb(get_box(ss, i))
                do n = 1, mgt%nc
                   select case ( mgt%dim)
                   case (2)
                      call minion_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,1), &
                                              fp(:,:,1,1), lo, mgt%ng, .true.)
                   case (3)
                      call minion_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,1), &
                                              fp(:,:,:,1), lo, mgt%ng, .true.)
                   end select
                end do
             end do

          case ( MG_SMOOTHER_MINION_FULL )

             call multifab_fill_boundary(uu, cross = lcross)
             do i = 1, mgt%nboxes
                if ( remote(ff, i) ) cycle
                up => dataptr(uu, i)
                fp => dataptr(ff, i)
                sp => dataptr(ss, i)
                mp => dataptr(mm, i)
                lo =  lwb(get_box(ss, i))
                do n = 1, mgt%nc
                   select case ( mgt%dim)
                   case (2)
                      call minion_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,1), &
                                              fp(:,:,1,1), lo, mgt%ng, .false.)
                   case (3)
                      call minion_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,1), &
                                              fp(:,:,:,1), lo, mgt%ng, .false.)
                   end select
                end do
             end do

          case ( MG_SMOOTHER_JACOBI )
             call multifab_fill_boundary(uu, cross = lcross)
             do i = 1, mgt%nboxes
                if ( remote(ff, i) ) cycle
                up => dataptr(uu, i)
                fp => dataptr(ff, i)
                sp => dataptr(ss, i)
                mp => dataptr(mm, i)
                lo =  lwb(get_box(ss, i))
                do n = 1, mgt%nc
                   select case ( mgt%dim)
                   case (2)
                      call jac_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,n), fp(:,:,1,n), &
                                           mp(:,:,1,1), mgt%ng)
                   case (3)
                      call jac_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,n), fp(:,:,:,n), &
                                           mp(:,:,:,1), mgt%ng)
                   end select
                end do
             end do
          case ( MG_SMOOTHER_GS_LEX )
             call multifab_fill_boundary(uu, cross = lcross)
             do i = 1, mgt%nboxes
                if ( remote(ff, i) ) cycle
                up => dataptr(uu, i)
                fp => dataptr(ff, i)
                sp => dataptr(ss, i)
                mp => dataptr(mm, i)
                lo =  lwb(get_box(ss, i))
                do n = 1, mgt%nc
                   select case ( mgt%dim)
                   case (2)
                      call gs_lex_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,n), &
                                              fp(:,:,1,n), mgt%ng)
                   case (3)
                      call gs_lex_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,n), &
                                              fp(:,:,:,n), mgt%ng)
                   end select
                end do
             end do
          case default
             call bl_error("MG_TOWER_SMOOTHER: no such smoother")
          end select

       end if ! if (mgt%dim > 1)

    else 
       if (lcross) then

        if (ff%dim == 1) then

             call multifab_fill_boundary(uu, cross = lcross)
             do i = 1, mgt%nboxes
                if ( remote(ff, i) ) cycle
                up => dataptr(uu, i)
                fp => dataptr(ff, i)
                sp => dataptr(ss, i)
                mp => dataptr(mm, i)
                lo =  lwb(get_box(ss, i))
                do n = 1, mgt%nc
                   select case ( mgt%dim)
                   case (1)
                      call nodal_line_solve_1d(sp(:,1,1,:), up(:,1,1,n), fp(:,1,1,n), &
                                               mp(:,1,1,1), lo, mgt%ng)
                   end select
                end do
             end do

        else

          ! k is the red-black parameter
          do k = 0, 1
             call multifab_fill_boundary(uu, cross = lcross)
             do i = 1, mgt%nboxes
                if ( remote(ff, i) ) cycle
                up => dataptr(uu, i)
                fp => dataptr(ff, i)
                sp => dataptr(ss, i)
                mp => dataptr(mm, i)
                lo =  lwb(get_box(ss, i))
                do n = 1, mgt%nc
                   select case ( mgt%dim)
                   case (1)
!                     call nodal_smoother_1d(mgt%omega, sp(:,1,1,:), up(:,1,1,n), fp(:,1,1,n), &
!                                            mp(:,1,1,1), lo, mgt%ng, k)
                      if (k.eq.0) &
                      call nodal_line_solve_1d(sp(:,1,1,:), up(:,1,1,n), fp(:,1,1,n), &
                                               mp(:,1,1,1), lo, mgt%ng)
                   case (2)
                      call nodal_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,n), &
                                             fp(:,:,1,n), mp(:,:,1,1), lo, mgt%ng, pmask, k)
                   case (3)
                      call nodal_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,n), &
                                             fp(:,:,:,n), mp(:,:,:,1), lo, mgt%ng, &
                                             mgt%uniform_dh, pmask, k)
                   end select
                end do
             end do
          end do
        end if
       else
          call multifab_fill_boundary(uu, cross = lcross)
          ! This value of k isn't used
          k = 0
          do i = 1, mgt%nboxes
             if ( remote(ff, i) ) cycle
             up => dataptr(uu, i)
             fp => dataptr(ff, i)
             sp => dataptr(ss, i)
             mp => dataptr(mm, i)
             lo =  lwb(get_box(ss, i))
             do n = 1, mgt%nc
                select case ( mgt%dim)
                case (1)
!                  call nodal_smoother_1d(mgt%omega, sp(:,1,1,:), up(:,1,1,n), fp(:,1,1,n), &
!                                         mp(:,1,1,1), lo, mgt%ng, k)
                   call nodal_line_solve_1d(sp(:,1,1,:), up(:,1,1,n), fp(:,1,1,n), &
                                            mp(:,1,1,1), lo, mgt%ng)
                case (2)
                   call nodal_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,n), fp(:,:,1,n), &
                                          mp(:,:,1,1), lo, mgt%ng, pmask, k)
                case (3)
                   call nodal_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,n), fp(:,:,:,n), &
                                          mp(:,:,:,1), lo, mgt%ng, mgt%uniform_dh, pmask, k)
                end select
             end do
          end do
       end if

       call multifab_internal_sync(uu)
    end if

    call destroy(bpt)

  end subroutine mg_tower_smoother

  subroutine mg_jacobi_smoother(mgt, lev, ss, uu, ff, mm)

    use bl_prof_module
    use mg_smoother_module, only: jac_smoother_1d, jac_smoother_2d, jac_smoother_3d

    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: ff
    type(multifab), intent(in) :: ss
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: lev
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: i, n, iter
    type(bl_prof_timer), save :: bpt
    logical :: lcross

    call build(bpt, "mgt_jacobi_smoother")

    lcross = ((ncomp(ss) == 5) .or. (ncomp(ss) == 7))

    if (mgt%skewed_not_set(lev)) then 
       do i = 1, mgt%nboxes
          if ( remote(mm,i) ) cycle
          mp => dataptr(mm, i)
          mgt%skewed(lev,i) = skewed_q(mp)
       end do
       mgt%skewed_not_set(lev) = .false.
    end if

    if ( .not. cell_centered_q(uu) ) then
       call bl_error('MG_JACOBI_SMOOTHER ONLY DESIGNED FOR CELL-CENTERED RIGHT NOW')
    end if

    do iter = 1, mgt%nu1
       call multifab_fill_boundary(uu, cross = lcross)
       do i = 1, mgt%nboxes
          if ( remote(ff, i) ) cycle
          up => dataptr(uu, i)
          fp => dataptr(ff, i)
          sp => dataptr(ss, i)
          mp => dataptr(mm, i)
          do n = 1, mgt%nc
             select case ( mgt%dim)
             case (1)
                call jac_smoother_1d(mgt%omega, sp(:,1,1,:), up(:,1,1,n), fp(:,1,1,n), &
                                     mp(:,1,1,1), mgt%ng)
             case (2)
                call jac_smoother_2d(mgt%omega, sp(:,:,1,:), up(:,:,1,n), fp(:,:,1,n), &
                                     mp(:,:,1,1), mgt%ng)
             case (3)
                call jac_smoother_3d(mgt%omega, sp(:,:,:,:), up(:,:,:,n), fp(:,:,:,n), &
                                     mp(:,:,:,1), mgt%ng)
             end select
          end do
       end do
    end do

    call destroy(bpt)

  end subroutine mg_jacobi_smoother

  subroutine mg_tower_prolongation(mgt, lev, uu, uu1)

    use bl_prof_module
    use mg_prolongation_module

    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(inout) :: uu, uu1
    integer, intent(in) :: lev
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    type(box) :: nbox, nbox1
    integer :: i, n, ir(mgt%dim)
    logical :: nodal_flag
    type(bl_prof_timer), save :: bpt

    call build(bpt, "mgt_prolongation")

    ir = 2

    nodal_flag = nodal_q(uu)

    if ( .not.nodal_flag ) then
       do i = 1, mgt%nboxes
          if ( remote(mgt%ff(lev), i) ) cycle
          fp => dataptr(uu,  i, get_box(uu,i))
          cp => dataptr(uu1, i, get_box(uu1,i))
          do n = 1, mgt%nc
             select case ( mgt%dim)
             case (1)
                call pc_c_prolongation(fp(:,1,1,n), cp(:,1,1,n), ir)
             case (2)
                call pc_c_prolongation(fp(:,:,1,n), cp(:,:,1,n), ir)
             case (3)
                call pc_c_prolongation(fp(:,:,:,n), cp(:,:,:,n), ir)
             end select
          end do
       end do
    else
       do i = 1, mgt%nboxes
          if ( remote(mgt%ff(lev), i) ) cycle
          nbox  = box_grow_n_f(get_box(uu,i),1,1)
          nbox1 = box_grow_n_f(get_box(uu1,i),1,1)
          fp => dataptr(uu,  i, nbox )
          cp => dataptr(uu1, i, nbox1)
          do n = 1, mgt%nc
             select case ( mgt%dim)
             case (1)
                call nodal_prolongation(fp(:,1,1,n), cp(:,1,1,n), ir)
             case (2)
                call nodal_prolongation(fp(:,:,1,n), cp(:,:,1,n), ir)
             case (3)
                call nodal_prolongation(fp(:,:,:,n), cp(:,:,:,n), ir)
             end select
          end do
       end do
    endif

    call destroy(bpt)

  end subroutine mg_tower_prolongation

  function mg_tower_converged(mgt, dd, uu, Ynorm) result(r)

    use itsol_module

    logical :: r
    type(mg_tower), intent(inout) :: mgt
    real(dp_t), intent(in) :: Ynorm
    type(multifab), intent(in) :: dd
    type(multifab), intent(in) :: uu
    r = itsol_converged(dd, uu, Ynorm, mgt%eps, mgt%abs_eps)
  end function mg_tower_converged

  recursive subroutine mg_tower_cycle(mgt, cyc, lev, ss, uu, rh, mm, nu1, nu2, gamma, &
                                      bottom_level)

    use bl_prof_module

    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(in) :: rh
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: ss
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: lev
    integer, intent(in) :: nu1, nu2
    integer, intent(inout) :: gamma
    integer, intent(in) :: cyc
    integer, intent(in), optional :: bottom_level
    integer :: i
    logical :: do_diag
    real(dp_t) :: nrm
    integer :: lbl
    logical :: nodal_flag
    type(bl_prof_timer), save :: bpt

    call build(bpt, "mgt_cycle")

    lbl = 1; if ( present(bottom_level) ) lbl = bottom_level

    do_diag = .false.; if ( mgt%verbose >= 4 ) do_diag = .true.

    nodal_flag = nodal_q(ss)

    call timer_start(mgt%tm(lev))

    if ( parallel_IOProcessor() .and. do_diag) &
       write(6,1000) lev

    if ( rh%dim == 1 ) then

       if (do_diag) then
          nrm = norm_inf(rh)
          if ( parallel_IOProcessor() ) &
             print *,'  DN: Norm before smooth         ',nrm
       end if

       call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)

       ! compute mgt%cc(lev) = ss * uu - rh
       call mg_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%uniform_dh)

       if (do_diag) then
          nrm = norm_inf(mgt%cc(lev))
          if ( parallel_IOProcessor() ) &
              print *,'  DN: Norm after  smooth          ',nrm
       end if

    else if ( lev == lbl ) then

       if (do_diag) then
          ! compute mgt%cc(lev) = ss * uu - rh
          call mg_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%uniform_dh)
          nrm = norm_inf(mgt%cc(lev))
          if ( parallel_IOProcessor() ) &
             print *,'  DN: Norm before bottom         ',nrm
       end if

       if (associated(mgt%bottom_mgt)) then
          call do_bottom_mgt(mgt, uu, rh)
       else
          call mg_tower_bottom_solve(mgt, lev, ss, uu, rh, mm)
       end if

       if ( cyc == MG_FCycle ) gamma = 1
       if (do_diag) then
          ! compute mgt%cc(lev) = ss * uu - rh
          call mg_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%uniform_dh)
          nrm = norm_inf(mgt%cc(lev))
          if ( parallel_IOProcessor() ) &
             print *,'  UP: Norm after  bottom         ',nrm
       end if

    else 

       if (do_diag) then
          nrm = norm_inf(rh)
          if ( parallel_IOProcessor() ) &
             print *,'  DN: Norm before smooth         ',nrm
       end if

       do i = 1, nu1
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do

       ! compute mgt%cc(lev) = ss * uu - rh
       call mg_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%uniform_dh)

       if (do_diag) then
          nrm = norm_inf(mgt%cc(lev))
          if ( parallel_IOProcessor() ) &
              print *,'  DN: Norm after  smooth         ',nrm
       end if

       call mg_tower_restriction(mgt, mgt%dd(lev-1), mgt%cc(lev), &
                                 mgt%mm(lev),mgt%mm(lev-1))
       ! HACK 
       if (nodal_flag) then 
          if (rh%dim .eq. 3) then
             call multifab_mult_mult_s(mgt%dd(lev-1),0.125_dp_t,mgt%dd(lev-1)%ng)
          else if (rh%dim .eq. 2) then
             call multifab_mult_mult_s(mgt%dd(lev-1),0.25_dp_t,mgt%dd(lev-1)%ng)
          end if
       end if
       call setval(mgt%uu(lev-1), zero, all = .TRUE.)

       do i = gamma, 1, -1
          call mg_tower_cycle(mgt, cyc, lev-1, mgt%ss(lev-1), mgt%uu(lev-1), &
                              mgt%dd(lev-1), mgt%mm(lev-1), nu1, nu2, gamma, bottom_level)
       end do
       ! uu  += cc, done, by convention, using the prolongation routine.

       call mg_tower_prolongation(mgt, lev, uu, mgt%uu(lev-1))

       if ( parallel_IOProcessor() .and. do_diag) &
          write(6,1000) lev

       if (do_diag) then
          ! compute mgt%cc(lev) = ss * uu - rh
          call mg_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%uniform_dh)
          nrm = norm_inf(mgt%cc(lev))
          if ( parallel_IOProcessor() ) then
             print *,'  UP: Norm after  interp         ',nrm
          end if
       end if

       do i = 1, nu2
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do

       if (do_diag) then
          ! compute mgt%cc(lev) = ss * uu - rh
          call mg_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%uniform_dh)
          nrm = norm_inf(mgt%cc(lev))
          if ( parallel_IOProcessor() ) &
             print *,'  UP: Norm after  smooth         ',nrm
       end if

       ! if at top of tower and doing an FCycle reset gamma
       if ( lev == mgt%nlevels .AND. cyc == MG_FCycle ) gamma = 2
    end if

    call timer_stop(mgt%tm(lev))

    call destroy(bpt)

1000 format('AT LEVEL ',i2)

  end subroutine mg_tower_cycle

  subroutine mini_cycle(mgt, lev, ss, uu, rh, mm, nu1, nu2)

    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(in) :: rh
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: ss
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: lev
    integer, intent(in) :: nu1, nu2
    integer :: i

    call timer_start(mgt%tm(lev))

    if ( lev == 1 ) then

       do i = 1, nu1
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do

    else 

       ! compute mgt%cc(lev) = ss * uu - rh
       call mg_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%uniform_dh)

       call mg_tower_restriction(mgt, mgt%dd(lev-1), mgt%cc(lev), &
                                 mgt%mm(lev),mgt%mm(lev-1))

       ! HACK 
       if (nodal_q(mgt%dd(lev-1))) then
          if (ss%dim .eq. 3) then
             call multifab_mult_mult_s(mgt%dd(lev-1),0.125_dp_t,mgt%dd(lev-1)%ng)
          else if (ss%dim .eq. 2) then
             call multifab_mult_mult_s(mgt%dd(lev-1),0.25_dp_t,mgt%dd(lev-1)%ng)
          end if
       end if

       call setval(mgt%uu(lev-1), zero, all = .TRUE.)

       call mg_tower_bottom_solve(mgt, lev-1, mgt%ss(lev-1), mgt%uu(lev-1), &
                                  mgt%dd(lev-1), mgt%mm(lev-1))

       call mg_tower_prolongation(mgt, lev, uu, mgt%uu(lev-1))

       do i = 1, nu2
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do

    end if

    call timer_stop(mgt%tm(lev))

  end subroutine mini_cycle

  subroutine mg_tower_solve(mgt, uu, rh, qq, num_iter, defect_history, defect_dirname, stat)

    use fabio_module
    use bl_prof_module

    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(inout) :: uu, rh
    integer, intent(out), optional :: stat
    real(dp_t), intent(out), optional :: qq(0:)
    integer, intent(out), optional :: num_iter
    character(len=*), intent(in), optional :: defect_dirname
    logical, intent(in), optional :: defect_history
    integer :: gamma
    integer :: it, cyc
    real(dp_t) :: ynorm, Anorm, nrm(3)
    integer :: n_qq, i_qq
    logical :: ldef
    type(bl_prof_timer), save :: bpt
    character(len=128) :: defbase

    call build(bpt, "mg_tower_solve")

    ldef = .false.; if ( present(defect_history) ) ldef = defect_history
    if ( ldef ) then
       if ( .not. present(defect_dirname) ) then
          call bl_error("MG_TOWER_SOLVE: defect_history but no defect_dirname")
       end if
    end if

    if ( present(stat) ) stat = 0
    if ( mgt%cycle_type == MG_FCycle ) then
       gamma = 2
    else
       gamma = mgt%gamma
    end if

    n_qq = 0
    i_qq = 0
    if ( present(qq) ) then
       n_qq = size(qq)
    end if
    cyc   = mgt%cycle_type
    Anorm = stencil_norm(mgt%ss(mgt%nlevels))
    ynorm = norm_inf(rh)

    ! compute mgt%dd(mgt%nlevels) = mgt%ss(mgt%nlevels) * uu - rh
    call mg_defect(mgt%ss(mgt%nlevels), &
                   mgt%dd(mgt%nlevels), rh, uu, mgt%mm(mgt%nlevels), mgt%uniform_dh)
    if ( i_qq < n_qq ) then
       qq(i_qq) = norm_l2(mgt%dd(mgt%nlevels))
       i_qq = i_qq + 1
    end if
    if ( ldef ) then
       write(unit = defbase, fmt='("def",I3.3)') 0
       call fabio_write(mgt%dd(mgt%nlevels), defect_dirname, defbase)
    end if

    if ( mgt%verbose > 0 ) then
       nrm(3) = norm_inf(mgt%dd(mgt%nlevels))
       nrm(1) = norm_inf(uu)
       nrm(2) = norm_inf(rh)
       if ( parallel_IOProcessor() ) then
          write(unit=*, &
               fmt='(i3,": Unorm= ",g15.8,", Rnorm= ",g15.8,", Ninf(defect) = ",g15.8, ", Anorm=",g15.8)') &
               0, nrm, Anorm
       end if
    end if
    if ( mg_tower_converged(mgt, mgt%dd(mgt%nlevels), uu, Ynorm) ) then
       if ( present(stat) ) stat = 0
       if ( mgt%verbose > 0 .AND. parallel_IOProcessor() ) then
          write(unit=*, fmt='("F90mg: MG finished at on input")') 
       end  if
       return
    end if
    do it = 1, mgt%max_iter
       call mg_tower_cycle(mgt, cyc, mgt%nlevels, mgt%ss(mgt%nlevels), &
                           uu, rh, mgt%mm(mgt%nlevels), mgt%nu1, mgt%nu2, gamma)
       ! compute mgt%cc(lev) = ss * uu - rh
       call mg_defect(mgt%ss(mgt%nlevels), &
                      mgt%dd(mgt%nlevels), rh, uu, mgt%mm(mgt%nlevels), mgt%uniform_dh)
       if ( mgt%verbose > 0 ) then
          nrm(1) = norm_inf(mgt%dd(mgt%nlevels))
          if ( parallel_IOProcessor() ) then
             write(unit=*, fmt='(i3,": Ninf(defect) = ",g15.8)') it, nrm(1)
          end if
       end if
       if ( i_qq < n_qq ) then
          qq(i_qq) = norm_l2(mgt%dd(mgt%nlevels))
          i_qq = i_qq + 1
       end if
       if ( ldef ) then
          write(unit = defbase,fmt='("def",I3.3)') it
          call fabio_write(mgt%dd(mgt%nlevels), defect_dirname, defbase)
       end if
       if ( mg_tower_converged(mgt, mgt%dd(mgt%nlevels), uu, Ynorm) ) exit
    end do
    if ( mgt%verbose > 0 .AND. parallel_IOProcessor() ) then
       write(unit=*, fmt='("F90mg: MG finished at ", i3, " iterations")') it
    end  if

    if ( it > mgt%max_iter ) then
       if ( present(stat) ) then
          stat = 1
       else
          call bl_error("MG_TOWER_SOLVE: failed to converge in max_iter iterations")
       end if
    end if

    if ( present(num_iter) ) num_iter = it

    call destroy(bpt)

  end subroutine mg_tower_solve

end module mg_module
