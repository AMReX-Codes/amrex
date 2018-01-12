module mg_module

  use multifab_module
  use cc_stencil_module
  use mg_tower_module
  use stencil_defect_module
  use stencil_types_module
  use bc_functions_module

  implicit none

  integer, save :: n_mg_iters = -1

  interface destroy
     module procedure mg_tower_destroy
  end interface

  private :: get_bottom_box_size, impose_physbc_cc_2d, impose_physbc_cc_3d

contains

  integer function amrex_f90mg_get_niters () bind(c,name='amrex_f90mg_get_niters')
    amrex_f90mg_get_niters = n_mg_iters
  end function amrex_f90mg_get_niters

  recursive subroutine mg_tower_build(mgt, la, pd, domain_bc, stencil_type, &
                                      nu1, nu2, nuf, nub, cycle_type, &
                                      smoother, &
                                      dh, &
                                      nc, ng, &
                                      max_nlevel, max_bottom_nlevel, min_width, &
                                      max_iter, abort_on_max_iter, eps, abs_eps, &
                                      bottom_solver, bottom_max_iter, bottom_solver_eps, &
                                      max_L0_growth, verbose, cg_verbose, nodal, use_hypre, &
                                      is_singular, ok_to_fix_singular, &
                                      the_bottom_comm, fancy_bottom_type, use_lininterp, ptype)
    use bl_IO_module
    use bl_prof_module

    type(mg_tower), intent(inout) :: mgt
    type(layout), intent(in   ) :: la
    type(box), intent(in) :: pd
    integer, intent(in) :: domain_bc(:,:)
    integer, intent(in) :: stencil_type

    integer, intent(in), optional :: nu1
    integer, intent(in), optional :: nu2
    integer, intent(in), optional :: nuf
    integer, intent(in), optional :: nub
    integer, intent(in), optional :: cycle_type
    integer, intent(in), optional :: smoother
    real(kind=dp_t), intent(in), optional :: dh(:)
    integer, intent(in), optional :: nc
    integer, intent(in), optional :: ng
    integer, intent(in), optional :: max_nlevel
    integer, intent(in), optional :: max_bottom_nlevel
    integer, intent(in), optional :: min_width
    integer, intent(in), optional :: max_iter
    logical, intent(in), optional :: abort_on_max_iter
    real(kind=dp_t), intent(in), optional :: eps
    real(kind=dp_t), intent(in), optional :: abs_eps
    integer, intent(in), optional :: bottom_solver
    integer, intent(in), optional :: bottom_max_iter
    real(kind=dp_t), intent(in), optional :: bottom_solver_eps
    real(kind=dp_t), intent(in), optional :: max_L0_growth
    integer, intent(in), optional :: verbose
    integer, intent(in), optional :: cg_verbose
    logical, intent(in), optional :: nodal(:)
    integer, intent(in), optional :: use_hypre
    logical, intent(in), optional :: is_singular
    logical, intent(in), optional :: ok_to_fix_singular
    integer, intent(in), optional :: the_bottom_comm
    integer, intent(in), optional :: fancy_bottom_type
    logical, intent(in), optional :: use_lininterp
    integer, intent(in), optional :: ptype

    integer         :: n, i, id, j, lo_grid, hi_grid, lo_dom, hi_dom, ng_for_res
    type(layout)    :: la1, la2
    type(box)       :: bounding_box
    type(boxarray)  :: ba
    logical         :: nodal_flag, eq_vol
    real(kind=dp_t) :: dvol, dvol_pd
    integer(kind=ll_t) :: vol, vol_pd
    type(bl_prof_timer), save :: bpt

    ! These are added to help build the bottom_mgt
    type(  layout)      :: old_coarse_la, new_coarse_la
    type(     box)      :: coarse_pd,bxs
    type(boxarray)      :: new_coarse_ba
    integer             :: bottom_box_size, success, communicator
    real(kind=dp_t)     :: coarse_dx(pd%dim)
    real(dp_t), pointer :: p(:,:,:,:)

    call bl_proffortfuncstart("mg_tower_build")
    call build(bpt, "mgt_build")

    ! Need to set this here because so many things below depend on it.
    mgt%dim = get_dim(la)

    ! Paste in optional arguments
    if ( present(ng)                ) mgt%ng                = ng
    if ( present(nc)                ) mgt%nc                = nc
    if ( present(max_nlevel)        ) mgt%max_nlevel        = max_nlevel
    if ( present(max_bottom_nlevel) ) mgt%max_bottom_nlevel = max_bottom_nlevel
    if ( present(max_iter)          ) mgt%max_iter          = max_iter
    if ( present(abort_on_max_iter) ) mgt%abort_on_max_iter = abort_on_max_iter
    if ( present(eps)               ) mgt%eps               = eps
    if ( present(abs_eps)           ) mgt%abs_eps           = abs_eps
    if ( present(smoother)          ) mgt%smoother          = smoother
    if ( present(nu1)               ) mgt%nu1               = nu1
    if ( present(nu2)               ) mgt%nu2               = nu2
    if ( present(nuf)               ) mgt%nuf               = nuf
    if ( present(nub)               ) mgt%nub               = nub
    if ( present(cycle_type)        ) mgt%cycle_type        = cycle_type
    if ( present(fancy_bottom_type) ) mgt%fancy_bottom_type = fancy_bottom_type
    if ( present(bottom_solver)     ) mgt%bottom_solver     = bottom_solver
    if ( present(bottom_solver_eps) ) mgt%bottom_solver_eps = bottom_solver_eps
    if ( present(bottom_max_iter)   ) mgt%bottom_max_iter   = bottom_max_iter
    if ( present(min_width)         ) mgt%min_width         = min_width
    if ( present(verbose)           ) mgt%verbose           = verbose
    if ( present(cg_verbose)        ) mgt%cg_verbose        = cg_verbose
    if ( present(use_hypre)         ) mgt%use_hypre         = use_hypre 
    if ( present(ok_to_fix_singular)) mgt%ok_to_fix_singular= ok_to_fix_singular
    if ( present(max_L0_growth)     ) mgt%max_L0_growth     = max_L0_growth 
    if ( present(use_lininterp)     ) mgt%use_lininterp     = use_lininterp
    if ( present(ptype)             ) mgt%ptype             = ptype

    if ( present(the_bottom_comm) ) then
       allocate(mgt%bottom_comm)
       mgt%bottom_comm = the_bottom_comm
    end if

    nodal_flag = .false.
    if ( present(nodal) ) then
       if ( all(nodal) ) then
          nodal_flag = .true.
       else if ( any(nodal) ) then
          call bl_error("mixed nodal/not nodal forbidden")
       end if
    end if

    mgt%stencil_type = stencil_type

    ! "lcross" is used in the call to multifab_fill_boundary -- when true it means
    !   that the corner cells do *not* need to be filled
    if (mgt%stencil_type .eq. CC_CROSS_STENCIL .or. &
        mgt%stencil_type .eq. HO_CROSS_STENCIL .or. &
        mgt%stencil_type .eq. ND_CROSS_STENCIL) then

        mgt%lcross = .true.

    else if (mgt%stencil_type .eq. ND_DENSE_STENCIL .or. &
             mgt%stencil_type .eq. ND_VATER_STENCIL .or. &
             mgt%stencil_type .eq. HO_DENSE_STENCIL) then

        mgt%lcross = .false.

    else 

       print *,'The stencil_type must be specified '
       print *, mgt%stencil_type,' is not a valid option'
       call bl_error("")

    end if

    if ( nodal_flag ) then
       mgt%ns = 1
    else
       if (mgt%stencil_type .eq. CC_CROSS_STENCIL) then
          ! Note that the cross stencils have 
          !      ns = 2*dm (lo/hi each direction) + 1 (center) + dm (extra used only at bc)
          if (mgt%dim .eq. 1) then
             mgt%ns = 4
          else if (mgt%dim .eq. 2) then
             mgt%ns = 7
          else if (mgt%dim .eq. 3) then
             mgt%ns = 10
          end if
       else if (mgt%stencil_type .eq. HO_CROSS_STENCIL) then
          if (mgt%dim .eq. 1) then
             call bl_error("HO_CROSS_STENCIL not supported in 1D")
          else if (mgt%dim .eq. 2) then
             mgt%ns = 9
          else if (mgt%dim .eq. 3) then
             mgt%ns = 13
          end if
       else if (mgt%stencil_type .eq. HO_DENSE_STENCIL) then
          if (mgt%dim .eq. 1) then
             call bl_error("HO_DENSE_STENCIL not supported in 1D")
          else if (mgt%dim .eq. 2) then
             mgt%ns = 25
          else if (mgt%dim .eq. 3) then
             mgt%ns = 61
          end if
       end if
    end if

    ng_for_res = 0; if ( nodal_flag ) ng_for_res = 1

    n = max_mg_levels(get_boxarray(la), nodal_flag, mgt%min_width)

    mgt%nlevels = min(n,mgt%max_nlevel) 

    n = mgt%nlevels

    allocate(mgt%cc(n), mgt%ff(n), mgt%dd(n), mgt%uu(n-1), mgt%ss(n), mgt%mm(n))
    allocate(mgt%pd(n),mgt%dh(mgt%dim,n))

    la1 = la

    do i = n, 1, -1
       call  multifab_build(mgt%cc(i), la1, mgt%nc, ng_for_res, nodal)
       call  multifab_build(mgt%ff(i), la1, mgt%nc, ng_for_res, nodal)
       call  multifab_build(mgt%dd(i), la1, mgt%nc, ng_for_res, nodal)

       ! This is the stencil of coefficients if doing a cell-centered solve
       ! This is the cell-centered coefficients if doing a nodal solve
       call  multifab_build(mgt%ss(i), la1, mgt%ns, ng_for_res, stencil = .true.)

       call setval(mgt%cc(i), zero, all = .TRUE.)
       call setval(mgt%ff(i), zero, all = .TRUE.)
       call setval(mgt%dd(i), zero, all = .TRUE.)
       !
       ! Set the stencil to zero; gotta do it by hand as multifab routines won't work.
       !
       !$omp parallel do private(j,p)
       do j = 1, nfabs(mgt%ss(i))
          p => dataptr(mgt%ss(i), j)
          p = zero
       end do
       !$omp end parallel do

       call imultifab_build(mgt%mm(i), la1, 1, 0, nodal)
       if ( i /= n ) &
          call multifab_build(mgt%uu(i), la1, mgt%nc, mgt%ng, nodal)
       
       mgt%pd(i) = coarsen(pd, 2**(n-i))
       if ( i > 1 ) call layout_build_coarse(la2, la1, (/(2,i=1,mgt%dim)/))
       la1 = la2
    end do

    mgt%uniform_dh = .true.
    ! These are real not integer so we have to be careful how we test on equality
    if ( present(dh) ) then
       mgt%dh(:,mgt%nlevels) = dh(:)
       select case ( mgt%dim )
       case (2)
          if (dh(1) > (ONE-1.0e-8_dp_t) * dh(2)  .and. &
              dh(1) < (ONE+1.0e-8_dp_t) * dh(2)) then
              mgt%uniform_dh = .true.
          else
              mgt%uniform_dh = .false.
          end if
       case (3)
          if (dh(1) > (ONE-1.0e-8_dp_t) * dh(2) .and. &
              dh(1) < (ONE+1.0e-8_dp_t) * dh(2) .and. & 
              dh(1) > (ONE-1.0e-8_dp_t) * dh(3) .and. &
              dh(1) < (ONE+1.0e-8_dp_t) * dh(3) ) then
              mgt%uniform_dh = .true.
          else
              mgt%uniform_dh = .false.
          end if
       end select
    else
       mgt%dh(:,mgt%nlevels) = ONE
    end if

    if ( (.not. mgt%uniform_dh) .and. nodal_flag ) &
       call bl_error("nodal solver does not support nonuniform dh")

    do i = mgt%nlevels-1, 1, -1
       mgt%dh(:,i) = mgt%dh(:,i+1)*TWO
    end do

    if ( .not. nodal_flag ) then
       allocate(mgt%skewed(mgt%nlevels,nfabs(mgt%cc(n))))
       allocate(mgt%skewed_not_set(mgt%nlevels))
       mgt%skewed_not_set = .true.
    end if
    !
    ! Save copy of domain BCs.
    !
    allocate(mgt%domain_bc(mgt%dim,1:2))

    mgt%domain_bc = domain_bc(:,1:2)

    if ( nodal_flag ) then
       !
       ! Set the face_type array to be BC_DIR or BC_NEU depending on domain_bc
       !
       allocate(mgt%face_type(nfabs(mgt%cc(n)),mgt%dim,2))
       mgt%face_type = BC_INT
       do id = 1,mgt%dim
          lo_dom = lwb(pd,id)
          hi_dom = upb(pd,id)
          do i = 1,nfabs(mgt%ss(mgt%nlevels))
             lo_grid =  lwb(get_box(mgt%ss(mgt%nlevels), i),id)
             if (lo_grid == lo_dom) mgt%face_type(i,id,1) = domain_bc(id,1)

             hi_grid = upb(get_box(mgt%ss(mgt%nlevels), i),id)
             if (hi_grid == hi_dom) mgt%face_type(i,id,2) = domain_bc(id,2)
          end do
       end do
    end if
    !
    ! Set both of these to false as default -- both must be true in order to subtract off sum(res)
    !
    mgt%bottom_singular    = .false.
    mgt%coeffs_sum_to_zero = .false.

    if ( present(is_singular) ) then
       mgt%bottom_singular = is_singular
    else

       !
       ! Do we cover the entire domain?
       ! Note that the volume is number of cells so it increments in units of 1, not dx*dy
       !
       ba = get_boxarray(get_layout(mgt%cc(mgt%nlevels)))
       
       if (int8_supported()) then
          vol    = boxarray_volume(ba)
          vol_pd = box_i8volume(pd)
          eq_vol = vol .eq. vol_pd
       else
          !
          ! Need these to be real, not int, so we can handle large numbers.
          !
          dvol    = boxarray_dvolume(ba)
          dvol_pd = box_dvolume(pd)
          eq_vol = abs(dvol-dvol_pd).lt.1.0e-2_dp_t
       end if

       ! If we cover the entire domain, then just test on the domain_bc values
       if ( eq_vol ) then
          mgt%bottom_singular = .true.
          do id = 1,mgt%dim
             if ( domain_bc(id,1) .eq. BC_DIR .or. domain_bc(id,2) .eq. BC_DIR ) &
                  mgt%bottom_singular = .false.
          end do
       else
          ! If we don't cover the entire domain, then this wont be singular
          mgt%bottom_singular = .false.
       end if
    end if

    ! if ( mgt%bottom_solver == 0 .and. .not. present(bottom_max_iter) ) mgt%bottom_max_iter = 20

    ! if only the bottom solver is 'solving' make sure that its eps is
    ! in effect
    if ( mgt%nlevels == 1 ) then
       ba = get_boxarray(get_layout(mgt%cc(1)))
       vol = boxarray_volume(ba)
       if (vol > 4**mgt%dim) &
         mgt%bottom_solver_eps = mgt%eps
    end if

    call destroy(bpt)

    ! If we're at a higher AMR level that coarsens to another mg_tower,
    !   or if you're here as part of a cc_applyop instead of a solve
    !   instead of coarsening within this one then don't bother
    !   creating the special bottom solver stuff

    if (mgt%nlevels == 1) mgt%bottom_solver = 1

    if (mgt%bottom_solver == 4 .and. mgt%use_hypre == 0) then

       ! This is the original layout
       old_coarse_la = get_layout(mgt%ss(1))

       bounding_box = bbox(get_boxarray(old_coarse_la));

       ! Does the boxarray fill the entire bounding box? If not then don't use this bottom solver
       ! also make sure there are at least 2**dm boxes or else you can't use the fancy bottom solver
       if (.not. contains(get_boxarray(old_coarse_la),bounding_box) &
            .or.  (nboxes(old_coarse_la) .lt. 2**mgt%dim) ) then

          mgt%bottom_solver = 1

          if ( parallel_IOProcessor() .and. verbose > 1 ) then
              print *,'F90mg: Do not use bottom_solver = 4 with this boxarray'
              print *,'F90mg: Using bottom_solver = 1 instead'
          end if

       else

           ! Get the old/new coarse problem domain
           coarse_pd = layout_get_pd(old_coarse_la)

           ! Get the new coarse boxarray and layout
           call box_build_2(bxs,bounding_box%lo(1:mgt%dim),bounding_box%hi(1:mgt%dim))
           call boxarray_build_bx(new_coarse_ba,bxs)
 
           ! compute the initial size of each of the boxes for the fancy
           ! bottom solver.  Each box will have bottom_box_size**dm cells
           success = 0
           call get_bottom_box_size(success,bottom_box_size,get_box(new_coarse_ba,1), &
                                    min_width, mgt%max_bottom_nlevel)

           if (success .eq. 1) then

               call boxarray_maxsize(new_coarse_ba,bottom_box_size)
               call layout_build_ba(new_coarse_la,new_coarse_ba,coarse_pd, &
                                    pmask=old_coarse_la%lap%pmask)
               call boxarray_destroy(new_coarse_ba)
               !
               ! Build a communicator on new_coarse_la%lap%prc.
               !
               communicator = parallel_create_communicator(new_coarse_la%lap%prc)

               if ( parallel_IOProcessor() .and. verbose > 1 ) then
                  print *,'F90mg: Coarse problem domain for bottom_solver = 4: '
                  print *,'   ... Bounding box is'
                  call print(bounding_box)
                  print *,'   ... Original boxes ',nboxes(old_coarse_la)
                  print *,'   ... New      boxes ',nboxes(new_coarse_la)
                  print *,'# cells on each side  ',bottom_box_size
               end if

               coarse_dx(:) = mgt%dh(:,1)

               allocate(mgt%bottom_mgt)
                  
               call bl_proffortfuncstop("mg_tower_build")
               call mg_tower_build(mgt%bottom_mgt, new_coarse_la, coarse_pd, &
                                   domain_bc, mgt%stencil_type, &
                                   dh = coarse_dx, &
                                   nc = mgt%nc, &
                                   ng = mgt%ng, &
                                   smoother = smoother, &
                                   nu1 = nu1, &
                                   nu2 = nu2, &
                                   nuf = nuf, &
                                   cycle_type = MG_VCycle, &
                                   bottom_solver = fancy_bottom_type, &
                                   bottom_max_iter = bottom_max_iter, &
                                   bottom_solver_eps = bottom_solver_eps, &
                                   max_iter = max_iter, &
                                   abort_on_max_iter = abort_on_max_iter, &
                                   max_nlevel = max_nlevel, &
                                   min_width = min_width, &
                                   eps = eps, &
                                   abs_eps = abs_eps, &
                                   max_L0_growth = max_L0_growth, &
                                   verbose = verbose, &
                                   cg_verbose = cg_verbose, &
                                   nodal = nodal, &
                                   the_bottom_comm = communicator)
               call bl_proffortfuncstart("mg_tower_build")
           else 

               mgt%bottom_solver = 1

               call boxarray_destroy(new_coarse_ba)

               if ( parallel_IOProcessor() .and. verbose > 1 ) then
                   print *,'F90mg: Unsuccessful in get_bottom_box_size'
                   print *,'F90mg: Setting bottom_solver to 1 instead'
               end if

           end if

       end if
    end if
    !
    ! We do this *after* the test on bottom_solver == 4 in case we redefine bottom_solver
    !    to be 1 or 2 in that test.
    !
    if ( nodal_flag .and. (mgt%bottom_solver == 1 .or. &
                           mgt%bottom_solver == 2 .or. &
                           mgt%bottom_solver == 3) ) then
       call build_nodal_dot_mask(mgt%nodal_mask,mgt%cc(1))
    end if

    call bl_proffortfuncstop("mg_tower_build")

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
    write(unit=un, fmt=*) 'cycle_type        = ', mgt%cycle_type
    call unit_skip(un, skip)
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
          do ii = 1,nfabs(mgt%cc(i))
             bb = get_box(mgt%cc(i),ii)
             if (mgt%dim == 1) then
                write(unit=un,fmt= '("  [",i4,"]: (",i4,") (",i4,")",i4 )') &
                     ii,bb%lo(1),bb%hi(1),bb%hi(1)-bb%lo(1)+1
             elseif (mgt%dim == 2) then
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

    logical      :: ldestroy_la
    type(layout) :: la
    integer      :: i

    ldestroy_la = .false.; if ( present(destroy_la) ) ldestroy_la = destroy_la

    la = get_layout(mgt%cc(mgt%nlevels))

    do i = 1, mgt%nlevels
       call multifab_destroy(mgt%cc(i))
       call multifab_destroy(mgt%ff(i))
       call multifab_destroy(mgt%dd(i))
       call multifab_destroy(mgt%ss(i))
       call imultifab_destroy(mgt%mm(i))
       if ( i /= mgt%nlevels ) then
          call multifab_destroy(mgt%uu(i))
       end if
    end do

    deallocate(mgt%cc, mgt%ff, mgt%dd, mgt%uu, mgt%mm, mgt%ss)
    deallocate(mgt%dh, mgt%pd)
    deallocate(mgt%domain_bc)

    if ( associated(mgt%face_type) ) deallocate(mgt%face_type)

    if ( associated(mgt%skewed) ) then
       deallocate(mgt%skewed)
       deallocate(mgt%skewed_not_set)
    end if

    if ( built_q(mgt%nodal_mask) ) call multifab_destroy(mgt%nodal_mask)

    if ( ldestroy_la ) call layout_destroy(la)

    if ( associated(mgt%bottom_comm) ) then
       call parallel_free_communicator(mgt%bottom_comm)
       deallocate(mgt%bottom_comm)
    end if

    if ( associated(mgt%bottom_mgt) ) then
       call mg_tower_destroy(mgt%bottom_mgt, .true.)
       deallocate(mgt%bottom_mgt)
    end if

  end subroutine mg_tower_destroy

  function max_mg_levels(ba, nodal_flag, min_size) result(r)

    type(boxarray), intent(in)           :: ba
    logical       , intent(in)           :: nodal_flag
    integer       , intent(in), optional :: min_size
    integer                              :: r

    integer, parameter :: rrr = 2
    type(box)          :: bx, bx1
    type(boxarray)     :: ba1
    real(kind=dp_t)    :: vol
    integer            :: i, rr, lmn, dm

    dm = get_dim(ba)

    ! lmn is the smallest we allow any side of a box to become
    lmn = 1; if ( present(min_size) ) lmn = min_size

    ! r keep track of the total number of mg levels
    r = 1

    ! rr keeps track of ref ratio between first and last mg levels we are currently testing
    rr = rrr

    ! create a temporary copy of the boxarray
    call boxarray_build_copy(ba1,ba)

    outer: do
       call boxarray_coarsen(ba1,rrr)
       vol = boxarray_volume(ba1)
       do i = 1, nboxes(ba)
          bx = get_box(ba,i)
          bx1 = coarsen(bx, rr)

          ! check to see if any side has become less than lmn
          if ( any(extent(bx1) < lmn) ) exit outer

          ! check to see if bx is evenly divisible by rr
          if ( bx /= refine(bx1, rr)  ) exit outer

          ! check to make sure the entire problem over all grids is not too small
          ! we prefer slightly larger 'bottom solves' for nodal problems since certain
          ! bc types do not have degrees of freedom on the domain boundary and you
          ! can end up with a 1 point bottom solve, even if the grid is 2**dm
          if (nodal_flag) then
             if ( vol <= 2**dm ) exit outer  
          else
             if ( vol < 2**dm ) exit outer
          end if
          
       end do
       rr = rr*rrr
       r  = r + 1
    end do outer

    call boxarray_destroy(ba1)

  end function max_mg_levels

  subroutine get_bottom_box_size(success,bottom_box_size,bxs,min_size,max_bottom_nlevel)

    integer  , intent(inout) :: success
    integer  , intent(  out) :: bottom_box_size
    type(box), intent(in   ) :: bxs
    integer  , intent(in   ) :: min_size
    integer  , intent(in   ) :: max_bottom_nlevel

    type(box) :: bx, bx1
    integer   :: rr
    integer   :: bottom_levs

    rr          = 2
    success     = 1
    bottom_levs = 1

    do
       bx  = bxs
       bx1 = coarsen(bx, rr)

       if ( any(extent(bx1) < min_size) ) exit

       if ( any(mod(extent(bx1),2) .eq. 1) ) then

          if ( all(mod(extent(bx1),3) .eq. 0) ) then
             ! test 3
             bottom_levs  = bottom_levs + 1
             rr = rr*3
             exit
          else if ( all(mod(extent(bx1),5) .eq. 0) ) then
             ! test 5
             bottom_levs  = bottom_levs + 1
             rr = rr*5
             exit
          else
             exit
          end if

       end if

       if ( bx /= refine(bx1, rr)  ) then
          !
          ! This means bx has an odd number of cells on one or more side
          ! we can only get here if rr=2, i.e., the first pass through this do loop
          ! we exit the do loop with bottom_levs=1 and will abort.
          !
          exit
       end if

       rr = rr*2

       bottom_levs = bottom_levs + 1
       !
       ! Max size of grids will be 2**max(bottom_levs,max_bottom_nlevel)
       !
       if (bottom_levs .eq. max_bottom_nlevel) exit

    end do

    bottom_box_size = rr

    if ( bottom_levs .eq. 1) success = 0

  end subroutine get_bottom_box_size

  subroutine do_bottom_mgt(mgt, uu, rh)

    use bl_prof_module

    type( mg_tower), intent(inout) :: mgt
    type( multifab), intent(inout) :: uu
    type( multifab), intent(in   ) :: rh

    type(bl_prof_timer), save :: bpt

    type( multifab ) :: bottom_uu
    type( multifab ) :: bottom_rh
    type( layout   ) :: la
    integer          :: mglev
    logical          :: do_diag

    do_diag = .false.; if ( mgt%verbose >= 4 ) do_diag = .true.

    call build(bpt, "do_bottom_mgt")

    if ( mgt%bottom_solver .ne. 4 ) &
       call bl_error("MG_TOWER_BOTTOM_SOLVE: must have bottom_solver == 4")

    mglev = mgt%bottom_mgt%nlevels

    la = get_layout(mgt%bottom_mgt%ss(mglev))

    call multifab_build(bottom_uu,la,1,nghost(uu),nodal_flags(uu))

    call setval(bottom_uu,ZERO,all=.true.)

    if (nodal_q(rh)) then
       call multifab_build(bottom_rh,la,1,1,nodal_flags(rh))
       call setval(bottom_rh,ZERO,all=.true.)
    else
       call multifab_build(bottom_rh,la,1,0,nodal_flags(rh))
    end if

    call multifab_copy_c(bottom_rh,1,rh,1,1,ng=0)

    call mg_tower_cycle(mgt%bottom_mgt, mgt%bottom_mgt%cycle_type, mglev, &
                        mgt%bottom_mgt%ss(mglev), bottom_uu, bottom_rh, &
                        mgt%bottom_mgt%mm(mglev), mgt%bottom_mgt%nu1, mgt%bottom_mgt%nu2)

    call multifab_copy_c(uu,1,bottom_uu,1,1,ng=0)

    call multifab_destroy(bottom_uu)
    call multifab_destroy(bottom_rh)

    call destroy(bpt)

  end subroutine do_bottom_mgt

  subroutine mg_tower_bottom_solve(mgt, lev, ss, uu, rh, mm, eps_in)

    use bl_prof_module
    use itsol_module, only: itsol_bicgstab_solve, itsol_cabicgstab_solve, itsol_cg_solve

    type( mg_tower), intent(inout) :: mgt
    type( multifab), intent(inout) :: uu
    type( multifab), intent(in) :: rh
    type( multifab), intent(in) :: ss
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: lev
    real(dp_t), intent(in), optional :: eps_in

    integer             :: i,stat,communicator,tag
    logical             :: singular_test,do_diag
    real(dp_t)          :: nrm, eps

    type(bl_prof_timer), save :: bpt

    call bl_proffortfuncstart("mg_tower_bottom_solve")
    call build(bpt, "mgt_bottom_solve")

    do_diag = .false.; if ( mgt%verbose >= 4 ) do_diag = .true.

    eps = mgt%bottom_solver_eps ; if ( present(eps_in) ) eps = eps_in

    if ( .not.nodal_q(rh) ) then
       singular_test =  mgt%bottom_singular .and. mgt%coeffs_sum_to_zero
    end if

    if ( do_diag ) then
       call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)
       nrm = norm_inf(mgt%cc(lev))
       if ( parallel_IOProcessor() ) &
          print *,' BOT: Norm before bottom         ',nrm
    end if

    if ( associated(mgt%bottom_comm) ) then
       communicator = mgt%bottom_comm
    else
       communicator = parallel_communicator()
    end if

    stat = 0

    select case ( mgt%bottom_solver )
    case (0)
       do i = 1, mgt%nuf
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do
    case (1)
       if ( nodal_q(rh) ) then
          call itsol_bicgstab_solve(ss, uu, rh, mm, &
                                    eps, mgt%bottom_max_iter, &
                                    mgt%cg_verbose, &
                                    mgt%stencil_type, mgt%lcross, &
                                    stat = stat, &
                                    singular_in = mgt%bottom_singular, &
                                    uniform_dh = mgt%uniform_dh,&
                                    nodal_mask = mgt%nodal_mask, &
                                    comm_in = communicator)
       else
          call itsol_bicgstab_solve(ss, uu, rh, mm, &
                                    eps, mgt%bottom_max_iter, &
                                    mgt%cg_verbose,  &
                                    mgt%stencil_type, mgt%lcross, &
                                    stat = stat, &
                                    singular_in = singular_test, &
                                    uniform_dh = mgt%uniform_dh, &
                                    comm_in = communicator)
       end if
       do i = 1, mgt%nub
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do
    case (2)
       if ( nodal_q(rh) ) then
          call itsol_cg_solve(ss, uu, rh, mm, &
                              eps, mgt%bottom_max_iter, mgt%cg_verbose, &
                              mgt%stencil_type, mgt%lcross, &
                              stat = stat, singular_in = mgt%bottom_singular, &
                              uniform_dh = mgt%uniform_dh, nodal_mask=mgt%nodal_mask)
       else
          call itsol_cg_solve(ss, uu, rh, mm, &
                              eps, mgt%bottom_max_iter, mgt%cg_verbose, &
                              mgt%stencil_type, mgt%lcross, &
                              stat = stat, singular_in = singular_test, &
                              uniform_dh = mgt%uniform_dh)
       end if
       do i = 1, mgt%nub
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do
    case (3)
       if ( nodal_q(rh) ) then
          call itsol_cabicgstab_solve(ss, uu, rh, mm, &
                                      eps, mgt%bottom_max_iter, &
                                      mgt%cg_verbose, &
                                      mgt%stencil_type, mgt%lcross, &
                                      stat = stat, &
                                      singular_in = mgt%bottom_singular, &
                                      uniform_dh = mgt%uniform_dh,&
                                      nodal_mask = mgt%nodal_mask, &
                                      comm_in = communicator)
       else
          call itsol_cabicgstab_solve(ss, uu, rh, mm, &
                                      eps, mgt%bottom_max_iter, &
                                      mgt%cg_verbose,  &
                                      mgt%stencil_type, mgt%lcross, &
                                      stat = stat, &
                                      singular_in = singular_test, &
                                      uniform_dh = mgt%uniform_dh, &
                                      comm_in = communicator)
       end if
       do i = 1, mgt%nub
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do

    case default
       call bl_error("MG_TOWER_BOTTOM_SOLVE: no such solver: ", mgt%bottom_solver)
    end select

    if (associated(mgt%bottom_comm)) then
       tag = parallel_tag(reset=.true.)
    end if

    if ( stat /= 0 ) then
       if ( parallel_IOProcessor() .and. mgt%verbose > 0 ) then
          call bl_warn("BREAKDOWN in bottom_solver: trying smoother")
       end if
       do i = 1, 20
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do
    end if

    if ( do_diag ) then
       call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)
       nrm = norm_inf(mgt%cc(lev))
       if ( parallel_IOProcessor() ) &
          print *,' BOT: Norm  after bottom         ',nrm
    end if

    call destroy(bpt)
    call bl_proffortfuncstop("mg_tower_bottom_solve")

  end subroutine mg_tower_bottom_solve

  subroutine mg_tower_restriction(mgt, crse, fine, mm_fine, mm_crse)

    use bl_prof_module
    use cc_restriction_module, only: cc_restriction_1d, cc_restriction_2d, cc_restriction_3d
    use nodal_restriction_module, only: nodal_restriction_1d, nodal_restriction_2d, &
                                        nodal_restriction_3d
    use impose_neumann_bcs_module

    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    type(imultifab), intent(in)   :: mm_fine
    type(imultifab), intent(in)   :: mm_crse
    type(mg_tower), intent(inout) :: mgt
    integer :: loc(mgt%dim), lof(mgt%dim)
    integer :: lom_fine(mgt%dim)
    integer :: lom_crse(mgt%dim)
    integer ::  lo(mgt%dim),  hi(mgt%dim)
    integer :: vlo(mgt%dim), vhi(mgt%dim)
    integer :: i, n, ir(mgt%dim)
    logical :: nodal_flag
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    integer        , pointer :: mp_fine(:,:,:,:)
    integer        , pointer :: mp_crse(:,:,:,:)

    integer :: mg_restriction_mode, ng
    type(bl_prof_timer), save :: bpt
    type(mfiter) :: mfi
    integer :: tlo(mgt%dim), thi(mgt%dim)
    type(box) :: tilebox
    logical :: do_tiling

    call bl_proffortfuncstart("mg_tower_restriction")
    call build(bpt, "mgt_restriction")

    ir = 2

    nodal_flag = nodal_q(crse)

    if ( nodal_flag ) then
       call multifab_fill_boundary(fine)
       call setval(crse, ZERO, all=.true.)
       mg_restriction_mode = 1
    
       !$omp parallel do private(i,fp,mp_fine,lof,lom_fine,ng,n)
       do i = 1, nfabs(fine)
          fp      => dataptr(fine,    i)
          mp_fine => dataptr(mm_fine, i)
          lof      = lwb(get_pbox(fine,    i))
          lom_fine = lwb(get_pbox(mm_fine, i))
          ng       = lom_fine(1) - lof(1)
          do n = 1, ncomp(fine)
             select case (mgt%dim)
             case (1)
                call impose_neumann_bcs_1d(fp(:,1,1,n),mp_fine(:,1,1,1),lom_fine,ng)
             case (2)
                call impose_neumann_bcs_2d(fp(:,:,1,n),mp_fine(:,:,1,1),lom_fine,ng)
             case (3)
                call impose_neumann_bcs_3d(fp(:,:,:,n),mp_fine(:,:,:,1),lom_fine,ng)
             end select
          end do
       end do
       !$omp end parallel do
    end if

    do_tiling = mgt%dim.eq.3 .and. .not.nodal_flag

    !$omp parallel default(none) &
    !$omp  private(i,mfi,tilebox,tlo,thi,cp,fp,mp_fine,mp_crse,loc,lof,lom_fine,lom_crse,lo,hi,vlo,vhi) &
    !$omp  shared(crse,fine,mm_fine,mm_crse,mgt,nodal_flag,ir,mg_restriction_mode,do_tiling)
    call mfiter_build(mfi, crse, tiling=do_tiling)
    do while(next_tile(mfi,i))

       tilebox = get_tilebox(mfi)
       tlo = lwb(tilebox)
       thi = upb(tilebox)

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
       vlo      = lo
       vhi      = hi

       do n = 1, mgt%nc
          select case ( mgt%dim )
          case (1)
             if ( .not. nodal_flag ) then
                call cc_restriction_1d(cp(:,1,1,n), loc, fp(:,1,1,n), lof, lo, hi, ir)
             else
                call nodal_restriction_1d(cp(:,1,1,n), loc, fp(:,1,1,n), lof, &
                                          mp_fine(:,1,1,1), lom_fine, &
                                          mp_crse(:,1,1,1), lom_crse, lo, hi, vlo, vhi, ir, .false., &
                                          mg_restriction_mode)
             end if
          case (2)
             if ( .not. nodal_flag ) then
                call cc_restriction_2d(cp(:,:,1,n), loc, fp(:,:,1,n), lof, lo, hi, ir)
             else
                call nodal_restriction_2d(cp(:,:,1,n), loc, fp(:,:,1,n), lof, &
                                          mp_fine(:,:,1,1), lom_fine, &
                                          mp_crse(:,:,1,1), lom_crse, lo, hi, vlo, vhi, ir, .false., &
                                          mg_restriction_mode)
             end if
          case (3)
             if ( .not. nodal_flag ) then
                call cc_restriction_3d(cp(:,:,:,n), loc, fp(:,:,:,n), lof, tlo, thi, ir)
             else
                call nodal_restriction_3d(cp(:,:,:,n), loc, fp(:,:,:,n), lof, &
                                          mp_fine(:,:,:,1), lom_fine, &
                                          mp_crse(:,:,:,1), lom_crse, lo, hi, vlo, vhi, ir, .false., &
                                          mg_restriction_mode)
             end if
          end select
       end do
    end do
    !$OMP END PARALLEL

    call destroy(bpt)
    call bl_proffortfuncstop("mg_tower_restriction")

  end subroutine mg_tower_restriction
  !
  ! In the following two "impose" routines we assume that any ghost cells
  ! covered by valid region have been filled properly by fill_boundary().
  ! Also, that "all" ghost cells have computable values.  We may copy some
  ! ghost cells only to overwrite them later with "correct" values.
  !
  subroutine impose_physbc_cc_2d(uu, mgt, lo, hi, dlo, dhi, ng)

    type(mg_tower),  intent(in   ) :: mgt
    integer,         intent(in   ) :: lo(:), hi(:), dlo(:), dhi(:), ng
    real(kind=dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:,:)

    if ( lo(1) == dlo(1) ) then

       if ( mgt%domain_bc(1,1) == BC_DIR ) then

          uu(lo(1)-1,:,1:mgt%nc) = -uu(lo(1),:,1:mgt%nc)

       else if ( mgt%domain_bc(1,1) == BC_NEU ) then

          uu(lo(1)-1,:,1:mgt%nc) =  uu(lo(1),:,1:mgt%nc)

       end if

       if ( mgt%domain_bc(1,1) == BC_DIR .or. mgt%domain_bc(1,1) == BC_NEU ) then
          !
          ! Corners
          !
          if ( lo(2) == dlo(2) ) then
             uu(lo(1)-1,lo(2)-1,1:mgt%nc) = uu(lo(1)-1,lo(2),1:mgt%nc)
          end if
          if ( hi(2) == dhi(2) ) then
             uu(lo(1)-1,hi(2)+1,1:mgt%nc) = uu(lo(1)-1,hi(2),1:mgt%nc)
          end if
       end if

    end if

    if ( hi(1) == dhi(1) ) then

       if ( mgt%domain_bc(1,2) == BC_DIR ) then

          uu(hi(1)+1,lo(2):hi(2),1:mgt%nc) = -uu(hi(1),lo(2):hi(2),1:mgt%nc)

       else if ( mgt%domain_bc(1,2) == BC_NEU ) then

          uu(hi(1)+1,lo(2):hi(2),1:mgt%nc) =  uu(hi(1),lo(2):hi(2),1:mgt%nc)

       end if

       if ( mgt%domain_bc(1,2) == BC_DIR .or. mgt%domain_bc(1,2) == BC_NEU ) then
          !
          ! Corners
          !
          if ( lo(2) == dlo(2) ) then
             uu(hi(1)+1,lo(2)-1,1:mgt%nc) = uu(hi(1)+1,lo(2),1:mgt%nc)
          end if
          if ( hi(2) == dhi(2) ) then
             uu(hi(1)+1,hi(2)+1,1:mgt%nc) = uu(hi(1)+1,hi(2),1:mgt%nc)
          end if
       end if

    end if
    !
    ! At this point all possible corner cells have been filled.
    !
    if ( lo(2) == dlo(2) ) then

       if ( mgt%domain_bc(2,1) == BC_DIR ) then

          uu(lo(1):hi(1),lo(2)-1,1:mgt%nc) = -uu(lo(1):hi(1),lo(2),1:mgt%nc)

       else if ( mgt%domain_bc(2,1) == BC_NEU ) then

          uu(lo(1):hi(1),lo(2)-1,1:mgt%nc) =  uu(lo(1):hi(1),lo(2),1:mgt%nc)

       end if

    end if

    if ( hi(2) == dhi(2) ) then

       if ( mgt%domain_bc(2,2) == BC_DIR ) then

          uu(lo(1):hi(1),hi(2)+1,1:mgt%nc) = -uu(lo(1):hi(1),hi(2),1:mgt%nc)

       else if ( mgt%domain_bc(2,2) == BC_NEU ) then

          uu(lo(1):hi(1),hi(2)+1,1:mgt%nc) =  uu(lo(1):hi(1),hi(2),1:mgt%nc)

       end if

    end if

  end subroutine impose_physbc_cc_2d

  subroutine impose_physbc_cc_3d(uu, mgt, lo, hi, dlo, dhi, ng)

    type(mg_tower),  intent(in   ) :: mgt
    integer,         intent(in   ) :: lo(:), hi(:), dlo(:), dhi(:), ng
    real(kind=dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)

    if ( lo(1) == dlo(1) ) then

       if ( mgt%domain_bc(1,1) == BC_DIR ) then

          uu(lo(1)-1,:,:,1:mgt%nc) = -uu(lo(1),:,:,1:mgt%nc)

       else if ( mgt%domain_bc(1,1) == BC_NEU ) then

          uu(lo(1)-1,:,:,1:mgt%nc) =  uu(lo(1),:,:,1:mgt%nc)

       end if

       if ( mgt%domain_bc(1,1) == BC_DIR .or. mgt%domain_bc(1,1) == BC_NEU ) then
          !
          ! Edges.
          !
          if ( lo(2) == dlo(2) ) then
             uu(lo(1)-1,lo(2)-1,lo(3):hi(3),1:mgt%nc) = uu(lo(1)-1,lo(2),lo(3):hi(3),1:mgt%nc)
          end if
          if ( hi(2) == dhi(2) ) then
             uu(lo(1)-1,hi(2)+1,lo(3):hi(3),1:mgt%nc) = uu(lo(1)-1,hi(2),lo(3):hi(3),1:mgt%nc)
          end if
          if ( lo(3) == dlo(3) ) then
             uu(lo(1)-1,lo(2):hi(2),lo(3)-1,1:mgt%nc) = uu(lo(1)-1,lo(2):hi(2),lo(3),1:mgt%nc)
          end if
          if ( hi(3) == dhi(3) ) then
             uu(lo(1)-1,lo(2):hi(2),hi(3)+1,1:mgt%nc) = uu(lo(1)-1,lo(2):hi(2),hi(3),1:mgt%nc)
          end if
          !
          ! Corners
          !
          if ( lo(2) == dlo(2) .and. lo(3) == dlo(3) ) then
             uu(lo(1)-1,lo(2)-1,lo(3)-1,1:mgt%nc) = uu(lo(1)-1,lo(2),lo(3),1:mgt%nc)
          end if
          if ( lo(2) == dlo(2) .and. hi(3) == dhi(3) ) then
             uu(lo(1)-1,lo(2)-1,hi(3)+1,1:mgt%nc) = uu(lo(1)-1,lo(2),hi(3),1:mgt%nc)
          end if
          if ( hi(2) == dhi(2) .and. lo(3) == dlo(3) ) then
             uu(lo(1)-1,hi(2)+1,lo(3)-1,1:mgt%nc) = uu(lo(1)-1,hi(2),lo(3),1:mgt%nc)
          end if
          if ( hi(2) == dhi(2) .and. hi(3) == dhi(3) ) then
             uu(lo(1)-1,hi(2)+1,hi(3)+1,1:mgt%nc) = uu(lo(1)-1,hi(2),hi(3),1:mgt%nc)
          end if
       end if

    end if

    if ( hi(1) == dhi(1) ) then

       if ( mgt%domain_bc(1,2) == BC_DIR ) then

          uu(hi(1)+1,:,:,1:mgt%nc) = -uu(hi(1),:,:,1:mgt%nc)

       else if ( mgt%domain_bc(1,2) == BC_NEU ) then

          uu(hi(1)+1,:,:,1:mgt%nc) =  uu(hi(1),:,:,1:mgt%nc)

       end if

       if ( mgt%domain_bc(1,2) == BC_DIR .or. mgt%domain_bc(1,2) == BC_NEU ) then
          !
          ! Edges.
          !
          if ( lo(2) == dlo(2) ) then
             uu(hi(1)+1,lo(2)-1,lo(3):hi(3),1:mgt%nc) = uu(hi(1)+1,lo(2),lo(3):hi(3),1:mgt%nc)
          end if
          if ( hi(2) == dhi(2) ) then
             uu(hi(1)+1,hi(2)+1,lo(3):hi(3),1:mgt%nc) = uu(hi(1)+1,hi(2),lo(3):hi(3),1:mgt%nc)
          end if
          if ( lo(3) == dlo(3) ) then
             uu(hi(1)+1,lo(2):hi(2),lo(3)-1,1:mgt%nc) = uu(hi(1)+1,lo(2):hi(2),lo(3),1:mgt%nc)
          end if
          if ( hi(3) == dhi(3) ) then
             uu(hi(1)+1,lo(2):hi(2),hi(3)+1,1:mgt%nc) = uu(hi(1)+1,lo(2):hi(2),hi(3),1:mgt%nc)
          end if
          !
          ! Corners.
          !
          if ( lo(2) == dlo(2) .and. lo(3) == dlo(3) ) then
             uu(hi(1)+1,lo(2)-1,lo(3)-1,1:mgt%nc) = uu(hi(1)+1,lo(2),lo(3),1:mgt%nc)
          end if
          if ( lo(2) == dlo(2) .and. hi(3) == dhi(3) ) then
             uu(hi(1)+1,lo(2)-1,hi(3)+1,1:mgt%nc) = uu(hi(1)+1,lo(2),hi(3),1:mgt%nc)
          end if
          if ( hi(2) == dhi(2) .and. lo(3) == dlo(3) ) then
             uu(hi(1)+1,hi(2)+1,lo(3)-1,1:mgt%nc) = uu(hi(1)+1,hi(2),lo(3),1:mgt%nc)
          end if
          if ( hi(2) == dhi(2) .and. hi(3) == dhi(3) ) then
             uu(hi(1)+1,hi(2)+1,hi(3)+1,1:mgt%nc) = uu(hi(1)+1,hi(2),hi(3),1:mgt%nc)
          end if
       end if

    end if
    !
    ! All possible corners should be done by this point. Only some faces & edges remain.
    !
    if ( lo(2) == dlo(2) ) then

       if ( mgt%domain_bc(2,1) == BC_DIR ) then

          uu(lo(1):hi(1),lo(2)-1,:,1:mgt%nc) = -uu(lo(1):hi(1),lo(2),:,1:mgt%nc)

       else if ( mgt%domain_bc(2,1) == BC_NEU ) then

          uu(lo(1):hi(1),lo(2)-1,:,1:mgt%nc) =  uu(lo(1):hi(1),lo(2),:,1:mgt%nc)

       end if

       if ( mgt%domain_bc(2,1) == BC_DIR .or. mgt%domain_bc(2,1) == BC_NEU ) then
          !
          ! Edges - only need to do lo & hi Z.
          !
          if ( lo(3) == dlo(3) ) then
             uu(lo(1):hi(1),lo(2)-1,lo(3)-1,1:mgt%nc) = uu(lo(1):hi(1),lo(2)-1,lo(3),1:mgt%nc)
          end if
          if ( hi(3) == dhi(3) ) then
             uu(lo(1):hi(1),lo(2)-1,hi(3)+1,1:mgt%nc) = uu(lo(1):hi(1),lo(2)-1,hi(3),1:mgt%nc)
          end if
       end if

    end if

    if ( hi(2) == dhi(2) ) then

       if ( mgt%domain_bc(2,2) == BC_DIR ) then

          uu(lo(1):hi(1),hi(2)+1,:,1:mgt%nc) = -uu(lo(1):hi(1),hi(2),:,1:mgt%nc)

       else if ( mgt%domain_bc(2,2) == BC_NEU ) then

          uu(lo(1):hi(1),hi(2)+1,:,1:mgt%nc) =  uu(lo(1):hi(1),hi(2),:,1:mgt%nc)

       end if

       if ( mgt%domain_bc(2,2) == BC_DIR .or. mgt%domain_bc(2,2) == BC_NEU ) then
          !
          ! Edges - only need to do lo & hi Z.
          !
          if ( lo(3) == dlo(3) ) then
             uu(lo(1):hi(1),hi(2)+1,lo(3)-1,1:mgt%nc) = uu(lo(1):hi(1),hi(2)+1,lo(3),1:mgt%nc)
          end if
          if ( hi(3) == dhi(3) ) then
             uu(lo(1):hi(1),hi(2)+1,hi(3)+1,1:mgt%nc) = uu(lo(1):hi(1),hi(2)+1,hi(3),1:mgt%nc)
          end if
       end if

    end if
    !
    ! All that remains to do are any Z faces.
    !
    if ( lo(3) == dlo(3) ) then

       if ( mgt%domain_bc(3,1) == BC_DIR ) then

          uu(lo(1):hi(1),lo(2):hi(2),lo(3)-1,1:mgt%nc) = -uu(lo(1):hi(1),lo(2):hi(2),lo(3),1:mgt%nc)

       else if ( mgt%domain_bc(3,1) == BC_NEU ) then

          uu(lo(1):hi(1),lo(2):hi(2),lo(3)-1,1:mgt%nc) =  uu(lo(1):hi(1),lo(2):hi(2),lo(3),1:mgt%nc)

       end if

    end if

    if ( hi(3) == dhi(3) ) then

       if ( mgt%domain_bc(3,2) == BC_DIR ) then

          uu(lo(1):hi(1),lo(2):hi(2),hi(3)+1,1:mgt%nc) = -uu(lo(1):hi(1),lo(2):hi(2),hi(3),1:mgt%nc)

       else if ( mgt%domain_bc(3,2) == BC_NEU ) then

          uu(lo(1):hi(1),lo(2):hi(2),hi(3)+1,1:mgt%nc) =  uu(lo(1):hi(1),lo(2):hi(2),hi(3),1:mgt%nc)

       end if

    end if

  end subroutine impose_physbc_cc_3d
  !
  ! Prolongate from mgt%uu(lev) -> uu.
  !
  subroutine mg_tower_prolongation(mgt, uu, lev)
    use bl_prof_module
    use mg_prolongation_module
    use impose_neumann_bcs_module

    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(inout) :: uu
    integer       , intent(in   ) :: lev

    real(kind=dp_t), pointer  :: fp(:,:,:,:), cp(:,:,:,:), up(:,:,:,:)
    integer,         pointer  :: mp(:,:,:,:)
    integer                   :: i, j, k, n, ng, ir(mgt%dim)
    integer                   :: lo(mgt%dim), hi(mgt%dim), lom(mgt%dim)
    integer                   :: loc(mgt%dim), lof(mgt%dim), dlo(mgt%dim), dhi(mgt%dim)
    type(box)                 :: dmn
    type(bl_prof_timer), save :: bpt

    call bl_proffortfuncstart("mg_tower_prolongation")
    call build(bpt, "mgt_prolongation")

    call bl_assert( mgt%nc == ncomp(uu)         , 'mg_tower_prolongation: ncomp')
    call bl_assert( mgt%nc == ncomp(mgt%uu(lev)), 'mg_tower_prolongation: ncomp')

    ir  = 2
    ng  = nghost(mgt%uu(lev))
    dmn = get_pd(get_layout(mgt%uu(lev)))
    dlo = lwb(dmn)
    dhi = upb(dmn)

    if ( .not. nodal_q(uu) ) then
       if ( mgt%use_lininterp ) then
          if ( ng > 0 ) then
             !
             ! Set up dirichlet/neumann boundaries so lininterp does right thing.
             ! If we don't have one ghost cell the interp routines will do a
             ! piecewise constant interp on the cells touching the grid boundary.
             !
             call multifab_fill_boundary(mgt%uu(lev))

             do n = 1, nfabs(uu)
                up => dataptr(mgt%uu(lev) ,n)
                lo =  lwb(get_ibox(mgt%uu(lev), n))
                hi =  upb(get_ibox(mgt%uu(lev), n))
                select case ( mgt%dim )
                case (2)
                   call impose_physbc_cc_2d(up(:,:,1,:),mgt,lo,hi,dlo,dhi,ng)
                case (3)
                   call impose_physbc_cc_3d(up(:,:,:,:),mgt,lo,hi,dlo,dhi,ng)
                end select
             end do
          end if
       end if

       !$OMP PARALLEL DO PRIVATE(i,n,loc,lof,lo,hi,fp,cp)
       do i = 1, nfabs(uu)
          loc =  lwb(get_pbox(mgt%uu(lev),i))
          lof =  lwb(get_pbox(uu, i))
          lo  =  lwb(get_ibox(uu, i))
          hi  =  upb(get_ibox(uu, i))
          fp  => dataptr(uu,         i)
          cp  => dataptr(mgt%uu(lev),i)
          do n = 1, mgt%nc
             select case ( mgt%dim )
             case (1)
                call cc_prolongation(fp(:,1,1,n), lof, cp(:,1,1,n), loc, lo, hi, ir, mgt%use_lininterp)
             case (2)
                call cc_prolongation(fp(:,:,1,n), lof, cp(:,:,1,n), loc, lo, hi, ir, mgt%use_lininterp, mgt%ptype)
             case (3)
                call cc_prolongation(fp(:,:,:,n), lof, cp(:,:,:,n), loc, lo, hi, ir, mgt%use_lininterp, mgt%ptype)
             end select
          end do
       end do
       !$OMP END PARALLEL DO

    else

      if ( using_nodal_cubic() .and. mgt%dim > 1 ) then
         call multifab_fill_boundary(mgt%uu(lev))

         do i = 1, nfabs(mgt%uu(lev))
            up  => dataptr(mgt%uu(lev), i)
            mp  => dataptr(mgt%mm(lev), i)
            lom =  lwb(get_ibox(mgt%mm(lev),i))
            lo  =  lwb(get_pbox(mgt%uu(lev),i))
            ng  =  lom(1) - lo(1)
            do n = 1, mgt%nc
             select case ( mgt%dim )
             case (2)
                call impose_neumann_bcs_2d(up(:,:,1,n),mp(:,:,1,1),lom,ng)
             case (3)
                call impose_neumann_bcs_3d(up(:,:,:,n),mp(:,:,:,1),lom,ng)
             end select
            end do
         end do
      end if

       !$OMP PARALLEL DO PRIVATE(i,n,loc,lof,lo,hi,fp,cp)
       do i = 1, nfabs(uu)
          loc =  lwb(get_pbox(mgt%uu(lev),i))
          lof =  lwb(get_pbox(uu, i))
          lo  =  lwb(get_ibox(uu, i))
          hi  =  upb(get_ibox(uu, i))
          fp  => dataptr(uu,         i)
          cp  => dataptr(mgt%uu(lev),i)
          do n = 1, mgt%nc
             select case ( mgt%dim )
             case (1)
                call nodal_prolongation_1d(fp(:,1,1,n), lof, cp(:,1,1,n), loc, lo, hi, ir)
             case (2)
                call nodal_prolongation_2d(fp(:,:,1,n), lof, cp(:,:,1,n), loc, lo, hi, ir)
             case (3)
                call nodal_prolongation_3d(fp(:,:,:,n), lof, cp(:,:,:,n), loc, lo, hi, ir)
             end select
          end do
       end do
       !$OMP END PARALLEL DO

       if ( using_nodal_cubic() .and. mgt%dim > 1 ) then
          !
          ! The [bi,tri]cubic interpolators don't preserve dirichlet BCs.
          !
          do n = 1, nfabs(uu)
             up  => dataptr(uu           ,n)
             mp  => dataptr(mgt%mm(lev+1),n)
             lo  =  lwb(get_ibox(uu, n))
             hi  =  upb(get_ibox(uu, n))

             select case ( mgt%dim )
             case (2)
                do j = lo(2),hi(2)
                   if ( bc_dirichlet(mp(lo(1),j,1,1),1,0) ) up(lo(1),j,1,1:mgt%nc) = ZERO
                   if ( bc_dirichlet(mp(hi(1),j,1,1),1,0) ) up(hi(1),j,1,1:mgt%nc) = ZERO
                end do

                do i = lo(1),hi(1)
                   if ( bc_dirichlet(mp(i,lo(2),1,1),1,0) ) up(i,lo(2),1,1:mgt%nc) = ZERO
                   if ( bc_dirichlet(mp(i,hi(2),1,1),1,0) ) up(i,hi(2),1,1:mgt%nc) = ZERO
                end do
             case (3)
                do k = lo(3),hi(3)
                   do j = lo(2),hi(2)
                      if ( bc_dirichlet(mp(lo(1),j,k,1),1,0) ) up(lo(1),j,k,1:mgt%nc) = ZERO
                      if ( bc_dirichlet(mp(hi(1),j,k,1),1,0) ) up(hi(1),j,k,1:mgt%nc) = ZERO
                   end do
                end do

                do k = lo(3),hi(3)
                   do i = lo(1),hi(1)
                      if ( bc_dirichlet(mp(i,lo(2),k,1),1,0) ) up(i,lo(2),k,1:mgt%nc) = ZERO
                      if ( bc_dirichlet(mp(i,hi(2),k,1),1,0) ) up(i,hi(2),k,1:mgt%nc) = ZERO
                   end do
                end do

                do j = lo(2),hi(2)
                   do i = lo(1),hi(1)
                      if ( bc_dirichlet(mp(i,j,lo(3),1),1,0) ) up(i,j,lo(3),1:mgt%nc) = ZERO
                      if ( bc_dirichlet(mp(i,j,hi(3),1),1,0) ) up(i,j,hi(3),1:mgt%nc) = ZERO
                   end do
                end do
             end select

          end do
          !
          ! Nor do they preserve the value of shared nodes (ones at grid boundaries).
          !
          call multifab_internal_sync(uu)
       end if

    endif

    call destroy(bpt)
    call bl_proffortfuncstop("mg_tower_prolongation")

  end subroutine mg_tower_prolongation

  function mg_tower_converged(mgt, dd, Ynorm) result(r)
    use itsol_module
    logical                       :: r
    type(mg_tower), intent(inout) :: mgt
    real(dp_t),     intent(in   ) :: Ynorm
    type(multifab), intent(in   ) :: dd
    r = itsol_converged(dd, Ynorm, mgt%eps, mgt%abs_eps)
  end function mg_tower_converged

  subroutine mg_tower_cycle(mgt,cyc,lev,ss,uu,rh,mm,nu1,nu2,bottom_level,bottom_solve_time)

    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(inout) :: rh
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: ss
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: lev
    integer, intent(in) :: nu1, nu2
    integer, intent(in) :: cyc
    integer, intent(in), optional :: bottom_level
    real(dp_t), intent(inout), optional :: bottom_solve_time

    ! Note: this used to depend on mgt%cycle_type, but now we explicitly use the cycle_type
    !       that is passed in in "cyc" so that we can mix the cycle types in a single solve.
    select case ( cyc )
        case(MG_VCycle)
            call mg_tower_v_cycle(mgt,cyc,lev,ss,uu,rh,mm,nu1,nu2,1,bottom_level,bottom_solve_time)
        case(MG_WCycle)
            call mg_tower_v_cycle(mgt,cyc,lev,ss,uu,rh,mm,nu1,nu2,2,bottom_level,bottom_solve_time)
        case(MG_FCycle)
            call mg_tower_fmg_cycle(mgt,cyc,lev,ss,uu,rh,mm,nu1,nu2,bottom_level,bottom_solve_time)
         case default
            call bl_error('mg_tower_cycle: unknown cycle_type: ', mgt%cycle_type)
    end select

  end subroutine mg_tower_cycle

  recursive subroutine mg_tower_fmg_cycle(mgt, cyc, lev, ss, uu, rh, mm, nu1, nu2,&
                                 bottom_level, bottom_solve_time)
    use fabio_module
    use bl_prof_module

    type(mg_tower),  intent(inout) :: mgt
    type(multifab),  intent(inout) :: rh
    type(multifab),  intent(inout) :: uu
    type(multifab),  intent(in   ) :: ss
    type(imultifab), intent(in   ) :: mm
    integer,         intent(in   ) :: lev
    integer,         intent(in   ) :: nu1, nu2
    integer,         intent(in   ) :: cyc
    integer,         intent(in   ), optional :: bottom_level
    real(dp_t),      intent(inout), optional :: bottom_solve_time

    logical                   :: do_diag
    real(dp_t)                :: nrm, stime
    integer                   :: lbl
    character(len=3)          :: number
    character(len=20)         :: filename

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mgt_f_cycle")
    call bl_proffortfuncstart("mg_tower_fmg_cycle")

    lbl = 1; if ( present(bottom_level) ) lbl = bottom_level

    do_diag = .false.; if ( mgt%verbose >= 4 ) do_diag = .true.

    call setval(uu,ZERO,all=.true.)


    if ( lev == lbl ) then
       stime = parallel_wtime()

       if ( associated(mgt%bottom_mgt) ) then
          if ( do_diag ) then
             call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)
             nrm = norm_inf(mgt%cc(lev))
             if ( parallel_IOProcessor() ) &
                print *,'  DN: Norm before bottom         ',nrm
          end if

          call do_bottom_mgt(mgt, uu, rh)

          if ( do_diag ) then
             call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)
             nrm = norm_inf(mgt%cc(lev))
             if ( parallel_IOProcessor() ) &
                print *,'  UP: Norm after  bottom         ',nrm
          end if
       else
          call mg_tower_bottom_solve(mgt, lev, ss, uu, rh, mm)
       end if

       if ( present(bottom_solve_time) ) &
            bottom_solve_time = bottom_solve_time + (parallel_wtime()-stime)
    else
       call mg_tower_restriction(mgt, mgt%dd(lev-1), rh, mgt%mm(lev),mgt%mm(lev-1))
  
       call mg_tower_fmg_cycle(mgt, cyc, lev-1, mgt%ss(lev-1), mgt%uu(lev-1), &
                      mgt%dd(lev-1), mgt%mm(lev-1), nu1, nu2, bottom_level, bottom_solve_time)

!       call print(mgt%uu(lev-1),'uu before prolongation')

       if ( .false. ) then
          !
          ! Some debugging code I want to keep around for a while.
          ! I don't want to have to recreate this all the time :-)
          !
          write(number,fmt='(i3.3)') lev
          filename = 'uu_before_p' // number
          call fabio_write(mgt%uu(lev-1), 'debug', trim(filename))
       end if
          
       call mg_tower_prolongation(mgt, uu, lev-1)

       if (lev == mgt%nlevels) then
            call mg_tower_v_cycle(mgt, MG_VCycle, lev, mgt%ss(lev), uu, &
                              rh, mgt%mm(lev), nu1, nu2, 1, bottom_level, bottom_solve_time)
       else
            call mg_tower_v_cycle(mgt, MG_VCycle, lev, mgt%ss(lev), mgt%uu(lev), &
                              rh, mgt%mm(lev), nu1, nu2, 1, bottom_level, bottom_solve_time)
       end if
    end if

    call destroy(bpt)
    call bl_proffortfuncstop("mg_tower_fmg_cycle")
   
  end subroutine mg_tower_fmg_cycle

  recursive subroutine mg_tower_v_cycle(mgt, cyc, lev, ss, uu, rh, mm, nu1, nu2, gamma, &
                                        bottom_level, bottom_solve_time)
    use bl_prof_module

    type(mg_tower),  intent(inout) :: mgt
    type(multifab),  intent(inout) :: rh
    type(multifab),  intent(inout) :: uu
    type(multifab),  intent(in   ) :: ss
    type(imultifab), intent(in   ) :: mm
    integer,         intent(in   ) :: lev
    integer,         intent(in   ) :: nu1, nu2
    integer,         intent(in   ) :: gamma
    integer,         intent(in   ) :: cyc
    integer,         intent(in   ), optional :: bottom_level
    real(dp_t),      intent(inout), optional :: bottom_solve_time

    integer    :: i,lbl
    logical    :: do_diag
    real(dp_t) :: nrm, stime

    type(bl_prof_timer), save :: bpt

    !!! call bl_proffortfuncstart("mg_tower_v_cycle")
    call build(bpt, "mgt_v_cycle")

    lbl = 1; if ( present(bottom_level) ) lbl = bottom_level

    do_diag = .false.; if ( mgt%verbose >= 4 ) do_diag = .true.

    if ( parallel_IOProcessor() .and. do_diag) &
         write(6,1000) lev

    if ( get_dim(rh) == 1 ) then
       if ( do_diag ) then
          nrm = norm_inf(rh)
          if ( parallel_IOProcessor() ) &
             print *,'  DN: Norm before smooth         ',nrm
       end if

       call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)

       call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)

       if ( do_diag ) then
          nrm = norm_inf(mgt%cc(lev))
          if ( parallel_IOProcessor() ) &
              print *,'  DN: Norm after  smooth          ',nrm
       end if
    else if ( lev == lbl ) then
       stime = parallel_wtime()

       if ( associated(mgt%bottom_mgt) ) then
          if ( do_diag ) then
             call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)
             nrm = norm_inf(mgt%cc(lev))
             if ( parallel_IOProcessor() ) &
                print *,'  DN: Norm before bottom         ',nrm
          end if

          call do_bottom_mgt(mgt, uu, rh)

          if ( do_diag ) then
             call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)
             nrm = norm_inf(mgt%cc(lev))
             if ( parallel_IOProcessor() ) &
                print *,'  UP: Norm after  bottom         ',nrm
          end if

       else
          call mg_tower_bottom_solve(mgt, lev, ss, uu, rh, mm)
       end if

       if ( present(bottom_solve_time) ) &
            bottom_solve_time = bottom_solve_time + (parallel_wtime()-stime)
    else 
       if ( do_diag ) then
          if (cyc == MG_FCycle) then
              call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)
              nrm = norm_inf(mgt%cc(lev))
           else
              nrm = norm_inf(rh)
           end if
       end if

       if ( do_diag .and. parallel_IOProcessor() ) &
          print *,'  DN: Norm before smooth         ',nrm

       do i = 1, nu1
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do

       call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)

       if ( do_diag ) then
          nrm = norm_inf(mgt%cc(lev))
          if ( parallel_IOProcessor() ) &
              print *,'  DN: Norm after  smooth         ',nrm
       end if

       call mg_tower_restriction(mgt, mgt%dd(lev-1), mgt%cc(lev), &
                                 mgt%mm(lev),mgt%mm(lev-1))

       call setval(mgt%uu(lev-1), zero, all = .TRUE.)

       do i = gamma, 1, -1
          !!! call bl_proffortfuncstop("mg_tower_v_cycle")
          call mg_tower_v_cycle(mgt, cyc, lev-1, mgt%ss(lev-1), mgt%uu(lev-1), &
                              mgt%dd(lev-1), mgt%mm(lev-1), nu1, nu2, gamma, bottom_level, bottom_solve_time)
          !!! call bl_proffortfuncstart("mg_tower_v_cycle")
       end do

       call mg_tower_prolongation(mgt, uu, lev-1)

       if ( parallel_IOProcessor() .and. do_diag) &
          write(6,1000) lev

       if ( do_diag ) then
          call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)
          nrm = norm_inf(mgt%cc(lev))
          if ( parallel_IOProcessor() ) then
             print *,'  UP: Norm after  interp         ',nrm
          end if
       end if

       do i = 1, nu2
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do

       if ( do_diag ) then
          call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)
          nrm = norm_inf(mgt%cc(lev))
          if ( parallel_IOProcessor() ) &
             print *,'  UP: Norm after  smooth         ',nrm
       end if

    end if

    call destroy(bpt)
    !!! call bl_proffortfuncstop("mg_tower_v_cycle")

1000 format('AT LEVEL ',i2)

  end subroutine mg_tower_v_cycle

  subroutine mini_cycle(mgt, lev, ss, uu, rh, mm, nu1, nu2)

    type(mg_tower),  intent(inout) :: mgt
    type(multifab),  intent(in   ) :: rh
    type(multifab),  intent(inout) :: uu
    type(multifab),  intent(in   ) :: ss
    type(imultifab), intent(in   ) :: mm
    integer,         intent(in   ) :: lev
    integer,         intent(in   ) :: nu1, nu2

    integer    :: i
    logical    :: do_diag
    real(dp_t) :: nrm

    call bl_proffortfuncstart("mini_cycle")

    do_diag = .false.; if ( mgt%verbose >= 4 ) do_diag = .true.
    !
    ! Always relax first at the level we come in at.
    !
    if ( do_diag ) then
       nrm = norm_inf(rh)
       if ( parallel_IOProcessor() ) &
          print *,'MINI_FINE: Norm before smooth         ',nrm
    end if

    do i = 1, nu1
       call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
    end do

    call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)

    if ( do_diag ) then
       nrm = norm_inf(mgt%cc(lev))
       if ( parallel_IOProcessor() ) &
           print *,'MINI_FINE: Norm after  smooth         ',nrm
    end if
    !
    ! If we are doing a mini V-cycle here, then we must compute the coarse residual, 
    ! relax at the next lower level, then interpolate the correction and relax again.
    !
    if ( lev > 1 ) then

       call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)

       call mg_tower_restriction(mgt, mgt%dd(lev-1), mgt%cc(lev), mgt%mm(lev),mgt%mm(lev-1))

       call setval(mgt%uu(lev-1), zero, all = .TRUE.)

       if ( do_diag ) then
          nrm = norm_inf(mgt%dd(lev-1))
          if ( parallel_IOProcessor() ) &
             print *,'MINI_CRSE: Norm before smooth         ',nrm
       end if

       do i = 1, nu1
          call mg_tower_smoother(mgt, lev-1, mgt%ss(lev-1), mgt%uu(lev-1), mgt%dd(lev-1), mgt%mm(lev-1))
       end do

       if ( do_diag ) then
          call compute_defect(mgt%ss(lev-1), mgt%cc(lev-1), mgt%dd(lev-1), mgt%uu(lev-1), mgt%mm(lev-1), &
                         mgt%stencil_type, mgt%lcross, mgt%uniform_dh)
          nrm = norm_inf(mgt%cc(lev-1))
          if ( parallel_IOProcessor() ) &
             print *,'MINI_CRSE: Norm  after smooth         ',nrm
       end if

       call mg_tower_prolongation(mgt, uu, lev-1)

       if ( do_diag ) then
          call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)
          nrm = norm_inf(mgt%cc(lev))
          if ( parallel_IOProcessor() ) &
             print *,'MINI_FINE: Norm before smooth         ',nrm
       end if

       do i = 1, nu2
          call mg_tower_smoother(mgt, lev, ss, uu, rh, mm)
       end do

       if ( do_diag ) then
          call compute_defect(ss, mgt%cc(lev), rh, uu, mm, mgt%stencil_type, mgt%lcross, mgt%uniform_dh)
          nrm = norm_inf(mgt%cc(lev))
          if ( parallel_IOProcessor() ) &
             print *,'MINI_FINE: Norm before smooth         ',nrm
       end if

    end if

    call bl_proffortfuncstop("mini_cycle")

  end subroutine mini_cycle

  subroutine mg_tower_smoother(mgt, lev, ss, uu, ff, mm)

    use cc_mg_tower_smoother_module   , only:    cc_mg_tower_smoother 
    use nodal_mg_tower_smoother_module, only: nodal_mg_tower_smoother

    integer        , intent(in   ) :: lev
    type( mg_tower), intent(inout) :: mgt
    type( multifab), intent(inout) :: uu
    type( multifab), intent(in   ) :: ff
    type( multifab), intent(in   ) :: ss
    type(imultifab), intent(in   ) :: mm

    if (nodal_q(ff)) then 
        call nodal_mg_tower_smoother(mgt, lev, ss, uu, ff, mm)
    else
        call cc_mg_tower_smoother(mgt, lev, ss, uu, ff, mm)
    end if

  end subroutine mg_tower_smoother

end module mg_module
