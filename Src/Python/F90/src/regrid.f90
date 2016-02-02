module regrid_module

  use iso_c_binding
  use ml_layout_module
  use multifab_module
  use make_new_grids_module
  use ml_restrict_fill_module
  use tag_boxes_module, only : tagging_needs_ghost_cells

  implicit none

  private

  public :: regrid

contains

  subroutine regrid(mla,phi,nlevs,max_levs,dx,the_bc_tower,amr_buf_width,max_grid_size,tag_boxes_cb)

    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: phi(:)
    integer        , intent(inout) :: nlevs
    real(dp_t)     , intent(in   ) :: dx(:)
    type(bc_tower) , intent(inout) :: the_bc_tower
    integer        , intent(in   ) :: amr_buf_width, max_grid_size, max_levs
    type(c_funptr) , intent(in   ) :: tag_boxes_cb

    ! local variables
    type(layout)      :: la_array(max_levs)
    type(multifab)    :: phi_orig(max_levs)
    type(multifab), allocatable :: phi_opt(:)
    type(ml_boxarray) :: mba

    integer :: dm,n,nl,n_buffer, nlevs_old, nc, ng

    logical :: new_grid, properly_nested, pmask(mla%dim), same_boxarray

    nlevs_old = nlevs

    dm = mla%dim
    pmask = mla%pmask

    nc = ncomp(phi(1))
    ng = nghost(phi(1))

    la_array(1) = mla%la(1)

    phi_orig(1) = phi(1)

    do n=2,nlevs
       ! create copies of the old data
       call multifab_build(phi_orig(n),mla%la(n),nc,ng)
       call multifab_copy(phi_orig(n),phi(n))

       ! get rid of the original multifab so we can create a new one
       ! with a new grid structure and the same name
       call multifab_destroy(phi(n))
    end do

    ! mba is big enough to hold max_levs levels
    ! even though we know we had nlevs last time, we might
    ! want more or fewer levels after regrid (if nlevs < max_levs)
    call ml_boxarray_build_n(mba,max_levs,dm)
    call ml_boxarray_set_ref_ratio(mba)

    ! copy the level 1 boxarray
    call copy(mba%bas(1),mla%mba%bas(1))

    ! set the problem domain at all levels
    mba%pd(1) = mla%mba%pd(1)
    do n=2,max_levs
       mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
    end do

    nl = 1
    new_grid = .true.

    ! this is the number of level n+1 buffer cells we require between levels
    ! n and n+2 for proper nesting
    n_buffer = 4

    do while ( (nl .lt. max_levs) .and. new_grid )

       if (tagging_needs_ghost_cells) then
          call multifab_fill_boundary(phi(nl))
          call multifab_physbc(phi(nl),1,1,nc,the_bc_tower%bc_tower_array(nl))
       end if

       ! determine whether we need finer grids based on tagging criteria
       ! if so, return new_grid=T and the la_array(nl+1)
       call make_new_grids(new_grid,la_array(nl),la_array(nl+1),phi(nl),dx(nl), &
                           amr_buf_width,mba%rr(nl,1),nl,max_grid_size,tag_boxes_cb)

       if (new_grid) then

          ! set the level nl+1 boxarray
          call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

          ! enforce proper nesting within the grid creation procedure
          if (nl .ge. 2) then

             ! Test on whether grids are already properly nested
             properly_nested = ml_boxarray_properly_nested(mba, n_buffer, pmask, &
                                                           max_fine_level=nl+1)

             if (.not. properly_nested) then

                do n = 2,nl
                   ! Delete multifabs so that we can rebuild them.
                   call destroy(phi(n))
                end do

                ! Change the layout at levels 2 through nl so the new grid
                ! structure is properly nested
                call enforce_proper_nesting(mba,la_array,max_grid_size)

                ! Loop over all the lower levels which we might have changed
                ! when we enforced proper nesting.
                do n = 2,nl

                   ! This makes sure the boundary conditions are properly defined everywhere
                   call bc_tower_level_build(the_bc_tower,n,la_array(n))

                   ! Rebuild the lower level data again if it changed.
                   call multifab_build(phi(n),la_array(n),nc,ng)

                   if (mla%nlevel .ge. n) then
                      same_boxarray = boxarray_same_q(get_boxarray(phi     (n)), &
                           &                          get_boxarray(phi_orig(n)))
                   else
                      same_boxarray = .false.
                   end if

                   if (.not. same_boxarray) then
                      ! first fill all refined cells by interpolating from coarse
                      ! data underneath...  no need to fill ghost cells
                      call fillpatch(phi(n),phi(n-1),0,mba%rr(n-1,:), &
                           the_bc_tower%bc_tower_array(n-1), &
                           the_bc_tower%bc_tower_array(n), &
                           1,1,1,nc,no_final_physbc_input=.true.)
                   end if

                   ! ... then overwrite with the original data at that level, if it existed
                   if (mla%nlevel .ge. n) then
                      call multifab_copy(phi(n),phi_orig(n))
                   end if

                end do

             end if ! if (.not. properly_nested)

          end if ! if (nl .ge. 2) then

          ! Define bc_tower at level nl+1.
          call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))

          ! Build the level nl+1 data only.
          call multifab_build(phi(nl+1),la_array(nl+1),nc,ng)

          if (mla%nlevel .ge. nl+1) then
             same_boxarray = boxarray_same_q(get_boxarray(phi     (nl+1)), &
                  &                          get_boxarray(phi_orig(nl+1)))
          else
             same_boxarray = .false.
          end if

          if (.not.same_boxarray) then
             ! first fill all refined cells by interpolating from coarse
             ! data underneath...
             ! no need to fill ghost cells
             call fillpatch(phi(nl+1),phi(nl),0,mba%rr(nl,:), &
                  the_bc_tower%bc_tower_array(nl), &
                  the_bc_tower%bc_tower_array(nl+1), &
                  1,1,1,nc,no_final_physbc_input=.true.)
          end if

          ! ... then overwrite with the original data at that level, if it existed
          if (mla%nlevel .ge. nl+1) then
             call multifab_copy(phi(nl+1),phi_orig(nl+1))
          end if

          nl = nl+1
          nlevs = nl

       end if

    end do  ! end do while ( (nl .lt. max_levs) .and. new_grid )

    nlevs = nl

    do n=2, nlevs_old
       call multifab_destroy(phi_orig(n))
    end do

    ! this destroys everything in mla except the coarsest layout
    call destroy(mla, keep_coarse_layout=.true.)

    call ml_layout_build_la_array(mla,la_array,mba,pmask,nlevs)
    call destroy(mba)

    ! We need to move data if a layout in la_array is not used in mla.
    ! We also need to destroy any unused layouts.
    allocate(phi_opt(nlevs))
    do n=1,nlevs
       if (mla%la(n) .ne. la_array(n)) then
          call multifab_build(phi_opt(n), mla%la(n), nc, ng)
          call multifab_copy(phi_opt(n), phi(n))
          call destroy(phi(n))
          call destroy(la_array(n))
          phi(n) = phi_opt(n)
       end if
    end do
    deallocate(phi_opt)

    ! This makes sure the boundary conditions are properly defined everywhere
    do n=1,nlevs
       call bc_tower_level_build(the_bc_tower,n,mla%la(n))
    end do

    call ml_restrict_and_fill(nlevs, phi, mla%mba%rr, the_bc_tower%bc_tower_array)

  end subroutine regrid

end module regrid_module
