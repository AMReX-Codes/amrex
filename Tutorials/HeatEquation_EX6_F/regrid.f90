module regrid_module

  use ml_layout_module
  use multifab_module
  use multifab_physbc_module
  use make_new_grids_module
  use ml_restriction_module
  use multifab_fill_ghost_module
  use bndry_reg_module

  implicit none

  private

  public :: regrid

contains

  subroutine regrid(mla,phi,flux,bndry_flx,nlevs,max_levs,dx,the_bc_tower,amr_buf_width,max_grid_size)

    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: phi(:)
    type(multifab) , intent(inout) :: flux(:,:)
    type(bndry_reg), intent(inout) :: bndry_flx(2:)
    integer        , intent(inout) :: nlevs, max_levs
    real(dp_t)     , intent(in   ) :: dx(:)
    type(bc_tower) , intent(inout) :: the_bc_tower
    integer        , intent(in   ) :: amr_buf_width, max_grid_size

    ! local variables
    type(layout)      :: la_array(max_levs)
    type(ml_layout)   :: mla_old
    type(ml_boxarray) :: mba
    type(multifab)    :: phi_orig(max_levs)

    integer :: i,dm,n,nl,n_buffer

    logical :: new_grid,properly_nested

    dm = mla%dim

    ! create a copy of the original mla
    call ml_layout_build(mla_old,mla%mba,mla%pmask)

    do n=1,nlevs
       ! create copies of the old data
       ! make sure to use mla_old since we will be destroying mla
       call multifab_build(phi_orig(n),mla_old%la(n),1,1)
       call multifab_copy_c(phi_orig(n),1,phi(n),1,1)

       ! get rid of the original multifab so we can create a new one
       ! with a new grid structure and the same name
       call multifab_destroy(phi(n))
    end do

    ! Destroy flux, phi_orig and bndry_reg before regridding -- we will create 
    !     new ones at the end
    do n = 1,nlevs
      do i = 1, dm
          call multifab_destroy(flux(n,i))
      end do
    end do

    do n = 2,nlevs
      call bndry_reg_destroy(bndry_flx(n))
    end do

    ! Get rid of the original mla so we can create a new one
    ! with a new grid structure and the same name
    call destroy(mla)

    ! mba is big enough to hold max_levs levels
    ! even though we know we had nlevs last time, we might 
    ! want more or fewer levels after regrid (if nlevs < max_levs)
    call ml_boxarray_build_n(mba,max_levs,dm)

    ! Tell mba about the ref_ratio between levels
    ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
    ! we use refinement ratio of 2 in every direction between all levels
    do n=2,max_levs
       mba%rr(n-1,:) = 2
    enddo
    
    ! copy the level 1 boxarray
    call copy(mba%bas(1),mla_old%mba%bas(1))

    ! set the problem domain at all levels
    mba%pd(1) = mla_old%mba%pd(1)
    do n=2,max_levs
       mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
    end do

    ! build the level 1 layout.
    call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),mla_old%pmask)

    ! This makes sure the boundary conditions are properly defined everywhere
    call bc_tower_level_build(the_bc_tower,1,la_array(1))

    ! build level 1 multifab
    call multifab_build(phi(1),la_array(1),1,1)
    
    ! copy level 1 data from original multifab
    call multifab_copy_c(phi(1),1,phi_orig(1),1,1)

    nl = 1
    new_grid = .true.

    ! this is the number of level n+1 buffer cells we require between levels
    ! n and n+2 for proper nesting
    n_buffer = 4

    do while ( (nl .lt. max_levs) .and. new_grid )

       ! need to fill ghost cells here in case we use them in tagging
       call multifab_fill_boundary(phi(nl))
       call multifab_physbc(phi(nl),1,1,1,the_bc_tower%bc_tower_array(nl))

       ! determine whether we need finer grids based on tagging criteria
       ! if so, return new_grid=T and the la_array(nl+1)
       call make_new_grids(new_grid,la_array(nl),la_array(nl+1),phi(nl),dx(nl), &
                           amr_buf_width,mba%rr(nl,1),nl,max_grid_size)

       if (new_grid) then

          ! set the level nl+1 boxarray
          call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

          ! enforce proper nesting within the grid creation procedure 
          if (nl .ge. 2) then

             ! Test on whether grids are already properly nested
             properly_nested = ml_boxarray_properly_nested(mba, n_buffer, mla_old%pmask, &
                                                           max_fine_level=nl+1)

             if (.not. properly_nested) then

                do n = 2,nl
                   ! Delete old multifabs so that we can rebuild them.
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
                   call multifab_build(phi(n),la_array(n),1,1)

                   ! first fill all refined cells by interpolating from coarse 
                   ! data underneath...
                   call fillpatch(phi(n),phi(n-1),phi(n)%ng,mba%rr(n-1,:), &
                                  the_bc_tower%bc_tower_array(n-1), &
                                  the_bc_tower%bc_tower_array(n), &
                                  1,1,1,1)

                   ! ... then overwrite with the original data at that level, if it existed
                   if (mla_old%nlevel .ge. n) then
                      call multifab_copy_c(phi(n),1,phi_orig(n),1,1)
                   end if

                end do

             end if ! if (.not. properly_nested)

          end if ! if (nl .ge. 2) then

          ! Define bc_tower at level nl+1.
          call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))

          ! Build the level nl+1 data only.
          call multifab_build(phi(nl+1),la_array(nl+1),1,1)

          ! first fill all refined cells by interpolating from coarse 
          ! data underneath...
          call fillpatch(phi(nl+1),phi(nl),phi(nl)%ng,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         1,1,1,1)

          ! ... then overwrite with the original data at that level, if it existed
          if (mla_old%nlevel .ge. nl+1) then
             call multifab_copy_c(phi(nl+1),1,phi_orig(nl+1),1,1)
          end if

          nl = nl+1
          nlevs = nl

       endif

    end do

    nlevs = nl

    ! Note: This build actually sets mla%la(n) = la_array(n) so we mustn't delete 
    !       la_array(n).  Doing the build this way means we don't have to re-create 
    !       all the multifabs because we have kept the same layouts.
    call build(mla,mba,la_array,mla_old%pmask,nlevs)

    ! This makes sure the boundary conditions are properly defined everywhere
    do n=1,nlevs
       call bc_tower_level_build(the_bc_tower,n,la_array(n))
    end do

    ! Create new flux, phi_orig and bndry_flx after regridding
    do n = 1,nlevs
      do i = 1, dm
          call multifab_build_edge(flux(n,i),mla%la(n),1,0,i)
      end do
    end do

    do n = 2,nlevs
      call bndry_reg_rr_build(bndry_flx(n),mla%la(n),mla%la(n-1),mla%mba%rr(n-1,:), &
                              ml_layout_get_pd(mla,n-1),other=.true.)
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(phi(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(phi(nlevs),1,1,1,the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(phi(n-1),phi(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(phi(n),phi(n-1),phi(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         1,1,1)

       enddo

    end if

    do n=1,mla_old%nlevel
       call destroy(phi_orig(n))
    end do

    call destroy(mba)
    call destroy(mla_old)

  end subroutine regrid

end module regrid_module
