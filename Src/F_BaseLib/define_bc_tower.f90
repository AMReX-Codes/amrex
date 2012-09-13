module define_bc_module

  use bl_types
  use ml_layout_module
  use bc_module

  implicit none

  type bc_level

     ! 1st index is the grid number (grid "0" corresponds to the entire problem domain)
     ! 2nd index is the direction (1=x, 2=y, 3=z)
     ! 3rd index is the side (1=lo, 2=hi)
     ! 4th index is the variable (only assuming 1 variable here)
     integer, pointer :: phys_bc_level_array(:,:,:) => Null()
     integer, pointer ::  adv_bc_level_array(:,:,:,:) => Null()
     integer, pointer ::  ell_bc_level_array(:,:,:,:) => Null()

  end type bc_level

  type bc_tower

     integer :: max_level_built = 0

     ! an array of bc_levels, one for each level of refinement
     type(bc_level), pointer :: bc_tower_array(:) => Null()

     ! 1st index is the direction (1=x, 2=y, 3=z)
     ! 2nd index is the side (1=lo, 2=hi)
     integer       , pointer :: domain_bc(:,:) => Null()

  end type bc_tower

  private

  interface build
     module procedure bc_tower_init
     module procedure bc_tower_level_build
  end interface build

  interface destroy
     module procedure bc_tower_destroy
  end interface destroy

  public :: bc_level, bc_tower, bc_tower_init, bc_tower_level_build, bc_tower_destroy, build, destroy

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine bc_tower_init(bct,num_levs,dm,phys_bc_in)

    type(bc_tower ), intent(  out) :: bct
    integer        , intent(in   ) :: num_levs
    integer        , intent(in   ) :: dm
    integer        , intent(in   ) :: phys_bc_in(:,:)

    allocate(bct%bc_tower_array(num_levs))
    allocate(bct%domain_bc(dm,2))

    bct%domain_bc(:,:) = phys_bc_in(:,:)

  end subroutine bc_tower_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine bc_tower_level_build(bct,n,la)

    type(bc_tower ), intent(inout) :: bct
    integer        , intent(in   ) :: n
    type(layout)   , intent(in   ) :: la

    integer :: ngrids,ncomp,dm

    ncomp = 1
    dm = layout_dim(la)

    if (associated(bct%bc_tower_array(n)%phys_bc_level_array)) then
      deallocate(bct%bc_tower_array(n)%phys_bc_level_array)
      deallocate(bct%bc_tower_array(n)%adv_bc_level_array)
      deallocate(bct%bc_tower_array(n)%ell_bc_level_array)
      bct%bc_tower_array(n)%phys_bc_level_array => NULL()
      bct%bc_tower_array(n)%adv_bc_level_array => NULL()
      bct%bc_tower_array(n)%ell_bc_level_array => NULL()
    end if

    ngrids = layout_nlocal(la)

    allocate(bct%bc_tower_array(n)%phys_bc_level_array(0:ngrids,dm,2))
    call phys_bc_level_build(bct%bc_tower_array(n)%phys_bc_level_array,la, &
                             bct%domain_bc)

    ! Here we allocate 1 component and set the default to be INTERIOR
    allocate(bct%bc_tower_array(n)%adv_bc_level_array(0:ngrids,dm,2,ncomp))
    call adv_bc_level_build(bct%bc_tower_array(n)%adv_bc_level_array, &
                            bct%bc_tower_array(n)%phys_bc_level_array)

    ! Here we allocate 1 component and set the default to be BC_INT
    allocate(bct%bc_tower_array(n)%ell_bc_level_array(0:ngrids,dm,2,ncomp))
    call ell_bc_level_build(bct%bc_tower_array(n)%ell_bc_level_array, &
                            bct%bc_tower_array(n)%phys_bc_level_array)

     bct%max_level_built = n

  end subroutine bc_tower_level_build

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine phys_bc_level_build(phys_bc_level,la_level,domain_bc)

    integer     , intent(inout) :: phys_bc_level(0:,:,:)
    integer     , intent(in   ) :: domain_bc(:,:)
    type(layout), intent(in   ) :: la_level

    ! local
    type(box) :: bx,pd
    integer :: d,i

    pd = layout_get_pd(la_level) 

    ! set to interior everywhere, then overwrite with physical bc's where appropriate
    phys_bc_level = INTERIOR

    ! grid 0 corresponds to the problem domain
    do d = 1,layout_dim(la_level)
       phys_bc_level(0,d,1) = domain_bc(d,1)
       phys_bc_level(0,d,2) = domain_bc(d,2)
    end do

    ! loop over individual grids
    do i = 1,layout_nlocal(la_level)    ! loop over grids
       bx = layout_get_box(la_level,global_index(la_level,i))  ! grab box associated with the grid
       do d = 1,layout_dim(la_level)    ! loop over directions
          ! if one side of a grid is a domain boundary, set the 
          ! physical boundary condition
          if (lwb(bx,d) == lwb(pd,d)) phys_bc_level(i,d,1) = domain_bc(d,1)
          if (upb(bx,d) == upb(pd,d)) phys_bc_level(i,d,2) = domain_bc(d,2)
       end do
    end do

  end subroutine phys_bc_level_build

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine adv_bc_level_build(adv_bc_level,phys_bc_level)

    integer  , intent(inout) ::  adv_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)

    integer :: dm
    integer :: igrid,d,lohi

    adv_bc_level = INTERIOR

    dm = size(adv_bc_level,dim=2)

    ! if the physical boundary conditions is something other than
    ! INTERIOR or PERIODIC, then overwrite the default value
    do igrid=0,size(adv_bc_level,dim=1)-1   ! loop over grids
    do d=1,dm                               ! loop over directions
    do lohi=1,2                             ! loop over lo/hi side

       if (phys_bc_level(igrid,d,lohi) == INLET) then

          call bl_error("define_bc_tower.f90: INLET not supported for this example")

       else if (phys_bc_level(igrid,d,lohi) == OUTLET) then

          adv_bc_level(igrid,d,lohi,1) = FOEXTRAP

       else if (phys_bc_level(igrid,d,lohi) == SYMMETRY) then

          call bl_error("define_bc_tower.f90: SYMMETRY not supported for this example")

       else if (phys_bc_level(igrid,d,lohi) == SLIP_WALL) then

          call bl_error("define_bc_tower.f90: SLIP_WALL not supported for this example")

       else if (phys_bc_level(igrid,d,lohi) == NO_SLIP_WALL) then

          adv_bc_level(igrid,d,lohi,1) = EXT_DIR

       end if

    end do
    end do
    end do

  end subroutine adv_bc_level_build

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ell_bc_level_build(ell_bc_level,phys_bc_level)

    integer  , intent(inout) ::  ell_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)

    integer :: dm
    integer :: igrid,d,lohi

    ell_bc_level = BC_INT
 
    dm = size(ell_bc_level,dim=2)

    ! if the physical boundary conditions is something other than
    ! INTERIOR, then overwrite the default value
    do igrid=0,size(ell_bc_level,dim=1)-1   ! loop over grids
    do d=1,dm                               ! loop over directions
    do lohi=1,2                             ! loop over lo/hi side

       if (phys_bc_level(igrid,d,lohi) == INLET) then

          call bl_error("define_bc_tower.f90: INLET not supported for this example")

       else if (phys_bc_level(igrid,d,lohi) == OUTLET) then

          ell_bc_level(igrid,d,lohi,1) = BC_NEU

       else if (phys_bc_level(igrid,d,lohi) == SYMMETRY) then

          call bl_error("define_bc_tower.f90: SYMMETRY not supported for this example")

       else if (phys_bc_level(igrid,d,lohi) == SLIP_WALL) then

          call bl_error("define_bc_tower.f90: SLIP_WALL not supported for this example")

       else if (phys_bc_level(igrid,d,lohi) == NO_SLIP_WALL) then

          ell_bc_level(igrid,d,lohi,1) = BC_DIR

       else if (phys_bc_level(igrid,d,lohi) == PERIODIC) then

          ell_bc_level(igrid,d,lohi,1) = BC_PER

       end if

    end do
    end do
    end do

  end subroutine ell_bc_level_build

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine bc_tower_destroy(bct)

    type(bc_tower), intent(inout) :: bct

    integer :: n

    do n = 1,bct%max_level_built
       deallocate(bct%bc_tower_array(n)%phys_bc_level_array)
       deallocate(bct%bc_tower_array(n)%adv_bc_level_array)
       deallocate(bct%bc_tower_array(n)%ell_bc_level_array)
       bct%bc_tower_array(n)%phys_bc_level_array => NULL()
       bct%bc_tower_array(n)%adv_bc_level_array => NULL()
       bct%bc_tower_array(n)%ell_bc_level_array => NULL()
    end do
    deallocate(bct%bc_tower_array)

    deallocate(bct%domain_bc)

  end subroutine bc_tower_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module define_bc_module
