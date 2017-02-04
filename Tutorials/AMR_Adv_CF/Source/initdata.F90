module initdata_module

  use my_amr_module
  use prob_module

  implicit none

  private

  public :: initdata

contains

  subroutine initdata ()
    if (len_trim(restart) .eq. 0) then
       call init_from_scratch();
    else
       call amrex_abort("init from checkpoint not implemented yet")
    end if
  end subroutine initdata

  subroutine init_from_scratch ()
    integer :: lev, finest, new_finest
    real(amrex_real) :: time
    type(amrex_boxarray), allocatable :: grids(:)
    type(amrex_distromap), allocatable :: dmaps(:)

    time = 0._amrex_real

    allocate(grids(0:amrex_max_level))
    allocate(dmaps(0:amrex_max_level))

    !
    ! Level 0
    !
    finest = 0
    call amrex_make_base_grids(grids(0))
    call amrex_distromap_build(dmaps(0), grids(0))

    call make_new_level(0, time, grids(0), dmaps(0))
    call init_level_data(0)

    !
    ! Fine levels
    !
    do while (finest < amrex_max_level)
!xxxxx       call amrex_make_new_grids(finest, time, new_finest, grids)
       if (new_finest .le. finest) exit
       finest = new_finest
       call amrex_distromap_build(dmaps(new_finest), grids(new_finest))

       call make_new_level(new_finest, time, grids(new_finest), dmaps(new_finest))
       call init_level_data(new_finest)
    end do

    do lev = 0, amrex_max_level
       call amrex_boxarray_destroy(grids(lev))
       call amrex_distromap_destroy(dmaps(lev))
    end do

  end subroutine init_from_scratch


  subroutine init_level_data (lev)
    integer, intent(in) :: lev

    type(amrex_geometry) :: geom
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:)

    geom = amrex_get_geometry(lev)

    call amrex_mfiter_build(mfi, phi_new(lev))

    do while (mfi%next())
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       call init_prob_data(lev, t_new(lev), bx%lo, bx%hi, phi, lbound(phi), ubound(phi), &
            geom%dx, geom%problo)
    end do

    call amrex_mfiter_destroy(mfi)

    call amrex_geometry_destroy(geom)
    
  end subroutine init_level_data

end module initdata_module
