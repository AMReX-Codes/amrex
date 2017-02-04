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
    type(amrex_boxarray) :: ba0
    type(amrex_distromap) :: dm0
    real(amrex_real) :: time

    time = 0._amrex_real

    !
    ! Level 0
    !
    call amrex_set_finest_level(0)

    call amrex_make_base_grids(ba0)
    call amrex_distromap_build(dm0, ba0)

    call make_new_level(0, time, ba0, dm0)
    
    call amrex_boxarray_destroy(ba0)
    call amrex_distromap_destroy(dm0)
 
    call init_level_data(0)

    call amrex_multifab_write(phi_new(0), "phi0")

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
