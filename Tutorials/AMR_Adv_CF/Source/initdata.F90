module initdata_module

  use my_amr_module
  use prob_module

  implicit none

  private

  public :: initdata

contains

  subroutine initdata ()
    if (len_trim(restart) .eq. 0) then
       call amrex_init_from_scratch(0.0_amrex_real)
!xxx       call average_down()
    else
       call amrex_abort("init from checkpoint not implemented yet")
    end if
  end subroutine initdata

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
