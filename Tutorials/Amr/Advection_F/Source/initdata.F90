module initdata_module

  use amrex_amr_module

  use my_amr_module, only : restart, plot_int
  use plotfile_module, only : writeplotfile
  use averagedown_module, only : averagedown

  implicit none

  private

  public :: initdata

contains

  subroutine initdata ()
    if (len_trim(restart) .eq. 0) then
       call amrex_init_from_scratch(0.0_amrex_real)
       call averagedown()

       call init_particle_data()
       
       if (plot_int .gt. 0) call writeplotfile
    else
       call amrex_abort("init from checkpoint not implemented yet")
    end if
  end subroutine initdata

  subroutine init_particle_data()
    use amr_data_module, only : pc, phi_new    
    use amrex_particlecontainer_module, only : amrex_particlecontainer_build
    use prob_module, only : init_part_data
    
    type(c_ptr)          :: amrcore = c_null_ptr    
    integer              :: lev
    type(amrex_mfiter)   :: mfi
    type(amrex_box)      :: bx
    
    amrcore = amrex_get_amrcore()
    call amrex_particlecontainer_build(pc, amrcore)       

    lev = 0
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.false.)
    do while(mfi%next())
       bx = mfi%tilebox()
       call init_part_data(pc, lev, mfi, bx%lo, bx%hi, amrex_geom(lev)%dx, amrex_problo)
    end do
    
    call pc%redistribute()

  end subroutine init_particle_data

end module initdata_module
