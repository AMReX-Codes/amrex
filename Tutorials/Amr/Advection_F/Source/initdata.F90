module initdata_module

  use amrex_amr_module

  use amrex_amrtracerparticlecontainer_module, only: amrex_amrtracerparticlecontainer_init, &
       amrex_init_particles_one_per_cell
  use my_amr_module, only : restart, plot_int
  use plotfile_module, only : writeplotfile
  use averagedown_module, only : averagedown

  implicit none

  private

  public :: initdata

contains

  subroutine initdata ()
    type(c_ptr) :: amrcore = c_null_ptr
    if (len_trim(restart) .eq. 0) then
       call amrex_init_from_scratch(0.0_amrex_real)
       call averagedown()

       amrcore = amrex_get_amrcore()
       call amrex_amrtracerparticlecontainer_init(amrcore)
       call amrex_init_particles_one_per_cell()

       if (plot_int .gt. 0) call writeplotfile
    else
       call amrex_abort("init from checkpoint not implemented yet")
    end if
  end subroutine initdata

end module initdata_module
