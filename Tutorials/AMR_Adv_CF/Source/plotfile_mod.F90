module plotfile_module

  use amrex_amr_module
  use my_amr_module, only : plot_file, phi_new, t_new, istep

  implicit none

  private

  public :: writeplotfile

contains

  subroutine writeplotfile ()

    integer :: nlevs
    character(len=127) :: name
    character(len=16)  :: current_step
    type(amrex_string) :: varname(1)
    
    if      (istep(0) .lt. 1000000) then
       write(current_step,fmt='(i5.5)') istep(0)
    else if (istep(0) .lt. 10000000) then
       write(current_step,fmt='(i6.6)') istep(0)
    else if (istep(0) .lt. 100000000) then
       write(current_step,fmt='(i7.7)') istep(0)
    else if (istep(0) .lt. 1000000000) then
       write(current_step,fmt='(i8.8)') istep(0)
    else
       write(current_step,fmt='(i15.15)') istep(0)
    end if
    name = trim(plot_file) // current_step

    nlevs = amrex_get_numlevels()

    call amrex_string_build(varname(1), "phi")

    call amrex_write_plotfile(name, nlevs, phi_new, varname, amrex_geom, t_new(0), istep, amrex_ref_ratio)

  end subroutine writeplotfile

end module plotfile_module

