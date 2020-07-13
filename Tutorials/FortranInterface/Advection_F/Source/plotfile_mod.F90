module plotfile_module

  use amrex_amr_module
  use amr_data_module, only : pc
  use my_amr_module, only : plot_file, phi_new, t_new, stepno

  implicit none

  private

  public :: writeplotfile

contains

  subroutine writeplotfile ()

    integer :: nlevs
    character(len=127) :: name
    character(len=16)  :: current_step
    type(amrex_string) :: varname(1)
    
    if      (stepno(0) .lt. 1000000) then
       write(current_step,fmt='(i5.5)') stepno(0)
    else if (stepno(0) .lt. 10000000) then
       write(current_step,fmt='(i6.6)') stepno(0)
    else if (stepno(0) .lt. 100000000) then
       write(current_step,fmt='(i7.7)') stepno(0)
    else if (stepno(0) .lt. 1000000000) then
       write(current_step,fmt='(i8.8)') stepno(0)
    else
       write(current_step,fmt='(i15.15)') stepno(0)
    end if
    name = trim(plot_file) // current_step

    nlevs = amrex_get_numlevels()

    call amrex_string_build(varname(1), "phi")

    call amrex_write_plotfile(name, nlevs, phi_new, varname, amrex_geom, &
         t_new(0), stepno, amrex_ref_ratio)

    call pc%write(name, "Tracer", .true.)
    
  end subroutine writeplotfile

end module plotfile_module

