module plotfile_module

  use amrex_amr_module
  use my_amr_module, only : plot_file, phi_new, t_new, stepno, ncomp
  
  use amr_data_module

  implicit none

  private

  public :: writeplotfile

contains

  subroutine writeplotfile ()

    type(amrex_multifab) :: plotmf(0:amrex_max_level)
    type(amrex_string), allocatable :: varname(:)
    integer, dimension(0:amrex_max_level) :: steps, rr
    integer :: ilev, nc
    character(len=127) :: name
    character(len=16)  :: current_step
    integer :: nlevs, icomp

    character(len=10), allocatable :: varnamedummy(:)
       
    nc = ncomp+6

    allocate(varname(nc))
    allocate(varnamedummy(ncomp))

    do icomp=1, ncomp
      write(varnamedummy(icomp),'("phi",I2.2)')icomp
    enddo

    do icomp=1, ncomp
       call amrex_string_build(varname(icomp), varnamedummy(icomp))
    enddo

    call amrex_string_build(varname(ncomp+1), "solution")
    call amrex_string_build(varname(ncomp+2), "rhs")
    call amrex_string_build(varname(ncomp+3), "exact_solution")
    call amrex_string_build(varname(ncomp+4), "error")
    call amrex_string_build(varname(ncomp+5), "acoef")
    call amrex_string_build(varname(ncomp+6), "bcoef")
   
    nlevs=amrex_get_finest_level()

      do ilev = 0, nlevs
          call amrex_multifab_build(plotmf(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, nc, 1)
   	do icomp=1,ncomp
          call plotmf(ilev) % copy(phi_new(ilev), icomp, icomp, 1, 0)
	enddo
        call plotmf(ilev) % copy(      solution(ilev), 1, ncomp+1, 1, 0)
        call plotmf(ilev) % copy(           rhs(ilev), 1, ncomp+2, 1, 0)
        call plotmf(ilev) % copy(exact_solution(ilev), 1, ncomp+3, 1, 0)
        call plotmf(ilev) % copy(      solution(ilev), 1, ncomp+4, 1, 0)
        call plotmf(ilev) % subtract(exact_solution(ilev),1,ncomp+4,1,0)
        call plotmf(ilev) % copy(acoef(ilev), 1, ncomp+5, 1, 0)
        call plotmf(ilev) % copy(bcoef(ilev), 1, ncomp+6, 1, 0)
    end do

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
	
    call amrex_write_plotfile(name, nlevs+1, plotmf, varname, amrex_geom, t_new(0), stepno, amrex_ref_ratio)

    ! let's not rely on finalizer, which is feature not all compilers support properly.
    do ilev = 0, nlevs
       call amrex_multifab_destroy(plotmf(ilev))
    end do
  end subroutine writeplotfile


end module plotfile_module
