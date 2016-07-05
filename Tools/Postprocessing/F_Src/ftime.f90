! simply print out the simulation time for a list of plotfiles.
program ftime

  use plotfile_module
  use filler_module
  use bl_IO_module

  implicit none

  integer f
  integer :: unit

  type(plotfile) :: pf
  character (len=256) :: fname, format_string

  integer :: narg, farg, ichomp

  ! default format
  format_string = '(A, 1x, g15.7)'

  ! parse any formatting arguments
  narg = command_argument_count()

  if ( narg .eq. 0 ) then
     print *, ''
     print *, 'Usage:'
     print *, '     ftime [ -f | --format <format_string> ] <pltfile_list>'
     print *, ''
     print *, 'Description:'
     print *, '     This program takes a whitespace-separated list of plotfiles and'
     print *, '     returns the time for each plotfile.'
     print *, ''
     print *, 'Options:'
     print *, '     -f | --format:    This allows the user to set a FORTRAN format statement'
     print *, '                       for the output of the plotfile name and time.'
     print *, ' '
     print *, '                       <format_string> should contain both a string and'
     print *, '                       floating point format specifier, in a standard'
     print *, '                       FORTRAN format such as "(A, 1x, ES14.7)" (without'
     print *, '                       the quotes).  The default is set to "(A, 1x, g15.7)".'
     print *, ''
     stop
  endif


  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-f', '--format')
        farg = farg + 1
        call get_command_argument(farg, value = fname)

        ichomp=index(trim(format_string), " ", back=.true.)
        format_string = format_string(:ichomp) // trim(fname) // ')'

     case default
        exit

     end select

     farg = farg + 1
  end do

  unit = unit_new()

  do f = farg, narg

     call get_command_argument(f, value = fname)
     call build(pf, fname, unit)

     write(*,format_string) trim(fname), pf%tm

     call destroy(pf)

  enddo

end program ftime
