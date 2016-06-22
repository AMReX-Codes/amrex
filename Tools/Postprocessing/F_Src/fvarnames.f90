! simply print out the list of variables stored in a plotfile
program fvarnames

  use plotfile_module
  use filler_module
  use bl_IO_module

  implicit none

  integer n, narg
  integer :: unit
  character (len=256) :: fname

  type(plotfile) :: pf

100 format (1x, i3, 3x, a)

  ! parse any formatting arguments
  narg = command_argument_count()

  if ( narg .eq. 0 ) then
     print *, ''
     print *, 'Usage:'
     print *, '     ftime plotfile'
     print *, ''
     print *, 'Description:'
     print *, '     This program takes a single plotfile and dumps out the list of variables'
     print *, ''
     stop
  endif

  unit = unit_new()

  call get_command_argument(1, value = fname)
  call build(pf, fname, unit)

  do n = 1, pf%nvars
     write(*,100) n, trim(pf%names(n))
  enddo

  call destroy(pf)

end program fvarnames
