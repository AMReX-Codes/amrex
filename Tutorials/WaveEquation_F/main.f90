program main

  use boxlib
  use parallel
  use multifab_module
  use bl_IO_module
  use ml_layout_module

  integer :: dim, nsteps, plot_int, n_cell, max_grid_size, verbose
  integer :: un, farg, narg
  integer, allocatable :: lo(:), hi(:)

  character(len=128) :: fname

  logical :: lexist, need_inputs

  type(ml_layout) :: mla

  type(box) :: bx

  type(multifab), allocatable :: state(:)

  namelist /probin/ dim, nsteps, plot_int, n_cell, max_grid_size, verbose

  ! if running in parallel, this will print out the number of MPI 
  ! processes and OpenMP threads
  call boxlib_initialize()

  ! default values
  dim           = 2
  nsteps        = 10
  plot_int      = 1
  n_cell        = 64
  max_grid_size = 32
  verbose       = 0

  ! read inputs file and overwrite any default values
  narg = command_argument_count()
  need_inputs = .true.
  farg = 1
  if ( need_inputs .AND. narg >= 1 ) then
     call get_command_argument(farg, value = fname)
     inquire(file = fname, exist = lexist )
     if ( lexist ) then
        farg = farg + 1
        un = unit_new()
        open(unit=un, file = fname, status = 'old', action = 'read')
        read(unit=un, nml = probin)
        close(unit=un)
        need_inputs = .false.
     end if
  end if

  allocate(lo(dim),hi(dim))

  lo(:) = 0
  hi(:) = n_cell-1






  allocate(state(1))






  deallocate(lo,hi)
  deallocate(state)

  call boxlib_finalize()

end program main
