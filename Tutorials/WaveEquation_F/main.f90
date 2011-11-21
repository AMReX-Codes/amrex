program main

  use boxlib
  use parallel
  use multifab_module
  use bl_IO_module
  use layout_module
  use init_data_module
  use write_plotfile_module
  use advance_module

  implicit none

  ! stuff you can set with the inputs file (otherwise use default values below)
  integer :: dim, nsteps, plot_int, n_cell, max_grid_size, verbose

  ! dummy indices using for reading in inputs file
  integer :: un, farg, narg
  logical :: need_inputs_file, found_inputs_file
  character(len=128) :: inputs_file_name

  integer, allocatable :: lo(:), hi(:)
  integer :: istep

  real(dp_t), allocatable :: prob_lo(:), prob_hi(:)
  real(dp_t) :: dx, dt, time, start_time, end_time
  
  logical, allocatable :: is_periodic(:)

  type(box)      :: bx
  type(boxarray) :: ba
  type(layout)   :: la
  type(multifab) :: data

  namelist /probin/ dim, nsteps, plot_int, n_cell, max_grid_size, verbose

  ! if running in parallel, this will print out the number of MPI 
  ! processes and OpenMP threads
  call boxlib_initialize()

  ! parallel_wtime() returns the number of wallclock-time seconds since
  ! the program began
  start_time = parallel_wtime()

  ! default values - will get overwritten by the inputs file
  dim           = 2
  nsteps        = 10
  plot_int      = 1
  n_cell        = 64
  max_grid_size = 32
  verbose       = 0

  ! read inputs file and overwrite any default values
  narg = command_argument_count()
  need_inputs_file = .true.
  farg = 1
  if ( need_inputs_file .AND. narg >= 1 ) then
     call get_command_argument(farg, value = inputs_file_name)
     inquire(file = inputs_file_name, exist = found_inputs_file )
     if ( found_inputs_file ) then
        farg = farg + 1
        un = unit_new()
        open(unit=un, file = inputs_file_name, status = 'old', action = 'read')
        read(unit=un, nml = probin)
        close(unit=un)
        need_inputs_file = .false.
     end if
  end if

  ! now that we have dim, we can allocate these
  allocate(lo(dim),hi(dim))
  allocate(is_periodic(dim))
  allocate(prob_lo(dim),prob_hi(dim))

  ! physical problem is a box on (-1,-1) to (1,1), periodic on all sides
  prob_lo(:) = -1.d0
  prob_hi(:) =  1.d0
  is_periodic(:) = .true.

  ! create a box from (0,0) to (n_cell-1,n_cell-1)
  lo(:) = 0
  hi(:) = n_cell-1
  bx = make_box(lo,hi)

  ! the grid spacing is the same in each direction
  dx = (prob_hi(1)-prob_lo(1)) / n_cell

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initialize the boxarray to be one single box
  call build(ba,bx)
  ! overwrite the boxarray to respect max_grid_size
  call boxarray_maxsize(ba,max_grid_size)

  ! build the layout, la
  call build(la,ba,pmask=is_periodic)

  ! build multifab with 2 components and 6 ghost cells
  call multifab_build(data,la,2,6)
  
  ! initialze data
  call init_data(data,dx,prob_lo)

  istep = 0
  time = 0.d0

  dt = 0.1d0*dx

  ! write out plotfile 0
  call write_plotfile(la,data,istep,dx,time,prob_lo,prob_hi)

  do istep=1,nsteps

     if (verbose .ge. 1) then
        ! we only want one processor to write to screen
        if ( parallel_IOProcessor() ) then
           print*,'Advancing time step',istep
        end if
     end if
     
     ! advance the data
     call advance(data,dx,dt)

     time = time + dt

     if (mod(istep,plot_int) .eq. 0 .or. istep .eq. nsteps) then
        ! write out plotfile
        call write_plotfile(la,data,istep,dx,time,prob_lo,prob_hi)
     end if

  end do

  ! make sure to destroy the multifab or you'll leak memory
  call destroy(data)

  deallocate(lo,hi,is_periodic,prob_lo,prob_hi)

  ! parallel_wtime() returns the number of wallclock-time seconds since
  ! the program began
  end_time = parallel_wtime()

  call boxlib_finalize()

  if ( parallel_IOProcessor() ) then
     print*,"Run time (s) =",end_time-start_time
  end if

end program main
