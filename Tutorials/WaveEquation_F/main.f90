program main

  ! These modules are all in BoxLib/Src/F_BaseLib
  use boxlib
  use multifab_module
  use bl_IO_module

  ! These modules are all in this directory
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
  real(dp_t) :: dx, dt, time, start_time, run_time, run_time_IOproc
  
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

  ! initialize the boxarray to be one single box
  call boxarray_build_bx(ba,bx)

  ! overwrite the boxarray to respect max_grid_size
  call boxarray_maxsize(ba,max_grid_size)

  ! build the layout, la
  ! the third argument is the problem domain, which in this case is bx
  call layout_build_ba(la,ba,bx,pmask=is_periodic)

  call destroy(ba)

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
  call destroy(la)

  deallocate(lo,hi,is_periodic,prob_lo,prob_hi)

  ! deallocate temporary boxarrays and communication mappings
  call layout_flush_copyassoc_cache()

  ! check for memory that should have been deallocated
  if ( parallel_IOProcessor() ) then
     print*, 'MEMORY STATS AT END OF PROGRAM'
     print*, ' '
  end if
  call print(multifab_mem_stats(),    "    multifab")
  call print(fab_mem_stats(),         "         fab")
  call print(boxarray_mem_stats(),    "    boxarray")
  call print(layout_mem_stats(),      "      layout")
  call print(boxassoc_mem_stats(),    "    boxassoc")
  call print(fgassoc_mem_stats(),     "     fgassoc")
  call print(syncassoc_mem_stats(),   "   syncassoc")
  call print(copyassoc_mem_stats(),   "   copyassoc")
  call print(fluxassoc_mem_stats(),   "   fluxassoc")

  ! parallel_wtime() returns the number of wallclock-time seconds since
  ! the program began
  run_time = parallel_wtime() - start_time

  ! collect run_time from each processor and store the maximum
  call parallel_reduce(run_time_IOproc, run_time, MPI_MAX, &
                       proc = parallel_IOProcessorNode())

  if ( parallel_IOProcessor() ) then
     print*,"Run time (s) =",run_time_IOproc
  end if

  call boxlib_finalize()

end program main
