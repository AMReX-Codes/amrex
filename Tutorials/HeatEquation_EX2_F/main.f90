program main

  use boxlib
  use multifab_module
  use bl_IO_module
  use layout_module
  use init_phi_module
  use write_plotfile_module
  use advance_module
  use define_bc_module

  implicit none

  ! stuff you can set with the inputs file (otherwise use default values below)
  integer :: dim, nsteps, plot_int, n_cell, max_grid_size
  integer :: bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi, bc_z_lo, bc_z_hi

  ! dummy indices using for reading in inputs file
  integer :: un, farg, narg
  logical :: need_inputs_file, found_inputs_file
  character(len=128) :: inputs_file_name

  integer, allocatable :: phys_bc(:,:)
  integer, allocatable :: lo(:), hi(:)
  integer :: istep,i

  real(dp_t), allocatable :: prob_lo(:), prob_hi(:)
  real(dp_t) :: dx, dt, time, start_time, run_time, run_time_IOproc
  
  logical, allocatable :: is_periodic(:)

  type(box)      :: bx
  type(boxarray) :: ba
  type(layout)   :: la
  type(multifab) :: phi

  type(bc_tower) :: the_bc_tower

  namelist /probin/ dim, nsteps, plot_int, n_cell, max_grid_size, &
       bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi, bc_z_lo, bc_z_hi

  ! if running in parallel, this will print out the number of MPI 
  ! processes and OpenMP threads
  call boxlib_initialize()

  ! parallel_wtime() returns the number of wallclock-time seconds since
  ! the program began
  start_time = parallel_wtime()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! default values - will get overwritten by the inputs file

  dim           = 2
  nsteps        = 1000
  plot_int      = 100
  n_cell        = 256
  max_grid_size = 64


  ! allowable options for this example are
  ! -1 = PERIODIC
  ! 11 = INLET        (phi = val    at boundary)
  ! 12 = OUTLET       (phi = extrap at boundary)
  ! 14 = SLIP_WALL    (dphi/dn=0    at boundary)
  ! 15 = NO_SLIP_WALL (dphi/dn=0    at boundary)
  bc_x_lo       = -1 ! PERIODIC
  bc_x_hi       = -1 ! PERIODIC
  bc_y_lo       = -1 ! PERIODIC
  bc_y_hi       = -1 ! PERIODIC
  bc_z_lo       = -1 ! PERIODIC
  bc_z_hi       = -1 ! PERIODIC

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  allocate(phys_bc(dim,2))
  allocate(is_periodic(dim))
  allocate(prob_lo(dim),prob_hi(dim))

  ! physical problem is a box on (-1,-1) to (1,1)
  prob_lo(:) = -1.d0
  prob_hi(:) =  1.d0

  ! build a single array to hold domain boundary conditions
  phys_bc(1,1) = bc_x_lo
  phys_bc(1,2) = bc_x_hi
  phys_bc(2,1) = bc_y_lo
  phys_bc(2,2) = bc_y_hi
  if (dim .eq. 3) then
     phys_bc(3,1) = bc_z_lo
     phys_bc(3,2) = bc_z_hi
  end if

  ! build an array indicating periodicity in each direction
  is_periodic(:) = .false.
  do i=1,dim
     if (phys_bc(i,1) .eq. -1 .and. phys_bc(i,2) .ne. -1) then
        call bl_error("Invalid BC's - both lo and hi need to be periodic")
     end if
     if (phys_bc(i,2) .eq. -1 .and. phys_bc(i,1) .ne. -1) then
        call bl_error("Invalid BC's - both lo and hi need to be periodic")
     end if
     if (phys_bc(i,1) .eq. -1 .and. phys_bc(i,2) .eq. -1) then
        is_periodic(i) = .true.
     end if
  end do

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

  ! build boundary conditions
  call bc_tower_init(the_bc_tower,1,dim,phys_bc)
  call bc_tower_level_build(the_bc_tower,1,la)

  ! build multifab with 1 component and 1 ghost cell
  call multifab_build(phi,la,1,1)
  
  ! initialze phi
  call init_phi(phi,dx,prob_lo,the_bc_tower)

  istep = 0
  time = 0.d0

  ! choose a time step with a diffusive CFL of 0.9
  dt = 0.9d0*dx**2/(2.d0*dim)

  ! write out plotfile 0
  if (plot_int .gt. 0) then
     call write_plotfile(la,phi,istep,dx,time,prob_lo,prob_hi)
  end if

  do istep=1,nsteps

     ! we only want one processor to write to screen
     if ( parallel_IOProcessor() ) then
        print*,'Advancing time step',istep,'with dt=',dt
     end if
     
     ! advance phi
     call advance(phi,dx,dt,the_bc_tower)

     time = time + dt

     if (plot_int .gt. 0) then
        if (mod(istep,plot_int) .eq. 0 .or. istep .eq. nsteps) then
           ! write out plotfile
           call write_plotfile(la,phi,istep,dx,time,prob_lo,prob_hi)
        end if
     end if

  end do

  ! make sure to destroy the multifab or you'll leak memory
  call destroy(phi)
  call destroy(la)
  call bc_tower_destroy(the_bc_tower)

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
