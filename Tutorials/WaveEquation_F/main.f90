program main

  use boxlib
  use parallel
  use multifab_module
  use bl_IO_module
  use ml_layout_module
  use init_data_module
  use write_plotfile_module
  use advance_module

  implicit none

  integer :: dim, nsteps, plot_int, n_cell, max_grid_size, verbose
  integer :: un, farg, narg, n_amr_levels, n, istep
  integer, allocatable :: lo(:), hi(:)

  real(dp_t), allocatable :: dx(:,:), prob_lo(:), prob_hi(:)
  real(dp_t) :: time, dt, start_time, end_time
  
  character(len=128) :: fname

  logical :: lexist, need_inputs
  logical, allocatable :: is_periodic(:)

  type(ml_layout) :: mla
  type(box) :: bx
  type(ml_boxarray) :: mba

  type(multifab), allocatable :: data(:)

  namelist /probin/ dim, nsteps, plot_int, n_cell, max_grid_size, verbose

  ! if running in parallel, this will print out the number of MPI 
  ! processes and OpenMP threads
  call boxlib_initialize()

  start_time = parallel_wtime()

  ! default values - will get overwritten by the inputs file
  dim           = 2
  nsteps        = 10
  plot_int      = 1
  n_cell        = 64
  max_grid_size = 32
  verbose       = 0

  ! not using AMR in this example
  n_amr_levels = 1


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
  allocate(is_periodic(dim))
  allocate(dx(n_amr_levels,dim))
  allocate(prob_lo(dim),prob_hi(dim))

  ! physical problem is a box on (-1,-1) to (1,1), periodic on all sides
  prob_lo(:) = -1.d0
  prob_hi(:) =  1.d0
  is_periodic(:) = .true.

  lo(:) = 0
  hi(:) = n_cell-1

  bx = make_box(lo,hi)

  dx(1,:) = (prob_hi(:)-prob_lo(:)) / n_cell

  ! build the ml_boxarray, mba
  ! initialize the number of levels and dimensionality
  call ml_boxarray_build_n(mba,n_amr_levels,dim)
  ! initialize the boxarray to be one single box
  call boxarray_build_bx(mba%bas(1),bx)
  ! overwrite the boxarray to respect max_grid_size
  call boxarray_maxsize(mba%bas(1),max_grid_size)
  ! initialze the problem domain
  mba%pd(1) = bx

  ! build the ml_layout, mla
  call ml_layout_restricted_build(mla,mba,n_amr_levels,is_periodic)

  allocate(data(n_amr_levels))
  do n=1,n_amr_levels
     ! build multifab with 2 components and 6 ghost cells
     call multifab_build(data(n),mla%la(n),2,6)
  end do

  ! initialze data
  call init_data(mla,dx,prob_lo,data)

  istep = 0
  time = 0.d0

  dt = 0.1d0*dx(1,1)

  ! write out plotfile 0
  call write_plotfile(mla,data,istep,dx,time,prob_lo,prob_hi)


  do istep=1,nsteps

     if (verbose .ge. 1) then
        ! we only want one processor to write to screen
        if ( parallel_IOProcessor() ) then
           print*,'Advancing time step',istep
        end if
     end if
     
     ! advance the data
     call advance(mla,dx,dt,data)

     time = time + dt

     if (mod(istep,plot_int) .eq. 0 .or. istep .eq. nsteps) then
        ! write out plotfile
        call write_plotfile(mla,data,istep,dx,time,prob_lo,prob_hi)
     end if

  end do

  do n=1,n_amr_levels
     ! make sure to destroy the multifab to avoid memory leaks
     call destroy(data(n))
  end do

  deallocate(lo,hi,dx,is_periodic,data,prob_lo,prob_hi)

  end_time = parallel_wtime()

  call boxlib_finalize()

  if ( parallel_IOProcessor() ) then
     print*,"Run time (s) =",end_time-start_time
  end if

end program main
