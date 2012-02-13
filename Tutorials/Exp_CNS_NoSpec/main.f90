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
  !
  ! We only support 3-D.
  !
  integer, parameter :: DM = 3
  !
  ! We need four grow cells.
  !
  integer, parameter :: NG  = 4
  !
  ! Stuff you can set with the inputs file (otherwise use default values below).
  !
  integer :: nsteps, plot_int, n_cell, max_grid_size

  integer :: un, farg, narg
  logical :: need_inputs_file, found_inputs_file
  character(len=128) :: inputs_file_name

  integer :: i, lo(DM), hi(DM), istep

  double precision :: prob_lo(DM), prob_hi(DM)
  double precision :: dx(DM), dt, time, start_time, end_time
  
  logical :: is_periodic(DM)

  type(box)      :: bx
  type(boxarray) :: ba
  type(layout)   :: la
  type(multifab) :: data

  namelist /probin/ nsteps, plot_int, n_cell, max_grid_size

  call boxlib_initialize()

  start_time = parallel_wtime()
  !
  ! Default values - will get overwritten by the inputs file.
  !
  nsteps        = 10
  plot_int      = 1
  n_cell        = 64
  max_grid_size = 32
  !
  ! Read inputs file and overwrite any default values.
  !
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
  !
  ! Physical problem is a box on (-1,-1) to (1,1), periodic on all sides.
  !
  prob_lo     = -1.d0
  prob_hi     =  1.d0
  is_periodic = .true.
  !
  ! Create a box from (0,0) to (n_cell-1,n_cell-1).
  !
  lo = 0
  hi = n_cell-1
  bx = make_box(lo,hi)
  !
  ! The grid spacing is the same in each direction.
  !
  do i = 1,DM
     dx(i) = (prob_hi(i)-prob_lo(i)) / n_cell
  end do
  !
  ! Initialize the boxarray to be one single box.
  !
  call boxarray_build_bx(ba,bx)
  !
  ! Overwrite the boxarray to respect max_grid_size.
  !
  call boxarray_maxsize(ba,max_grid_size)

  call layout_build_ba(la,ba,pmask=is_periodic)

  call destroy(ba)

  ! build multifab with 2 components and 6 ghost cells
  call multifab_build(data,la,2,6)
  
  call init_data(data,dx,prob_lo)

  istep = 0
  time  = 0.d0
  dt    = 0.1d0*dx(1)

  call write_plotfile(la,data,istep,dx,time,prob_lo,prob_hi)

  do istep=1,nsteps

     if (parallel_IOProcessor()) then
        print*,'Advancing time step',istep
     end if
     
     call advance(data,dt,dx)

     time = time + dt

     if (mod(istep,plot_int) .eq. 0 .or. istep .eq. nsteps) then
        call write_plotfile(la,data,istep,dx,time,prob_lo,prob_hi)
     end if

  end do

  call destroy(data)
  call destroy(la)

  end_time = parallel_wtime()

  call boxlib_finalize()

  if ( parallel_IOProcessor() ) then
     print*,"Run time (s) =",end_time-start_time
  end if

end program main
