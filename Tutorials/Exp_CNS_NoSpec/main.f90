program main

  use boxlib
  use parallel
  use multifab_module
  use bl_IO_module
  use layout_module
  use init_data_module
  use write_plotfile_module
  use advance_module
  use bl_prof_module

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
  ! We have five components (no species).
  !
  integer, parameter :: NC  = 5

  integer            :: nsteps, plot_int, n_cell, max_grid_size
  integer            :: un, farg, narg
  logical            :: need_inputs_file, found_inputs_file
  character(len=128) :: inputs_file_name
  integer            :: i, lo(DM), hi(DM), istep
  double precision   :: prob_lo(DM), prob_hi(DM), cfl, eta, alam
  double precision   :: dx(DM), dt, time, start_time, end_time
  logical            :: is_periodic(DM)
  type(box)          :: bx
  type(boxarray)     :: ba
  type(layout)       :: la
  type(multifab)     :: U

  type(bl_prof_timer), save :: bpt, bpt_init_data, bpt_sleep

  !
  ! What's settable via an inputs file.
  !
  namelist /probin/ nsteps, plot_int, n_cell, max_grid_size, cfl, eta, alam

  call boxlib_initialize()
  call bl_prof_initialize(on = .true.)

  call build(bpt, "bpt_main")


  start_time = parallel_wtime()

  call build(bpt_sleep, "bpt_sleep")
  call Sleep(5)
  call destroy(bpt_sleep)

  !
  ! Namelist default values -- overwritable via inputs file.
  !
  nsteps        = 100
  plot_int      = 10
  n_cell        = 32
  max_grid_size = 32
  cfl           = 0.5d0
  eta           = 1.8d-4 ! Diffusion coefficient.
  alam          = 1.5d2  ! Diffusion coefficient.
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

  if ( parallel_IOProcessor() ) then
     write(6,probin)
  end if
  !
  ! Physical problem is a box on (-1,-1) to (1,1), periodic on all sides.
  !
  prob_lo     = -0.1d0
  prob_hi     =  0.1d0
  is_periodic = .true.
  !
  ! Create a box from (0,0) to (n_cell-1,n_cell-1).
  !
  lo = 0
  hi = n_cell-1
  bx = make_box(lo,hi)

  do i = 1,DM
     dx(i) = (prob_hi(i)-prob_lo(i)) / n_cell
  end do

  call boxarray_build_bx(ba,bx)

  call boxarray_maxsize(ba,max_grid_size)

  call layout_build_ba(la,ba,boxarray_bbox(ba),pmask=is_periodic)

  call destroy(ba)

  call multifab_build(U,la,NC,NG)
  
  call build(bpt_init_data, "bpt_init_data")
  call init_data(U,dx,prob_lo,prob_hi)
  call destroy(bpt_init_data)

  istep = 0
  time  = 0.d0

  if (plot_int > 0) then
     call write_plotfile(U,istep,dx,time,prob_lo,prob_hi)
  end if

  do istep=1,nsteps

     if (parallel_IOProcessor()) then
        print*,'Advancing time step',istep,'time = ',time
     end if
     
     call advance(U,dt,dx,cfl,eta,alam)

     time = time + dt

     if (plot_int > 0) then
        if (mod(istep,plot_int) .eq. 0 .or. istep .eq. nsteps) then
           call write_plotfile(U,istep,dx,time,prob_lo,prob_hi)
        end if
     end if

  end do

  call destroy(U)
  call destroy(la)

  end_time = parallel_wtime()

  if ( parallel_IOProcessor() ) then
     print*,"Run time (s) =",end_time-start_time
  end if


  call destroy(bpt)
  call bl_prof_glean("bl_prof_report")
  call bl_prof_finalize()


  call boxlib_finalize()

end program main
