program main

  use boxlib
  use parallel
  use multifab_module
  use bl_IO_module
  use layout_module
  use init_phi_module
  use write_plotfile_module
  use advance_module
  use define_bc_module

  implicit none

  ! stuff you can set with the inputs file (otherwise use default values below)
  integer :: nlevs, dim, nsteps, plot_int, n_cell, max_grid_size
  integer :: bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi, bc_z_lo, bc_z_hi

  ! dummy indices using for reading in inputs file
  integer :: un, farg, narg
  logical :: need_inputs_file, found_inputs_file
  character(len=128) :: inputs_file_name

  ! will be allocated with dim components
  integer       , allocatable :: lo(:), hi(:)
  logical       , allocatable :: is_periodic(:)
  real(dp_t)    , allocatable :: prob_lo(:), prob_hi(:)

  ! will be allocated with (dim,2) components
  integer       , allocatable :: phys_bc(:,:)

  ! will be allocated with nlevs components
  real(dp_t)    , allocatable :: dx(:)
  type(multifab), allocatable :: phi(:)

  integer    :: istep,i,n,n_cell_level
  real(dp_t) :: dt,time,start_time,end_time
  
  type(box)         :: bx
  type(ml_boxarray) :: mba
  type(ml_layout)   :: mla

  type(bc_tower) :: the_bc_tower

  namelist /probin/ nlevs, dim, nsteps, plot_int, n_cell, max_grid_size, &
       bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi, bc_z_lo, bc_z_hi

  ! if running in parallel, this will print out the number of MPI 
  ! processes and OpenMP threads
  call boxlib_initialize()

  ! parallel_wtime() returns the number of wallclock-time seconds since
  ! the program began
  start_time = parallel_wtime()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! default values - will get overwritten by the inputs file

  nlevs         = 3
  dim           = 2
  nsteps        = 1000
  plot_int      = 100
  n_cell        = 256
  max_grid_size = 64

  ! allowable options for this example are
  ! -1 = PERIODIC
  ! 12 = OUTLET (dphi/dn=0 at boundary)
  ! 15 = NO_SLIP_WALL (wall with fixed phi=1)
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
  allocate(is_periodic(dim))
  allocate(prob_lo(dim),prob_hi(dim))
  allocate(phys_bc(dim,2))

  ! now that we have nlevs, we can allocate these
  allocate(dx(nlevs))
  allocate(phi(nlevs))

  ! physical problem is a box on (-1,-1) to (1,1)
  prob_lo(:) = -1.d0
  prob_hi(:) =  1.d0

  ! put all the domain boundary conditions into phys_bc
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

  ! tell mba how many levels and dimensionality of problem
  call ml_boxarray_build_n(mba,nlevs,dim)

  ! create a box from (0,0) to (n_cell-1,n_cell-1)
  lo(:) = 0
  hi(:) = n_cell-1
  bx = make_box(lo,hi)

  ! tell mba about the ref_ratio between levels
  ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
  ! we use refinement ratio of 2 in every direction between all levels
  do n=2,nlevs
     mba%rr(n-1,:) = 2
  enddo

  ! tell mba about the problem domain at each level
  mba%pd(1) = bx
  do n=2,nlevs
     mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
  enddo

  ! initialize the boxarray at level 1 to be one single box
  call boxarray_build_bx(mba%bas(1),bx)

  ! overwrite the boxarray at level 1 to respect max_grid_size
  call boxarray_maxsize(mba%bas(1),max_grid_size)

  ! now build the boxarray at other levels
  n_cell_level = n_cell
  do n=2,nlevs

     ! length of the problem domain at this level
     n_cell_level = n_cell_level * mba%rr(n-1,1)

     ! logic to refine the central n_cell x n_cell regio
     lo(:) = n_cell_level/2-n_cell/2
     hi(:) = n_cell_level/2+n_cell/2-1
     bx = make_box(lo,hi)

     ! initialize the boxarray at level n to be one single box
     call boxarray_build_bx(mba%bas(n),bx)

     ! overwrite the boxarray at level n to respect max_grid_size
     call boxarray_maxsize(mba%bas(n),max_grid_size)

  end do

  ! build the ml_layout, mla
  call ml_layout_build(mla,mba,is_periodic)

  ! don't need this anymore - free up memory
  call destroy(mba)

  ! the grid spacing is the same in each direction
  dx(1) = (prob_hi(1)-prob_lo(1)) / n_cell
  do n=2,nlevs
     dx(n) = dx(n-1) / mla%mba%rr(n-1,1)
  end do

  ! build boundary conditions
  call bc_tower_init(the_bc_tower,nlevs,dim,phys_bc)
  do n=1,nlevs
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! build multifab with 1 component and 1 ghost cell
  do n=1,nlevs
     call multifab_build(phi(n),mla%la(n),1,1)
  end do
  
  ! initialze phi
  call init_phi(mla,phi,dx,prob_lo,the_bc_tower)

  istep = 0
  time = 0.d0

  ! choose a time step with a diffusive CFL of 0.9
  dt = 0.9d0*dx(nlevs)**2/(2.d0*dim)

  ! write out plotfile 0
  call write_plotfile(mla,phi,istep,dx,time,prob_lo,prob_hi)

  do istep=1,nsteps

     ! we only want one processor to write to screen
     if ( parallel_IOProcessor() ) then
        print*,'Advancing time step',istep,'with dt=',dt
     end if
     
     ! advance phi
     call advance(mla,phi,dx,dt,the_bc_tower)

     time = time + dt

     if (mod(istep,plot_int) .eq. 0 .or. istep .eq. nsteps) then
        ! write out plotfile
        call write_plotfile(mla,phi,istep,dx,time,prob_lo,prob_hi)
     end if

  end do

  ! make sure to destroy the multifab or you'll leak memory
  do n=1,nlevs
     call destroy(phi(n))
  end do
  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

  deallocate(lo,hi,is_periodic,prob_lo,prob_hi)

  ! parallel_wtime() returns the number of wallclock-time seconds since
  ! the program began
  end_time = parallel_wtime()

  call boxlib_finalize()

  if ( parallel_IOProcessor() ) then
     print*,"Run time (s) =",end_time-start_time
  end if

end program main
