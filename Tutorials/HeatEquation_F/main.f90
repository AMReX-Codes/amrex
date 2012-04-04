program main

  use boxlib
  use parallel
  use multifab_module
  use bl_IO_module
  use layout_module
  use initialize_module
  use multifab_fill_ghost_module
  use define_bc_module
  use init_data_module
  use write_plotfile_module
  use advance_module
  use make_flux_module

  use ml_restriction_module , only : ml_cc_restriction

  implicit none

  ! stuff you can set with the inputs file (otherwise use default values below)
  integer :: nlevs, dim, nsteps, plot_int, n_cell, max_grid_size, verbose

  ! dummy indices using for reading in inputs file
  integer :: un, farg, narg
  logical :: need_inputs_file, found_inputs_file
  character(len=128) :: inputs_file_name

  integer, allocatable :: lo(:), hi(:), ref_ratio(:,:)
  integer :: i,n,grow_fac(10)
  integer :: istep
  logical :: test

  real(dp_t), allocatable :: prob_lo(:), prob_hi(:)
  real(dp_t) :: dt, time, start_time, end_time
  real(dp_t) :: coef

  type(bc_tower) :: the_bc_tower 
  
  real(dp_t), allocatable :: dx(:)

  logical   , allocatable :: pmask(:)
  logical   , allocatable :: is_neumann(:)

  type(box)     , allocatable :: bx(:)
  type(boxarray), allocatable :: ba(:)
  type(box)     , allocatable :: pd(:)
  type(layout)  , allocatable :: la(:)
  type(multifab), allocatable :: data(:)
  type(multifab), allocatable :: flux(:,:)

  namelist /probin/ nlevs, dim, nsteps, plot_int, n_cell, max_grid_size, verbose

  ! if running in parallel, this will print out the number of MPI 
  ! processes and OpenMP threads
  call boxlib_initialize()

  ! parallel_wtime() returns the number of wallclock-time seconds since
  ! the program began
  start_time = parallel_wtime()

  ! default values - will get overwritten by the inputs file
  nlevs         = 1
  dim           = 2
  nsteps        = 10
  plot_int      = 1
  n_cell        = 64
  max_grid_size = 32
  verbose       = 0
  nlevs         = 1

  coef = 1.d0

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

  allocate(bx(nlevs))
  allocate(ba(nlevs))
  allocate(pd(nlevs))
  allocate(la(nlevs))
  allocate(data(nlevs))
  allocate(flux(nlevs,dim))
  allocate(dx(nlevs))

  ! now that we have dim, we can allocate these
  allocate(lo(dim),hi(dim))
  allocate(pmask(dim),is_neumann(dim))
  allocate(prob_lo(dim),prob_hi(dim))
  allocate(ref_ratio(nlevs-1,dim))

  ! physical problem is a box on (-1,-1) to (1,1), periodic on all sides
  prob_lo(:) = 0.d0
  prob_hi(:) =  1.d0
  pmask(1) = .true.
  pmask(2) = .true.

  ref_ratio(:,:) = 2

  ! create a box from (0,0) to (n_cell-1,n_cell-1)
  lo(:) = 0
  hi(:) = n_cell-1
  bx(1) = make_box(lo,hi)
  pd(1) = bx(1)

  ! the grid spacing is the same in each direction
  dx(1) = (prob_hi(1)-prob_lo(1)) / n_cell
  write(*,*) 'Prob_lo / Prob_hi ',prob_lo(1), prob_hi(1)
  write(*,*) '           nlevs  ',nlevs
  write(*,*) ' Level 0   n_cell ',n_cell

  ! These boxes cover the entire domain
  do n = 2,nlevs
     bx(n) = refine(bx(n-1),ref_ratio(n-1,:))
     pd(n) = bx(n)
  end do

  if (nlevs > 1) &
     grow_fac(2) = n_cell * ref_ratio(1,1) / 4
     ! grow_fac(2) = 3 * n_cell * ref_ratio(1,1) / 8
  if (nlevs > 2) &
     grow_fac(3) = n_cell * ref_ratio(1,1) * ref_ratio(2,1) * 3 / 8

  ! These boxes are nested and centered at the center of the domain
  !  (If you comment out the next lines then the boxes will cover the entire domain)
  do n = 2,nlevs
     bx(n) = grow(bx(n),-grow_fac(n))
  end do
  !  (End of part to comment out)

  ! dx gets smaller by ref_ratio at each higher level
  do n = 2,nlevs
     dx(n) = dx(n-1) / dble(ref_ratio(n-1,1))
  end do

  ! Initialize the boxarray at each level to be one single box, 
  !   then divide it into pieces no larger than max_grid_size
  do n = 1,nlevs
     print *,'At level ',n, 'we have box ' 
     call print(bx(n))

     call boxarray_build_bx(ba(n),bx(n))

     ! overwrite the boxarray to respect max_grid_size
     call boxarray_maxsize(ba(n),max_grid_size)
  end do

  ! build the layout, la
  do n = 1,nlevs
     call print(pd(n),'PD AT N')
     call layout_build_ba(la(n),ba(n),pd(n),pmask=pmask)
     call destroy(ba(n))
  end do

  call initialize_bc(dim,the_bc_tower,nlevs,pmask)
  do n = 1,nlevs
     call bc_tower_level_build(the_bc_tower,n,la(n))
  end do

  ! build multifab with 1 component and 6 ghost cells
  do n = 1,nlevs
     call multifab_build(data(n),la(n),1,1)
  end do

  ! initialize data
  do n = 1,nlevs
     call init_data(n,data(n),dx(n),prob_lo)
  end do

  istep  = 0
  time   = 0.d0

  dt = 0.1d0 * dx(nlevs)**2

  ! write out plotfile 0
  call write_plotfile(nlevs,ref_ratio(:,1),la,data,istep,dx,time,prob_lo,prob_hi)

  ! Build the flux arrays
  do n = 1,nlevs
     do i = 1,dim
        call multifab_build_edge( flux(n,i),la(n),nc=1,ng=0,dir=i)
        call setval( flux(n,i),1.d20, all=.true.)
     end do
  end do

  do istep = 1,nsteps

     if ( verbose .eq. 1 .and. (parallel_IOProcessor()).and.(mod(istep,1000) .eq. 0) ) then
        print*,'Advancing time step ',istep,' :Time = ', istep*dt
     end if

     ! Fill the ghost cells of all levels > 1
     do n = 2,nlevs
        call multifab_fill_ghost_cells(data(n),data(n-1), &
                                       nghost(data(n)),ref_ratio(n-1,:), &
                                       the_bc_tower%bc_tower_array(n-1), &
                                       the_bc_tower%bc_tower_array(n  ), &
                                       1,1,1)
     end do

     ! Construct the fluxes
     call make_fluxes(nlevs,data,flux,dx,ref_ratio)

     ! Advance the data given the fluxes
     call advance(nlevs,data,flux,dx,dt,coef)

     ! Restrict the cell-centered fine grid solution onto the coarse grid
     do n = nlevs,2,-1
         call ml_cc_restriction(data(n-1),data(n),ref_ratio(n-1,:))
     end do

     time = time + dt

     if (mod(istep,plot_int) .eq. 0 .or. istep .eq. nsteps) then
        ! write out plotfile
        test = contains_nan(data(n)%fbs(1),1,1)
        print *,'Testing data n = ',n,' test = ',test
        call write_plotfile(nlevs,ref_ratio(:,1),la,data,istep,dx,time,prob_lo,prob_hi)
     end if

  end do

  ! make sure to destroy the multifab or you'll leak memory
  do n = 1,nlevs
     call destroy(data(n))
     call destroy(la(n))
  end do

  ! Destroy the flux arrays
  do n = 1,nlevs
     do i = 1,dim
        call destroy(flux(n,i))
     end do
  end do

  deallocate(flux)
  deallocate(data,dx,la,ba,bx)
  deallocate(lo,hi,pmask,prob_lo,prob_hi)

  ! parallel_wtime() returns the number of wallclock-time seconds since
  ! the program began
  end_time = parallel_wtime()

  call boxlib_finalize()

  if ( parallel_IOProcessor() ) then
     print*,"Run time (s) =",end_time-start_time
  end if

end program main
