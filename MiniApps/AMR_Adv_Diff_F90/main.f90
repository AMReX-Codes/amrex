program main
 
  use multifab_module
  use bl_IO_module
  use ml_layout_module
  use init_phi_module
  use write_plotfile_module
  use advance_module
  use define_bc_module
  use bndry_reg_module
  use make_new_grids_module
  use regrid_module

  ! This is where we keep the problem-specific values
  use prob_module

  implicit none

  ! stuff you can set with the inputs file (otherwise use default values below)
  integer    :: max_levs, dim, nsteps, plot_int, n_cell, max_grid_size
  integer    :: amr_buf_width, cluster_minwidth, cluster_blocking_factor
  real(dp_t) :: cluster_min_eff
  integer    :: regrid_int
  integer    :: bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi, bc_z_lo, bc_z_hi

  ! subcycling
  logical    :: do_subcycling
  integer    :: num_substeps

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

  ! will be allocated with max_levs components
  real(dp_t)    , allocatable :: dx(:), dt(:)
  type(multifab), allocatable :: phi_old(:)
  type(multifab), allocatable :: phi_new(:)

  ! logical multifab indicating which cells are refined
  ! will be allocated with max_levs components

  integer    :: istep,i,n,nl,nlevs
  logical    :: new_grid
  real(dp_t) :: time,start_time,run_time,run_time_IOproc
  real(dp_t) :: dt_adv, dt_diff
  
  type(box)         :: bx
  type(ml_boxarray) :: mba
  type(ml_layout)   :: mla
  type(layout), allocatable :: la_array(:)

  ! Boundary register to hold the fine flux contributions at the coarse resolution
  type(bndry_reg), allocatable :: bndry_flx(:)

  type(bc_tower) :: the_bc_tower

  namelist /probin/ max_levs, dim, nsteps, plot_int, n_cell, max_grid_size, amr_buf_width, &
       cluster_minwidth, cluster_blocking_factor, cluster_min_eff, regrid_int, &
       bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi, bc_z_lo, bc_z_hi, do_subcycling

  ! if running in parallel, this will print out the number of MPI 
  ! processes and OpenMP threads
  call boxlib_initialize()

  ! parallel_wtime() returns the number of wallclock-time seconds since
  ! the program began
  start_time = parallel_wtime()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! default values - will get overwritten by the inputs file

  max_levs      = 3
  dim           = 2
  nsteps        = 1000
  plot_int      = 100
  n_cell        = 256
  max_grid_size = 64
  amr_buf_width = 2

  cluster_minwidth = 16
  cluster_blocking_factor = 8
  cluster_min_eff = 0.7d0
  regrid_int = 4

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

  do_subcycling = .true.
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

  ! now that we have max_levs, we can allocate these
  allocate(dx(max_levs))
  allocate(dt(max_levs))
  allocate(phi_old(max_levs))
  allocate(phi_new(max_levs))
  allocate(la_array(max_levs))

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

  call cluster_set_minwidth(cluster_minwidth)
  call cluster_set_blocking_factor(cluster_blocking_factor)
  call cluster_set_min_eff(cluster_min_eff)

  ! store amr refinement ratio in ml_boxarray_module
  ! we use refinement ratio of 2 in every direction between all levels
  call amr_ref_ratio_init(max_levs,dim,2)

  ! tell mba about max_levs and dimensionality of problem
  call ml_boxarray_build_n(mba,max_levs,dim)
  ! tell mba about the ref_ratio between levels
  call ml_boxarray_set_ref_ratio(mba)

  ! physical problem is a box on (-1,-1) to (1,1)
  prob_lo(:) = -1.d0
  prob_hi(:) =  1.d0

  ! set grid spacing at each level
  ! the grid spacing is the same in each direction
  dx(1) = (prob_hi(1)-prob_lo(1)) / n_cell
  do n=2,max_levs
     dx(n) = dx(n-1) / mba%rr(n-1,1)
  end do

  ! tell the_bc_tower about max_levs, dim, and phys_bc
  call bc_tower_init(the_bc_tower,max_levs,dim,phys_bc)

  ! create a box from (0,0) to (n_cell-1,n_cell-1)
  lo(:) = 0
  hi(:) = n_cell-1
  bx = make_box(lo,hi)

  ! tell mba about the problem domain at each level
  mba%pd(1) = bx
  do n=2,max_levs
     mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
  enddo

  ! initialize the boxarray at level 1 to be one single box
  call boxarray_build_bx(mba%bas(1),bx)

  ! overwrite the boxarray at level 1 to respect max_grid_size
  call boxarray_maxsize(mba%bas(1),max_grid_size)

  ! build the level 1 layout
  call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),is_periodic)

  ! build the level 1 multifab with 1 component and 2 ghost cells
  call multifab_build(phi_new(1),la_array(1),1,2)

  ! define level 1 of the_bc_tower
  call bc_tower_level_build(the_bc_tower,1,la_array(1))

  ! initialize phi_new on level 1
  ! Here we choose to initialize phi_new because at the beginning of each time step
  !      we will copy phi_new into phi_old
  call init_phi_on_level(phi_new(1),dx(1),prob_lo,the_bc_tower%bc_tower_array(1))

  nl = 1
  new_grid = .true.

  do while ( (nl .lt. max_levs) .and. (new_grid) )

     ! determine whether we need finer grids based on tagging criteria
     ! if so, return new_grid=T and the la_array(nl+1)
     call make_new_grids(new_grid,la_array(nl),la_array(nl+1),phi_new(nl),dx(nl), &
                         amr_buf_width,mba%rr(nl,1),nl,max_grid_size)
     
     if (new_grid) then

        ! tell mba about the finer level boxarray
        call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

        ! Build the level nl+1 data
        call multifab_build(phi_new(nl+1),la_array(nl+1),1,2)
        
        ! define level nl+1 of the_bc_tower
        call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))
            
        ! initialize phi on level nl+1
        call init_phi_on_level(phi_new(nl+1),dx(nl+1),prob_lo,the_bc_tower%bc_tower_array(nl+1))

        ! increment current level counter
        nl = nl+1

     endif

  end do

  ! the current number of levels in the simulation is nlevs, not necessarily max_levs
  nlevs = nl

  ! destroy phi - we are going to build it again using the new multilevel
  ! layout after we have tested and reconfigured the grids due to proper nesting
  do n=1,nlevs
     call multifab_destroy(phi_new(n))
  end do

  if (nlevs .ge. 3) then
     ! check for proper nesting
     call enforce_proper_nesting(mba,la_array,max_grid_size)
  end if

  do n = 1,nlevs
     call destroy(la_array(n))
  end do

  call ml_layout_restricted_build(mla,mba,nlevs,is_periodic)

  call destroy(mba)

  ! this makes sure the boundary conditions are properly defined everywhere
  do n = 1,nlevs
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! Note that we need two ghost cells because we construct fluxes for the advective fluxes
  ! If we were doing pure explicit diffusion we would use just one ghost cell.
  do n=1,nlevs
     call multifab_build(phi_new(n),mla%la(n),1,2)
     call multifab_build(phi_old(n),mla%la(n),1,2)
  end do

  call init_phi(mla,phi_new,dx,prob_lo,the_bc_tower)

  istep = 0
  time = 0.d0

  if (.not.do_subcycling) then

     ! Choose a time step with a local advective CFL of 0.9 
     dt_adv  = 0.9d0*dx(max_levs) / max(abs(uadv),abs(vadv))

     ! Choose a time step with a local diffusive CFL of 0.9 
     dt_diff = 0.9d0*dx(max_levs)**2/(2.d0*dim*mu)

     ! Set the dt for all levels based on the criterion at the finest level
     ! We must use the most restrictive of the two time steps constraints
     dt(:) = min(dt_adv,dt_diff)

  else

     ! num_substeps is the number of fine steps (at level n+1) for each coarse
     ! (at level n) time step in subcycling algorithm. In other words, it is
     ! the time refinement ratio, which depends on the problem type being solved.

     ! Choose a time step with a local advective CFL of 0.9 
     dt_adv  = 0.9d0*dx(1) / max(abs(uadv),abs(vadv))

     ! Choose a time step with a local diffusive CFL of 0.9 
     dt_diff = 0.9d0*dx(1)**2/(2.d0*dim*mu)

     ! We must use the most restrictive of the two time steps constraints
     dt(1) = min(dt_adv,dt_diff)

     ! In the current implementation, num_substeps must be set to 
     !    2 (for an explicit hyperbolic problem) or 
     !    4 (for an explicit parabolic  problem) 
     num_substeps = 4

     do n = 2, max_levs 
        dt(n) = dt(n-1) / dble(num_substeps)
     end do
  endif

  ! Write plotfile at t=0
  call write_plotfile(mla,phi_new,istep,dx,time,prob_lo,prob_hi)

  ! Build the bndry_reg multifabs which will store the flux information from the
  !       fine grids at the coarse resolution.
  allocate(bndry_flx(2:nlevs))
  do n = 2, nlevs
     call flux_reg_build(bndry_flx(n),mla%la(n),mla%la(n-1),mla%mba%rr(n-1,:), &
          ml_layout_get_pd(mla,n-1),nc=1)
  end do

  do istep=1,nsteps

     ! regrid
     if ( istep > 1 .and. max_levs > 1 .and. regrid_int > 0 .and. &
          (mod(istep-1,regrid_int) .eq. 0) ) then

        ! Destroy phi_old before regridding
        do n = 1,nlevs
          call multifab_destroy(phi_old(n))
        end do

        do n=2,nlevs
           call destroy(bndry_flx(n))       
        end do

        call regrid(mla,phi_new, &
                    nlevs,max_levs,dx,the_bc_tower,amr_buf_width,max_grid_size)

        do n = 2, nlevs
           call flux_reg_build(bndry_flx(n),mla%la(n),mla%la(n-1),mla%mba%rr(n-1,:), &
                ml_layout_get_pd(mla,n-1),nc=1)
        end do

        ! Create new phi_old after regridding
        do n = 1,nlevs
          call multifab_build(phi_old(n),mla%la(n),1,2)
        end do

     end if

     ! we only want one processor to write to screen
     if ( parallel_IOProcessor() ) then
        print*,'  '
        print*,' STEP = ', istep, ' TIME = ',time, ' DT = ',dt(1)
     end if

     ! Advance phi by one coarse time step
     call advance(mla,phi_old,phi_new,bndry_flx,dx,dt,the_bc_tower,do_subcycling,num_substeps)

     time = time + dt(1)

     ! Write plotfile after every plot_int time steps, or if the simulation is done.
     if (plot_int .gt. 0) then
        if (mod(istep,plot_int) .eq. 0 .or. istep .eq. nsteps) then
           call write_plotfile(mla,phi_new,istep,dx,time,prob_lo,prob_hi)
        end if
     end if

  end do

  ! Make sure to destroy the multifabs or you'll leak memory
  do n = 1,nlevs
     call destroy(phi_old(n))
     call destroy(phi_new(n))
  end do

  call bc_tower_destroy(the_bc_tower)

  do n = 2, nlevs
     call destroy(bndry_flx(n))
  end do
  deallocate(bndry_flx)

  call destroy(mla)

  deallocate(lo,hi,is_periodic,prob_lo,prob_hi)

  ! deallocate temporary boxarrays and communication mappings
  call layout_flush_copyassoc_cache ()

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

  ! We put this here for debugging to make sure our print statements get seen.
  call flush(6)

  call boxlib_finalize()

end program main
