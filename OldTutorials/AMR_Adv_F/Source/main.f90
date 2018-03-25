program main
 
  use bl_IO_module
  use ml_layout_module
  use init_phi_module
  use write_plotfile_module
  use advance_module
  use define_bc_module
  use set_velocity_module
  use bndry_reg_module
  use make_new_grids_module
  use regrid_module

  implicit none

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
  type(multifab), allocatable :: velocity(:,:)

  ! logical multifab indicating which cells are refined
  ! will be allocated with max_levs components

  integer    :: istep,i,n,nl,nlevs
  logical    :: new_grid
  real(dp_t) :: time,start_time,run_time,run_time_IOproc
  real(dp_t) :: vmax
  
  type(box)         :: bx
  type(ml_boxarray) :: mba
  type(ml_layout)   :: mla
  type(layout), allocatable :: la_array(:)

  ! Boundary register to hold the fine flux contributions at the coarse resolution
  type(bndry_reg), allocatable :: bndry_flx(:)

  ! holds description of boundary conditions
  type(bc_tower) :: the_bc_tower

  ! for problems with subcycling,
  ! num_substeps is the number of fine steps (at level n+1) for each coarse
  ! (at level n) time step in subcycling algorithm. In other words, it is
  ! the time refinement ratio, which depends on the problem type being solved.
  ! In the current implementation, num_substeps must be set to 
  !    2 (for an explicit hyperbolic problem)
  integer, parameter :: num_substeps = 2

  ! namelist parameters you can set with the inputs file
  ! (otherwise use default values below)
  ! see description in comments below
  integer    :: dim, max_levs, plot_int, n_cell, max_grid_size, nsteps, regrid_int
  real(dp_t) :: cfl,stop_time
  logical    :: do_subcycling
  integer    :: bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi, bc_z_lo, bc_z_hi
  real(dp_t) :: prob_lo_x, prob_lo_y, prob_lo_z, prob_hi_x, prob_hi_y, prob_hi_z
  integer    :: amr_buf_width, cluster_minwidth, cluster_blocking_factor
  real(dp_t) :: cluster_min_eff

  ! dummy indices using for reading in inputs file
  integer :: un, farg, narg
  logical :: need_inputs_file, found_inputs_file
  character(len=128) :: inputs_file_name

  namelist /probin/ dim                      ! dimensionality of problem
  namelist /probin/ max_levs                 ! total number of AMR levels
  namelist /probin/ plot_int                 ! plotfile interval
  namelist /probin/ n_cell                   ! number of cells on each side
  namelist /probin/ max_grid_size            ! max number of cells on each side
  namelist /probin/ nsteps                   ! number of time steps
  namelist /probin/ cfl                      ! advective cfl
  namelist /probin/ stop_time                ! simulation stop time
  namelist /probin/ regrid_int               ! how often to regrid
  namelist /probin/ do_subcycling            ! use subcycling
  ! allowable boundary condition options currently are
  ! -1 = PERIODIC
  namelist /probin/ bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi, bc_z_lo, bc_z_hi
  namelist /probin/ prob_lo_x, prob_lo_y, prob_lo_z, prob_hi_x, prob_hi_y, prob_hi_z
  namelist /probin/ amr_buf_width            ! number of buffer cells between successive levels
  namelist /probin/ cluster_minwidth         ! minimimum size of an AMR grid
  namelist /probin/ cluster_blocking_factor  ! grids must be an integer multiple
  namelist /probin/ cluster_min_eff          ! larger value means more, smaller, tighter grids

  ! if running in parallel, this will print out the number of MPI 
  ! processes and OpenMP threads
  call boxlib_initialize()

  ! parallel_wtime() returns the number of wallclock-time seconds since
  ! the program began
  start_time = parallel_wtime()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! default values - will get overwritten by the inputs file

  dim           = 2
  max_levs      = 3
  plot_int      = 1
  n_cell        = 64
  max_grid_size = 16
  nsteps        = 10
  cfl           = 0.9d0
  stop_time     = 2.d0

  regrid_int = 4
  do_subcycling = .true.

  bc_x_lo = -1 ! PERIODIC
  bc_x_hi = -1 ! PERIODIC
  bc_y_lo = -1 ! PERIODIC
  bc_y_hi = -1 ! PERIODIC
  bc_z_lo = -1 ! PERIODIC
  bc_z_hi = -1 ! PERIODIC

  prob_lo_x = 0.d0
  prob_lo_y = 0.d0
  prob_lo_z = 0.d0
  prob_hi_x = 1.d0
  prob_hi_y = 1.d0
  prob_hi_z = 1.d0

  amr_buf_width = 2
  cluster_minwidth = 16
  cluster_blocking_factor = 8
  cluster_min_eff = 0.7d0
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
  allocate(velocity(max_levs,dim))
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

  ! physical problem is a box on (0,0) to (1,1)
  prob_lo(1) = prob_lo_x
  prob_hi(1) = prob_hi_x
  prob_lo(2) = prob_lo_y
  prob_hi(2) = prob_hi_y
  if (dim .gt. 2) then
      prob_lo(3) = prob_lo_z
      prob_hi(3) = prob_hi_z
  end if

  ! set grid spacing at each level
  ! the grid spacing is the same in each direction
  dx(1) = (prob_hi(1)-prob_lo(1)) / n_cell
  do n=2,max_levs
     dx(n) = dx(n-1) / mba%rr(n-1,1)
  end do

  ! We initialize dt to very large here because we will not let dt grow 
  !    by more than 10% when we compute it in each time step
  dt(:) = HUGE(1.0) 

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

  ! build the level 1 multifab with 1 component and 3 ghost cells
  call multifab_build(phi_new(1),la_array(1),1,3)

  ! define level 1 of the_bc_tower
  call bc_tower_level_build(the_bc_tower,1,la_array(1))

  ! initialize phi_new on level 1
  ! Here we choose to initialize phi_new because at the beginning of each time step
  !      we will copy phi_new into phi_old
  call init_phi_on_level(phi_new(1),dx(1),prob_lo)

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
        call multifab_build(phi_new(nl+1),la_array(nl+1),1,3)
        
        ! define level nl+1 of the_bc_tower
        call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))

        ! initialize phi on level nl+1
        call init_phi_on_level(phi_new(nl+1),dx(nl+1),prob_lo)

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

  ! build ml_layout with nlevs levels, NOT max_levs levels
  call ml_layout_restricted_build(mla,mba,nlevs,is_periodic)

  call destroy(mba)

  ! this makes sure the boundary conditions are properly defined everywhere
  do n = 1,nlevs
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! Note that we need two ghost cells because we construct fluxes for the advective fluxes
  ! If we were doing pure explicit diffusion we would use just one ghost cell.
  do n=1,nlevs
     call multifab_build(phi_new(n),mla%la(n),1,3)
     call multifab_build(phi_old(n),mla%la(n),1,3)
     do i=1,dim
        call multifab_build_edge(velocity(n,i),mla%la(n),1,1,i)
     end do
  end do

  call init_phi(mla,phi_new,dx,prob_lo,the_bc_tower)

  istep = 0
  time = 0.d0

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

     ! Compute dt using the CFL condition and making sure we don't go past stop_time
     call compute_dt(velocity,time,dt)

     ! regrid
     if ( istep > 1 .and. max_levs > 1 .and. regrid_int > 0 .and. &
          (mod(istep-1,regrid_int) .eq. 0) ) then

        ! Destroy phi_old and velocity before regridding
        do n = 1,nlevs
          call multifab_destroy(phi_old(n))
          do i = 1, dim
             call multifab_destroy(velocity(n,i))
          end do
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

        ! Create new phi_old and velocity after regridding
        do n = 1,nlevs
          call multifab_build(phi_old(n),mla%la(n),1,3)
          do i=1,dim
             call multifab_build_edge(velocity(n,i),mla%la(n),1,1,i)
          end do
        end do

     end if

     ! we only want one processor to write to screen
     if ( parallel_IOProcessor() ) then
        print*,'  '
        print*,' STEP = ', istep, ' TIME = ',time, ' DT = ',dt(1)
     end if

     ! Advance phi by one coarse time step
     call advance(mla,phi_old,phi_new,velocity,bndry_flx,dx,dt,time,the_bc_tower, &
                  do_subcycling,num_substeps)

     time = time + dt(1)

     ! Write plotfile after every plot_int time steps, or if the simulation is done.
     if (plot_int .gt. 0) then
        if (mod(istep,plot_int) .eq. 0 .or. istep .eq. nsteps .or. time .ge. stop_time) then
           call write_plotfile(mla,phi_new,istep,dx,time,prob_lo,prob_hi)
        end if
     end if

     if (time .ge. stop_time) exit

  end do

  ! Make sure to destroy the multifabs or you'll leak memory
  do n = 1,nlevs
     call destroy(phi_new(n))
     call destroy(phi_old(n))
     do i = 1, dim
        call multifab_destroy(velocity(n,i))
     end do
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
  flush(6)

  call boxlib_finalize()

contains

  subroutine compute_dt(vel,time_n,dt)

      type(multifab) , intent(inout) :: vel(:,:)
      real(dp_t)     , intent(in   ) :: time_n
      real(dp_t)     , intent(inout) :: dt(:)

      ! Local variables
      real(dp_t)                     :: dt_old(mla%nlevel)
      real(dp_t)                     :: cfl_grow_fac

      ! We don't let dt grow by more than 10% over the previous dt
      cfl_grow_fac = 1.1d0

      ! Here we set the velocity so we can find the max vel to use in setting dt
      call set_velocity(mla,vel,dx,time_n)

      vmax = -HUGE(1.d0)
      do n = 1, nlevs
         do i=1,dim
            vmax = max(vmax,norm_inf(vel(n,i)))
         end do
      end do

      ! Store the old dt so we make sure we don't grow too fast
      dt_old(1:nlevs) = dt(1:nlevs)

      if (.not.do_subcycling) then

         ! Choose a time step with a local advective CFL of "cfl" (set in the inputs file)
         ! Set the dt for all levels based on the criterion at the finest level
         dt(:) = cfl*dx(nlevs) / vmax

         ! We don't let dt grow by more than the factor of cfl_grow_fac over the previous dt
         dt(:) = min(dt(:), cfl_grow_fac * dt_old(:))

         ! Make sure we don't overshoot stop_time
         if (time+dt(1) .gt. stop_time) &
            dt(:) = stop_time - time

      else

         ! Choose a time step with a local advective CFL of "cfl" (set in the inputs file)
         dt(1) = cfl*dx(1) / vmax

         ! We don't let dt grow by more than the factor of cfl_grow_fac over the previous dt
         dt(1) = min(dt(1), cfl_grow_fac * dt_old(1))

         ! Make sure we don't overshoot stop_time
         if (time+dt(1) .gt. stop_time) &
             dt(1) = stop_time - time
    
         ! Set the dt for all levels based on refinement ratio
         do n = 2, nlevs 
            dt(n) = dt(n-1) / dble(num_substeps)
         end do

      endif

  end subroutine compute_dt

end program main
