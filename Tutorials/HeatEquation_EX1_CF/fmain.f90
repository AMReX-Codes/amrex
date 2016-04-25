
subroutine fmain () bind(c)

  use boxlib_module

  use init_phi_module, only : init_phi
  use advance_module, only : advance

  implicit none

  integer :: n_cell, max_grid_size, nsteps, plot_int
  integer, parameter :: ncomp = 1, nghost = 1  ! one component, one ghost
  integer :: istep
  double precision :: dt, time
  type(ParmParse) :: pp
  type(Box) :: domain
  type(BoxArray) :: bs
  type(Geometry) :: geom
  type(MultiFab) :: new_phi, old_phi

  ! ParmParse is way of reading inputs from the inputs file
  ! "get" means it must be set in the inputs file, whereas
  ! "query" means it may not may not be in the inputs file
  call parmparse_build(pp)

  call pp%get("n_cell", n_cell);  ! # of cells in each dimension
  call pp%get("nsteps", nsteps) ! # of steps

  max_grid_size = 32   ! default max grid size
  call pp%query("max_grid_size", max_grid_size);

  plot_int = -1 ! default to no plotfiles
  call pp%query("plot_int", plot_int);
  
  ! Define a single box covering the domain
  domain = Box((/0,0,0/), (/n_cell-1, n_cell-1, n_cell-1/))

  ! Initialize the boxarray "bs" from the single box "bx"
  call boxarray_build(bs, domain)

  ! Break up boxarray "bs" into chunks no larger than "max_grid_size" along a direction
  call bs%maxSize(max_grid_size)

  ! This defines a Geometry object.
  call geometry_build(geom, domain)

  ! Build data multifabs
  call multifab_build(new_phi, bs, ncomp, nghost)
  call multifab_build(old_phi, bs, ncomp, nghost)

  ! Intialize data
  call init_phi(new_phi, geom)

  istep = 0
  time = 0.d0

  ! choose a time step with a diffusive CFL of 0.9
  dt = 0.9d0*geom%dx(1)**2/(2.d0*bl_num_dims)

  do istep = 1, nsteps

     if ( parallel_IOProcessor() ) then
        print*,'Advancing time step',istep,'with dt=',dt
     end if

     ! Swap the guts of multifabs so we don't have to allocate and de-allocate data
     call multifab_swap(new_phi, old_phi)

     ! advance phi
     call advance(old_phi, new_phi, geom, dt)

     time = time + dt

  end do

end subroutine fmain
