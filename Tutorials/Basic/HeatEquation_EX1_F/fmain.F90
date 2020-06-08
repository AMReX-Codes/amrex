
program main

  use amrex_base_module

  use init_phi_module, only : init_phi
  use advance_module, only : advance

  implicit none

  integer :: n_cell, max_grid_size, nsteps, plot_int
  integer, parameter :: ncomp = 1, nghost = 1  ! one component, one ghost
  integer :: istep
  real(amrex_real) :: dt, coeff, time
  type(amrex_parmparse) :: pp
  type(amrex_box) :: domain
  type(amrex_boxarray)  :: ba
  type(amrex_distromap) :: dm
  type(amrex_geometry)  :: geom
  type(amrex_multifab)  :: new_phi, old_phi

  call amrex_init()

  ! amrex_parmparse is way of reading inputs from the inputs file
  ! "get" means it must be set in the inputs file, whereas
  ! "query" means it may not may not be in the inputs file
  call amrex_parmparse_build(pp)

  call pp%get("n_cell", n_cell);  ! # of cells in each dimension
  call pp%get("nsteps", nsteps) ! # of steps

  max_grid_size = 32   ! default max grid size
  call pp%query("max_grid_size", max_grid_size);

  plot_int = -1 ! default to no plotfiles
  call pp%query("plot_int", plot_int);
  
  call amrex_parmparse_destroy(pp)

  ! Define a single box covering the domain
  domain = amrex_box((/0,0,0/), (/n_cell-1, n_cell-1, n_cell-1/))

  ! Initialize the boxarray "ba" from the single box "bx"
  call amrex_boxarray_build(ba, domain)

  ! Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
  call ba%maxSize(max_grid_size)

  ! Build a DistributionMapping for the boxarray
  call amrex_distromap_build(dm, ba)
  
  ! This defines a amrex_geometry object.
  call amrex_geometry_build(geom, domain)

  ! Build data multifabs
  call amrex_multifab_build(new_phi, ba, dm, ncomp, nghost)
  call amrex_multifab_build(old_phi, ba, dm, ncomp, nghost)

  call amrex_distromap_destroy(dm)
  call amrex_boxarray_destroy(ba)

  ! Intialize data
  call init_phi(new_phi, geom)

  istep = 0
  time = 0.d0

  ! choose a time step with a diffusive CFL of 0.9

#if (AMREX_SPACEDIM==1)
  coeff = 1./(geom%dx(1)*geom%dx(1))
#elif (AMREX_SPACEDIM==2)
  coeff = 1./(geom%dx(1)*geom%dx(1)) + 1./(geom%dx(2)*geom%dx(2))
#else
  coeff = 1./(geom%dx(1)*geom%dx(1)) + 1./(geom%dx(2)*geom%dx(2)) + 1./(geom%dx(3)*geom%dx(3))
#endif

  dt = 0.9d0/(2.d0*coeff)

  do istep = 1, nsteps

     if ( amrex_parallel_IOProcessor() ) then
        print*,'Advancing time step',istep,'with dt=',dt
     end if

     ! Swap the guts of multifabs so we don't have to allocate and de-allocate data
     call amrex_multifab_swap(new_phi, old_phi)

     ! advance phi
     call advance(old_phi, new_phi, geom, dt)

     time = time + dt

  end do

  call amrex_multifab_destroy(new_phi)
  call amrex_multifab_destroy(old_phi)

  call amrex_geometry_destroy(geom)

  call amrex_finalize()

end program main

