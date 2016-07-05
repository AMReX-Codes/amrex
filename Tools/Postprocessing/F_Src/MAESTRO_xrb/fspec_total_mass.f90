! This routine sums the total mass on the 2d grid for all species.
!
! This is a modification of the tutorial/fspeciesmass.f90 routine.

program fspec_total_mass

  use bl_space, only: MAX_SPACEDIM
  use bl_constants_module, only: ZERO, HALF
  use bl_IO_module
  use plotfile_module
  use network

  implicit none

  ! pf is a plotfile object that will contain both the meta-data
  ! describing the layout of the boxes that make up the AMR levels and
  ! the actual data itself (although, not necessarily all of the data
  ! at any one time, to save memory).
  type(plotfile) pf

  integer :: unit
  integer :: f, i, j, ii, jj, c
  real(kind=dp_t) :: xx, yy
  integer :: rr, r1

  real(kind=dp_t) :: dx(MAX_SPACEDIM)

  ! the pointer p will point to the data for a single box in the
  ! computational domain.  Regardless of the dimensionality of the
  ! simulation, p will always have 4 indices (3 space + 1 component).
  real(kind=dp_t), pointer :: p(:,:,:,:)

  real(kind=dp_t) :: spec_mass(nspec)
  real(kind=dp_t) :: dV

  integer :: dens_comp, X_comp(nspec), rhoX_comp(nspec)

  ! imask is an array that will be used to determine if we've already
  ! made use of data at a particular physical location (e.g. does a
  ! fine patch cover a coarser region -- if so, we want to use the
  ! finer-gridded data).  It should have the same dimensionality as
  ! the dataset.
  logical, allocatable :: imask(:,:)

  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  integer :: narg, farg
  character(len=256) :: fname, spec_name

  call network_init()

  unit = unit_new()

  narg = command_argument_count()
  farg = 1

  print *, "# time ", short_spec_names

  !---------------------------------------------------------------------------
  ! loop over the plotfiles for processing
  !---------------------------------------------------------------------------
  do f = farg, narg

     ! get the name of the plotfile
     call get_command_argument(f, value = fname)

     !------------------------------------------------------------------------
     ! build the plotfile object that contains the data describing
     ! the plotfile
     !------------------------------------------------------------------------
     call build(pf, fname, unit)

     ! make sure we are 2-d
     if (pf%dim /= 2) then
        print *, "ERROR: not a 2-d dataset"
        stop
     endif

     !------------------------------------------------------------------------
     ! figure out the variable indices from the plotfile
     !------------------------------------------------------------------------

     ! MAESTRO and CASTRO store the species information differently.  Look
     ! for either rho_X(spec_name) or X(spec_name) and the density.
  
     ! density
     dens_comp = plotfile_var_index(pf, "density")
  

     do c = 1, nspec
        spec_name = short_spec_names(c)
        ! species -- first try the Maestro naming convention
        X_comp(c) = plotfile_var_index(pf, "X(" // trim(spec_name) // ")")
        ! Castro naming convection
        rhoX_comp = plotfile_var_index(pf, "rho_" // trim(spec_name))
     enddo

     ! make sure we have either rho * X or (rho and X)
     if (any(rhoX_comp < 0) .and. (dens_comp < 0 .and. any(X_comp < 0))) then
        print *, "ERROR: variables not defined"
        stop
     endif

     !------------------------------------------------------------------------
     ! get the domain information and initialize the mask
     !------------------------------------------------------------------------

     ! get dx for the coarse level.  
     dx = plotfile_get_dx(pf, 1)


     ! get the index bounds for the finest level.  
     ! Note, lo and hi are ZERO-based indicies
     flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
     fhi = upb(plotfile_get_pd_box(pf, pf%flevel))


     ! imask is a flag for each zone at the finest possible level.
     ! imask = .true. means that we have not output data at this
     ! physical location from any level of data.  Once we output data
     ! at this location, imask is set to .false.  As we loop over
     ! levels, we will compare to the finest level index space to
     ! determine if we've already output here
     allocate(imask(flo(1):fhi(1),flo(2):fhi(2)))
     imask(:,:) = .true.


     !------------------------------------------------------------------------
     ! loop over the data
     !------------------------------------------------------------------------

     ! We start at the finest grid, and if we haven't already used
     ! data from that physical location (according to imask), add to
     ! the species mass
     spec_mass = ZERO


     ! r1 is the refinement factor between the current level grid
     ! spacing and the FINEST level
     r1  = 1

     do i = pf%flevel, 1, -1

        ! rr is the factor between the COARSEST level grid spacing and
        ! the current level.  Here we are assuming the same refinement
        ! in each coordinate direction.
        rr = product(pf%refrat(1:i-1,1))

        ! loop over all the boxes at this level of refinment
        do j = 1, nboxes(pf, i)
        
           ! read in the data 1 patch at a time
           call fab_bind(pf, i, j)

           ! get the integer bounds of the current box, in terms of this
           ! level's index space
           lo = lwb(get_box(pf, i, j))
           hi = upb(get_box(pf, i, j))

           ! get a pointer to the current patch
           p => dataptr(pf, i, j)

        
           ! loop over all of the zones in the patch.  Here, xx is the
           ! x-coordinate of the zone center and yy is the y-coordinate of
           ! the zone center.
           do jj = lo(2), hi(2)
              yy = (jj + HALF)*dx(2)/rr
                 
              do ii = lo(1), hi(1)
                 xx = (ii + HALF)*dx(1)/rr

                 
                 ! Convert the cell-centered indices at the current
                 ! level into the corresponding RANGE on the finest
                 ! level, and test if we've used data in any of those
                 ! locations.  If we haven't then we use this level's
                 ! data.
                 if ( any(imask(ii*r1:(ii+1)*r1-1, &
                                jj*r1:(jj+1)*r1-1) ) ) then

                    ! weight the zone's data by its size
                    dV = (dx(1)/rr)*(dx(2)/rr)

                    ! species mass is { rho * X * dV } 
                    if (all(rhoX_comp > 0)) then
                       do c = 1,nspec
                          spec_mass(c) = spec_mass(c) + &
                               p(ii,jj,1,rhoX_comp(c)) * dV
                       enddo
                    else
                       do c =1,nspec
                          spec_mass(c) = spec_mass(c) + &
                               p(ii,jj,1,dens_comp) * &
                               p(ii,jj,1,X_comp(c)) * dV
                       enddo
                    endif


                    ! mark this range of the domain as used in the mask
                    imask(ii*r1:(ii+1)*r1-1, &
                          jj*r1:(jj+1)*r1-1) = .false.
                 end if

              end do
           enddo

           ! delete this box's data to free memory
           call fab_unbind(pf, i, j)
                
        end do

        ! adjust r1 (refinement factor between the current level grid
        ! spacing and the FINEST level) for the next lowest level
        if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
     end do

     !------------------------------------------------------------------------
     ! output the current file's data
     !------------------------------------------------------------------------
     print *, pf%tm, spec_mass

     !------------------------------------------------------------------------
     ! clean-up
     !------------------------------------------------------------------------
     call destroy(pf)
     deallocate(imask)

  enddo  ! end loop over files

end program fspec_total_mass
