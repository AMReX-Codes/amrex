! A simple routine to add up the total mass of species X in a certain
! region on the grid.
!
! This routine is a modification of the routine in
! ../tutorial/fspeciesmass2d.f90
!

program fspeciesmass

  use bl_space, only: MAX_SPACEDIM
  use bl_constants_module, only: ZERO, HALF
  use bl_IO_module
  use plotfile_module

  implicit none

  type(plotfile) pf

  integer :: unit
  integer :: f, i, j, ii, jj
  real(kind=dp_t) :: xx, yy
  integer :: rr, r1

  real(kind=dp_t), dimension(MAX_SPACEDIM) :: dx, plo, phi
  real(kind=dp_t), pointer :: p(:,:,:,:)
  real(kind=dp_t) :: spec_mass, dV
  real(kind=dp_t) :: lo_coord, hi_coord, ulo_coord, uhi_coord
  logical :: lo_check, hi_check
  integer :: dim

  integer :: dens_comp, X_comp, rhoX_comp

  logical, allocatable :: imask(:,:)

  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  integer :: narg, farg
  character(len=256) :: fname, spec_name


  unit = unit_new()

  !---------------------------------------------------------------------------
  ! parse the commandline optionns using the f2kcli module
  !---------------------------------------------------------------------------

  narg = command_argument_count()

  ! first process any options -- loop over the arguments and if it
  ! matches one of the commandline options, interpret its value.

  ! we expect to have the commandline options come first, followed by
  ! any plotfiles we wish to process

  ! set the defaults
  spec_name = "C12"
  ulo_coord = -1
  uhi_coord = -1
  lo_check = .false.
  hi_check = .false.

  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-s','--species')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) spec_name

     case ('-a','--above')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) ulo_coord
        lo_check = .true.

     case ('-b','--below')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) uhi_coord
        hi_check = .true.

     case default
        exit

     end select
     farg = farg + 1
  end do

  ! any remaining commandline options (farg to narg) are the plotfiles we 
  ! want to process

  if (farg > narg) then
     print *, " "
     print *, "Dump out the mass of a particular species in a section of the"
     print *, 'grid.'
     print *, " "
     print *, "usage:"
     print *, "   fspeciesmass [--species X -a lo_coord -b hi_coord] plotfile"
     print *, " "
     print *, "args:"
     print *, "   [-s|--species 'X']   use the species named 'X'. defaults to"
     print *, "                        'C12'"
     print *, '   [-a|--above lo_coord] calculate the mass only for r >= lo_coord'
     print *, '   [-b|--below hi_coord] calculate the mass only for r <= hi_coord'
     print *, '   the default behaviour is to use the full grid.  the "-a" and'
     print *, '   "-b" options can be combined'
     print *, " "
     stop
  end if

  ! sanity check
  if (ulo_coord >= 0 .and. uhi_coord >= 0 .and. ulo_coord > uhi_coord) &
       call bl_error("Error: specified lo_coord > hi_coord")
  if (lo_check .and. ulo_coord < 0) &
       call bl_error("Error: specified lo_coord < 0")
  if (hi_check .and. uhi_coord < 0) &
       call bl_error("Error: specified hi_coord < 0")

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

     dim = plotfile_dim(pf)

     ! make sure we are 2-d
     if (dim/= 2) then
        print *, "ERROR: not a 2-d dataset"
        stop
     endif

     plo = pf%plo
     phi = pf%phi

     ! sanity checks
     lo_coord = plo(dim)
     hi_coord = phi(dim)
     if (ulo_coord >=0) lo_coord = ulo_coord
     if (uhi_coord >=0) hi_coord = uhi_coord

     if (lo_coord < plo(dim) .or. lo_coord > phi(dim)) &
          call bl_error("Error: specified lo_coord outside of domain")
     if (hi_coord < plo(dim) .or. hi_coord > phi(dim)) &
          call bl_error("Error: specified hi_coord outside of domain")

     !------------------------------------------------------------------------
     ! figure out the variable indices from the plotfile
     !------------------------------------------------------------------------

     ! density
     dens_comp = plotfile_var_index(pf, "density")
  
     ! species -- first try the Maestro naming convention
     X_comp = plotfile_var_index(pf, "X(" // trim(spec_name) // ")")

     ! Castro naming convection
     rhoX_comp = plotfile_var_index(pf, "rho_" // trim(spec_name))
     
     ! make sure we have either rho * X or (rho and X)
     if (rhoX_comp < 0 .and. (dens_comp < 0 .and. X_comp < 0)) then
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
              yy = plo(2) + (jj + HALF)*dx(2)/rr

              if (yy < lo_coord .or. yy > hi_coord) cycle
                 
              do ii = lo(1), hi(1)
                 xx = plo(1) + (ii + HALF)*dx(1)/rr

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
                    if (rhoX_comp > 0) then
                       spec_mass = spec_mass + &
                            p(ii,jj,1,rhoX_comp) * dV
                    else
                       spec_mass = spec_mass + &
                            p(ii,jj,1,dens_comp) * p(ii,jj,1,X_comp) * dV
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

end program fspeciesmass
