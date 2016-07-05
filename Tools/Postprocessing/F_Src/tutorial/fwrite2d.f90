! A simple routine to read in a plotfile and write it out again by

!
! This particular routine is written for 2-d datasets only.  Extending
! it to 3-d (or 1-d) is straightfoward.

program fwrite2d

  use bl_space, only: MAX_SPACEDIM
  use bl_constants_module, only: ZERO, HALF
  use bl_IO_module
  use bl_error_module
  use plotfile_module
  use multifab_module

  implicit none

  ! pf is a plotfile object that will contain both the meta-data
  ! describing the layout of the boxes that make up the AMR levels and
  ! the actual data itself (although, not necessarily all of the data
  ! at any one time, to save memory).
  type(plotfile) pf

  type(multifab), allocatable :: mf_array(:)
  type(     fab)              :: fb
  type(layout)                :: la
  type(boxarray)              :: ba
  type(list_box)              :: bl
  type(box)                   :: bx,pd
  integer       , allocatable :: ref_ratio(:)
  real(dp_t)                  :: prob_lo(2)
  real(dp_t)                  :: prob_hi(2)
  real(dp_t)                  :: time

  integer :: unit
  integer :: f, i, j

  real(kind=dp_t) :: dx(MAX_SPACEDIM)

  ! the pointer p will point to the data for a single box in the
  ! computational domain.  Regardless of the dimensionality of the
  ! simulation, p will always have 4 indices (3 space + 1 component).
  real(kind=dp_t), pointer ::  p(:,:,:,:)
  real(kind=dp_t), pointer :: mp(:,:,:,:)

  integer :: dens_comp, X_comp, rhoX_comp
  integer :: nlevs

  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  integer :: narg, farg
  character(len=256) :: fname, spec_name

  ! This is the directory name we will write the multifab to
  character(len=6) :: dirname

  ! These are the names of the variables we will write to the new plotfile
  character(len=20), allocatable :: plot_names(:)

  dirname = 'newplt'

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


  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('--species')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) spec_name

     case default
        exit

     end select
     farg = farg + 1
  end do

  allocate(plot_names(1))
  plot_names(1) = spec_name

  ! any remaining commandline options (farg to narg) are the plotfiles we 
  ! want to process

  if (farg > narg) then
     print *, " "
     print *, "Dump out the mass of a particular species on the grid"
     print *, "Works only with 2-d datasets."
     print *, " "
     print *, "usage:"
     print *, "   fspeciesmass [--species 'X'] plotfile"
     print *, " "
     print *, "args:"
     print *, "   [--species 'X']   use the species named 'X'"
     print *, " "
     stop
  end if

  !---------------------------------------------------------------------------
  ! Only process one plotfile
  !---------------------------------------------------------------------------
  f = farg

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

  ! Figure out how many levels in this plotfile
  nlevs = plotfile_nlevels(pf)

  allocate(mf_array(nlevs))
  allocate(ref_ratio(nlevs-1))

  print *,'NLEVS ',nlevs
  print *,'PF%FLEVEL ',pf%flevel

  if (nlevs .ne. pf%flevel) &
      call bl_error('nlevs not equal to flevel')

  ! get the index bounds for the finest level.  
  ! Note, lo and hi are ZERO-based indicies
  flo(1:pf%dim) = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi(1:pf%dim) = upb(plotfile_get_pd_box(pf, pf%flevel))


  !------------------------------------------------------------------------
  ! Build the new multifab
  !------------------------------------------------------------------------

  ! First define ref_ratio
  do i = 1, pf%flevel-1
     ref_ratio(i) = pf%refrat(i,1)
  end do

  ! Now loop over levels and grids to define the boxes so we can build the multifab
  do i = 1, pf%flevel

     ! loop over all the boxes at this level of refinment
     do j = 1, nboxes(pf, i)
        call push_back(bl,get_box(pf,i,j))
     end do

     call build(ba,bl)
     call layout_build_ba(la,ba,plotfile_get_pd_box(pf,1))
     call print(ba,'BA')

     ! Destroy this list and boxarray so we can start over at the new level
     call destroy(bl)
     call destroy(ba)
     

     ! Create a new multifab with 0 ghost cells and 1 component
     call multifab_build(mf_array(i),la,1,0)

  end do

  !------------------------------------------------------------------------
  ! loop over the data to "fill" the multifab
  !------------------------------------------------------------------------

  do i = pf%flevel, 1, -1

        ! loop over all the boxes at this level of refinment
        do j = 1, nboxes(pf, i)
        
           ! read in the data 1 patch at a time
           call fab_bind(pf, i, j)

           ! get the integer bounds of the current box, in terms of this
           ! level's index space
           lo(1:pf%dim) = lwb(get_box(pf, i, j))
           hi(1:pf%dim) = upb(get_box(pf, i, j))

           ! get a pointer to the current patch
           p => dataptr(pf, i, j)

           ! get a pointer to the multifab's same-sized box
           mp => dataptr(mf_array(i), j, get_box(mf_array(i),j))

           ! "fill" the data in the multifab by assigning the pointer
           mp(:,:,1,1) = p(:,:,1,dens_comp)
        
           ! delete this box's data to free memory
           call fab_unbind(pf, i, j)
                
        end do

  end do

  ! Define the problem domain at the coarsest level
  pd = plotfile_get_pd_box(pf,1)

  ! Define dx at the coarsest level
  dx(1:pf%dim) = plotfile_get_dx(pf, 1)

  ! Define problo and probhi at the coarsest level
  prob_lo = pf%plo
  prob_hi = pf%phi

  ! Get time from pf
  time = plotfile_time(pf)

  call fabio_ml_multifab_write_d(mf_array, ref_ratio, dirname, plot_names, &
                                 pd, prob_lo, prob_hi, time, dx(1:pf%dim))


  !------------------------------------------------------------------------
  ! clean-up
  !------------------------------------------------------------------------
  call destroy(pf)

end program fwrite2d
