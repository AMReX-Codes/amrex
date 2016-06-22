! Process a 1-d sedov problem to produce rho, u, and p as a
! function of r, for comparison to the analytic solution.

program fextract1d

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use sort_d_module

  implicit none

  type(plotfile) pf
  integer :: unit
  integer :: i, j, ii
  integer :: rr, r1
  integer :: uno

  integer :: nbins
  real(kind=dp_t), allocatable :: r(:)
  real(kind=dp_t) :: maxdist, x_maxdist, y_maxdist
  real(kind=dp_t) :: xctr, yctr

  real(kind=dp_t) :: dx(MAX_SPACEDIM)
  real(kind=dp_t) :: dx_fine

  integer :: index

  real(kind=dp_t), pointer :: p(:,:,:,:)

  real(kind=dp_t), allocatable :: dens_bin(:), vel_bin(:), pres_bin(:)
! esm [
  real(kind=dp_t), allocatable :: rad_bin(:)
! esm ]

  integer :: dens_comp, xmom_comp, pres_comp
! esm [
  integer :: rad_comp
! esm ]

  logical, allocatable :: imask(:)
  integer :: nspec, lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  character(len=256) :: slicefile
  character(len=256) :: pltfile

  logical :: lexist
  integer :: narg, farg
  character(len=256) :: fname

  unit = unit_new()
  uno =  unit_new()

  ! set the defaults
  slicefile = ''
  pltfile  = ''

  narg = command_argument_count()

  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-p', '--pltfile')
        farg = farg + 1
        call get_command_argument(farg, value = pltfile)

     case ('-s', '--slicefile')
        farg = farg + 1
        call get_command_argument(farg, value = slicefile)

     case default
        exit

     end select
     farg = farg + 1
  end do

  if ( len_trim(pltfile) == 0 .OR. len_trim(slicefile) == 0 ) then
     print *, "usage: fextract args"
     print *, "args [-p|--pltfile]   plotfile   : plot file directory (required)"
     print *, "     [-s|--slicefile] slice file : slice file          (required)"
     stop
  end if

  print *, 'pltfile   = "', trim(pltfile), '"'
  print *, 'slicefile = "', trim(slicefile), '"'

  call build(pf, pltfile, unit)

  do i = 1, pf%flevel
     call fab_bind_level(pf, i)
  end do


  ! figure out the variable indices
  
  ! density
  dens_comp = -1
  do i = 1, pf%nvars
     if (pf%names(i) == "density") then
        dens_comp = i
        exit
     endif
  enddo
  if (dens_comp < 0) call bl_error("density not found")


  ! x-momentum
  xmom_comp = -1
  do i = 1, pf%nvars
     if (pf%names(i) == "xmom") then
        xmom_comp = i
        exit
     endif
  enddo
  if (xmom_comp < 0) call bl_error("xmom not found")

  ! pressure
  pres_comp = -1
  do i = 1, pf%nvars
     if (pf%names(i) == "pressure") then
        pres_comp = i
        exit
     endif
  enddo
  if (pres_comp < 0) call bl_error("pres not found")

  ! rad
  rad_comp = -1
  do i = 1, pf%nvars
     if (pf%names(i) == "rad") then
        rad_comp = i
        exit
     endif
  enddo
  if (rad_comp < 0) call bl_error("rad not found")

  ! get the index bounds and dx for the coarse level.  Note, lo and hi are
  ! ZERO based indicies
  lo = lwb(plotfile_get_pd_box(pf, 1))
  hi = upb(plotfile_get_pd_box(pf, 1))

  dx = plotfile_get_dx(pf, 1)

! esm [
  write(6,*) 'Checkpoint 1'
  write(6,*) 'lo = ', lo
  write(6,*) 'hi = ', hi
  write(6,*) 'pf%plo = ', pf%plo
  write(6,*) 'pf%phi = ', pf%phi
! esm ]


  ! get the index bounds for the finest level
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

  ! compute the size of the radially-binned array -- we'll do it to
  ! the furtherest corner of the domain
! esm  maxdist = abs(pf%phi(1))
! esm [
  maxdist = abs(pf%phi(1) - pf%plo(1))
! esm ]
  dx_fine = minval(plotfile_get_dx(pf, pf%flevel))

  nbins = int(maxdist/dx_fine)

  allocate(r(0:nbins-1))
  do i = 0, nbins-1
! esm     r(i) = (dble(i) + HALF)*dx_fine
! esm [
     r(i) = pf%plo(1) + (dble(i) + HALF)*dx_fine
! esm ]
  enddo

  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here
  allocate(imask(flo(1):fhi(1)))
  imask(:) = .true.

  ! allocate storage for the data 
  allocate(dens_bin(0:nbins-1))
  allocate( vel_bin(0:nbins-1))
  allocate(pres_bin(0:nbins-1))
! esm [
  allocate(rad_bin(0:nbins-1))
! esm ]

  dens_bin(:) = ZERO
  vel_bin(:) = ZERO
  pres_bin(:) = ZERO
! esm [
  rad_bin(:) = ZERO
! esm ]

  ! loop over the data, starting at the finest grid, and if we haven't
  ! already store data in that grid location (according to imask),
  ! store it.  We'll put it in the correct order later.

  ! r1 is the factor between the current level grid spacing and the
  ! FINEST level
  r1  = 1

! esm [
  write(6,*) 'Checkpoint 2'
! esm ]

  do i = pf%flevel, 1, -1

     ! rr is the factor between the COARSEST level grid spacing and
     ! the current level
     rr = product(pf%refrat(1:i-1,1))

! esm [
  write(6,*) 'Checkpoint 3; i = ', i
! esm ]
     do j = 1, nboxes(pf, i)

        lo(:) = 1
        hi(:) = 1        
        lo = lwb(get_box(pf, i, j))
        hi = upb(get_box(pf, i, j))

        ! get a pointer to the current patch
        p => dataptr(pf, i, j)
        
        ! loop over all of the zones in the patch.  Here, we convert
        ! the cell-centered indices at the current level into the
        ! corresponding RANGE on the finest level, and test if we've
        ! stored data in any of those locations.  If we haven't then
        ! we store this level's data and mark that range as filled.
        do ii = lo(1), hi(1)

              if ( any(imask(ii*r1:(ii+1)*r1-1) ) ) then

                 index = ii * r1

                 dens_bin(index:index+(r1-1)) = p(ii,1,1,dens_comp)

                 vel_bin(index:index+(r1-1)) = &
                       abs(p(ii,1,1,xmom_comp)) / p(ii,1,1,dens_comp) 

                 pres_bin(index:index+(r1-1)) = p(ii,1,1,pres_comp) 

! esm [
                 rad_bin(index:index+(r1-1)) = p(ii,1,1,rad_comp) 
! esm ]

                 imask(ii*r1:(ii+1)*r1-1) = .false.
                 
              end if

        enddo
     end do

     ! adjust r1 for the next lowest level
     if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
  end do

! esm [
  write(6,*) 'Checkpoint 4'
! esm ]

1000 format("#",100(a24,1x))
1001 format(1x, 100(g24.12,1x))

  ! slicefile
  open(unit=uno, file=slicefile, status = 'replace')

  ! write the header
  write(uno,1000) "x", "density", "velocity", "pressure"

  ! write the data in columns
  do i = 0, nbins-1
     ! Use this to protect against a number being xx.e-100
     !   which will print without the "e"
     if (abs(dens_bin(i)) .lt. 1.d-99)  dens_bin(i) = 0.d0
     if (abs( vel_bin(i)) .lt. 1.d-99)   vel_bin(i) = 0.d0
     if (abs(pres_bin(i)) .lt. 1.d-99)  pres_bin(i) = 0.d0
! esm [
     if (abs( rad_bin(i)) .lt. 1.d-99)   rad_bin(i) = 0.d0
! esm ]
! esm     write(uno,1001) r(i), dens_bin(i), vel_bin(i), pres_bin(i)
! esm [
     write(uno,1001) r(i), dens_bin(i), vel_bin(i), pres_bin(i), rad_bin(i)
! esm ]
  end do

  close(unit=uno)

  do i = 1, pf%flevel
     call fab_unbind_level(pf, i)
  end do

  call destroy(pf)

  deallocate(r)

end program fextract1d
