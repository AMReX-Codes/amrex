! Process a 2-d gaussian radiation pulse

program fgaussianpulse

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module

  implicit none

  type(plotfile) pf
  integer :: unit
  integer :: i, j, ii, jj
  real(kind=dp_t) :: xx, yy
  integer :: rr, r1
  integer :: uno

  integer :: nbins
  real(kind=dp_t), allocatable :: r(:), rl(:)
  real(kind=dp_t) :: maxdist, x_maxdist, y_maxdist
  real(kind=dp_t) :: xctr, yctr

  real(kind=dp_t) :: dx(MAX_SPACEDIM)
  real(kind=dp_t) :: dx_fine

  real(kind=dp_t) :: r_zone
  integer :: index

  real(kind=dp_t), pointer :: p(:,:,:,:)

  integer, allocatable :: ncount(:)
  real(kind=dp_t), allocatable :: rad_bin(:)

  integer :: rad_comp

  logical, allocatable :: imask(:,:)
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  character(len=256) :: slicefile
  character(len=256) :: pltfile

  integer :: narg, farg
  character(len=256) :: fname

  unit = unit_new()
  uno =  unit_new()


  ! set the defaults
  slicefile = ''
  pltfile  = ''

  xctr = HALF
  yctr = HALF


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

     case ('--xctr')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) xctr

     case ('--yctr')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) yctr

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
  
  ! radiation
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


  ! get the index bounds for the finest level
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))


  ! compute the size of the radially-binned array -- we'll do it to
  ! the furtherest corner of the domain
  x_maxdist = max(abs(pf%phi(1) - xctr), abs(pf%plo(1) - xctr))
  y_maxdist = max(abs(pf%phi(2) - yctr), abs(pf%plo(2) - yctr))
  
  maxdist = sqrt(x_maxdist**2 + y_maxdist**2)

  dx_fine = minval(plotfile_get_dx(pf, pf%flevel))
  nbins = int(maxdist/dx_fine)

  allocate(r(0:nbins-1))
  allocate(rl(0:nbins))

  do i = 0, nbins-1
     r(i) = (dble(i) + HALF)*dx_fine
     rl(i) = dble(i)*dx_fine
  enddo
  rl(nbins) = dble(nbins)*dx_fine


  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here
  allocate(imask(flo(1):fhi(1),flo(2):fhi(2)))
  imask(:,:) = .true.


  ! allocate storage for the data 
  allocate(rad_bin(0:nbins-1))
  allocate(  ncount(0:nbins-1))

  ncount(:) = 0
  rad_bin(:) = ZERO

  ! loop over the data, starting at the finest grid, and if we haven't
  ! already store data in that grid location (according to imask),
  ! store it.  

  ! r1 is the factor between the current level grid spacing and the
  ! FINEST level
  r1  = 1

  do i = pf%flevel, 1, -1

     ! rr is the factor between the COARSEST level grid spacing and
     ! the current level
     rr = product(pf%refrat(1:i-1,1))

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
        do jj = lo(2), hi(2)
           yy = (jj + HALF)*dx(2)/rr

           do ii = lo(1), hi(1)
              xx = (ii + HALF)*dx(1)/rr

              if ( any(imask(ii*r1:(ii+1)*r1-1, &
                             jj*r1:(jj+1)*r1-1) ) ) then
                 
                 r_zone = sqrt((xx-xctr)**2 + (yy-yctr)**2)

                 index = r_zone/dx_fine

                 ! weight the zone's data by its size
                 rad_bin(index) = rad_bin(index) + &
                      p(ii,jj,1,rad_comp)*r1**2

                 ncount(index) = ncount(index) + r1**2

                 imask(ii*r1:(ii+1)*r1-1, &
                       jj*r1:(jj+1)*r1-1) = .false.
                 
              end if

           end do
        enddo

     end do

     ! adjust r1 for the next lowest level
     if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
  end do

  ! normalize
  do i = 0, nbins-1
     if (ncount(i) /= 0) then
        rad_bin(i) = rad_bin(i)/ncount(i)
     endif
  enddo

1000 format("#",100(a24,1x))
1001 format(1x, 100(g24.13,1x))

  ! slicefile
  open(unit=uno, file=slicefile, status = 'replace')

  ! write the header
  write(uno,1000) "x", "rad"

  ! write the data in columns
  do i = 0, nbins-1
     ! Use this to protect against a number being xx.e-100 
     !   which will print without the "e"
     if (abs(rad_bin(i)) .lt. 1.d-99) rad_bin(i) = 0.d0
     write(uno,1001) r(i), rad_bin(i)
  end do

  close(unit=uno)

  do i = 1, pf%flevel
     call fab_unbind_level(pf, i)
  end do

  call destroy(pf)

end program fgaussianpulse
