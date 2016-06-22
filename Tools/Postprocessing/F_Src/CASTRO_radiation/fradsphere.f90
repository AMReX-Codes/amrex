! Print out the radiation quantities at a specified distance from the
! origin, for a 1-d CASTRO run.  This is geared toward the radiating 
! sphere problem.

program fradsphere

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
  integer :: ung
  integer :: cnt, max_points, nvs

  real(kind=dp_t) :: dx(MAX_SPACEDIM)
  real(kind=dp_t), pointer :: p(:,:,:,:)
  real(kind=dp_t), allocatable :: sv(:,:)
  integer, allocatable :: isv(:)
  logical, allocatable :: imask(:)
  integer :: nspec, lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)
  integer :: dim

  character(len=256) :: pltfile
  character(len=256) :: groupfile
  integer :: indslsh
  integer :: iloc, jloc, kloc
  integer :: narg, farg
  character(len=256) :: fname

  character(len=1) :: dirstr

  real(kind=dp_t) :: rmin, rmax
  real(kind=dp_t) :: radius

  integer :: ngroups
  real(kind=dp_t), allocatable :: nu_groups(:), dnu_groups(:)

  integer :: irad_begin
  integer :: idx_obs

  character(len=256) :: header_line
  integer :: ipos

  unit = unit_new()
  ung =  unit_new()

  pltfile  = ''
  groupfile = "group_structure.dat"

  radius = 0.d0

  narg = command_argument_count()

  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-p', '--pltfile')
        farg = farg + 1
        call get_command_argument(farg, value = pltfile)

     case ('-g', '--groupfile')
        farg = farg + 1
        call get_command_argument(farg, value = groupfile)

     case ('-r', '--radius')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) radius

     case default
        exit

     end select
     farg = farg + 1
  end do

  if ( len_trim(pltfile) == 0 ) then
     print *, " "
     print *, "Print out the radiation quantities at a specified distance from"
     print *, "the origin.  This is written for the 1-d radiating sphere problem."
     print *, " "
     print *, "./fradsphere -p plotfile -r radius -g groupfile"
     print *, " "
     print *, "Here groupfile is the file containing the group structure information"
     print *, "as output by Castro (usually group_structure.dat)."
     stop
  end if


  print *, 'pltfile   = "', trim(pltfile), '"'


  call build(pf, pltfile, unit)

  nvs = pf%nvars

  ! find the index of the first radiation group's energy
  irad_begin = -1
  irad_begin = plotfile_var_index(pf, "rad0")

  if (irad_begin < 0) then
     call bl_error("ERROR: radiation field not found in plotfile")
  endif

  ! get the index bounds and dx for the coarse level.  Note, lo and hi are
  ! ZERO based indicies
  lo(:) = ZERO
  hi(:) = ZERO

  lo = lwb(plotfile_get_pd_box(pf, 1))
  hi = upb(plotfile_get_pd_box(pf, 1))

  dx = plotfile_get_dx(pf, 1)


  ! get the physical domain size
  rmin = pf%plo(1)
  rmax = pf%phi(1)

  if (radius < rmin .or. radius > rmax) then
     call bl_error("ERROR: specified observer radius outside of domain")
  endif

  print *, 'rmin = ', rmin
  print *, 'rmax = ', rmax

  dim = pf%dim
  if (dim /= 1) call bl_error("ERROR: fradsphere only works for dim = 1")

  ! get the index bounds for the finest level
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

  cnt = 0

  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here

  allocate(imask(flo(1):fhi(1)))
  max_points = fhi(1) - flo(1) + 1

  imask = .true.

  ! allocate storage for the data -- we will assume that we are
  ! completely refined, and allocate enough storage for that,
  ! although, we won't really need all of this
  !
  ! note: sv(:,1) will be the coordinate information.  
  ! the variables will be stored in sv(:,2:nvs+1)
  allocate(sv(max_points,nvs+1), isv(max_points))

  ! loop over the data, starting at the finest grid, and if we haven't
  ! already store data in that grid location (according to imask),
  ! store it.  We'll put it in the correct order later.

  cnt = 0

  ! r1 is the factor between the current level grid spacing and the
  ! FINEST level
  r1  = 1

  do i = pf%flevel, 1, -1

     ! rr is the factor between the COARSEST level grid spacing and
     ! the current level
     rr = product(pf%refrat(1:i-1,1))

     do j = 1, nboxes(pf, i)

        ! read in the data 1 patch at a time
        call fab_bind(pf, i, j)

        lo(:) = 1
        hi(:) = 1
        lo = lwb(get_box(pf, i, j))
        hi = upb(get_box(pf, i, j))

        p => dataptr(pf, i, j)

        do ii = lo(1), hi(1)
           if ( any(imask(ii*r1:(ii+1)*r1-1) ) ) then
              cnt = cnt + 1
                    
              sv(cnt,1) = rmin + (ii + HALF)*dx(1)/rr
              sv(cnt,2:) = p(ii,1,1,:)
                    
              imask(ii*r1:(ii+1)*r1-1) = .false.
           end if
        end do

        call fab_unbind(pf, i, j)
     end do

     ! adjust r1 for the next lowest level
     if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
  end do


  ! sort the data based on the coordinates
  call sort(sv(1:cnt,1),isv(1:cnt))

  dirstr = "r"


  ! open up the group file and read in the group information
  open(unit=ung,file=groupfile)

  ! read in the number of groups
  read(ung,fmt='(a256)') header_line
  ipos = index(header_line, "=") + 1
  read (header_line(ipos:),*) ngroups

  allocate (nu_groups(ngroups), dnu_groups(ngroups))

  ! skip the next header line
  read(ung,*) header_line

  ! read in the group centers and weights
  do i = 1, ngroups
     read(ung,*) nu_groups(i), dnu_groups(i)
  enddo

  close(ung)


1000 format("#",100(a24,1x))
1001 format(1x, 100(g24.12,1x))


  ! find the index corresponding to the desired observer radius
  idx_obs = -1
  do i = 1,cnt-1
     if (radius >= sv(isv(i),1) .and. radius < sv(isv(i+1),1)) then
        idx_obs = i
        exit
     endif
  enddo

  if (idx_obs == -1) call bl_error("ERROR: radius not found in domain")

  ! output all the radiation energies

  write(*,1000) "group name", "group center energy", &
       "E_rad(nu)*dnu (erg/cm^3)", "E_rad(nu) (erg/cm^3/Hz)"
  do i = 1, ngroups
     write (*,1001) pf%names(irad_begin-1+i), &               
                    nu_groups(i), &  
                    sv(isv(idx_obs),1+irad_begin-1+i), &      
                    sv(isv(idx_obs),1+irad_begin-1+i)/dnu_groups(i)
  enddo

  call destroy(pf)

end program fradsphere
