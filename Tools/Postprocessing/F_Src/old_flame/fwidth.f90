! compute the width of the RT mixed region.
! here, we compute the average of x*(.5-x), where x is the
! carbon mass fraction.

! in this version, we can cut the domain into strips when doing the
! averaging, so we can deal with very large domains without running
! out of memory.  `./fwidth.exe --cut 4 plt*' will use 4 strips (cut
! in the x-direction) and build the finely gridded subdomain one 
! strip at a time when doing the averaging.

program fdata
  implicit none

  call process_snapshots

end program fdata

subroutine process_snapshots
  use plotfile_module
  use filler_module
! use utility_module
  use bl_IO_module

  implicit none

  integer, parameter :: NCOLOR = 256
  integer :: plo(2), phi(2)
  integer f, i, j, farg
  integer :: n
  integer :: index_store

  integer, parameter :: nc = 1
  integer :: index_comp(nc)
  character (len=64) :: compname(nc)
  logical :: found

  integer :: jmin, jmax, j2min, j2max
  real(kind=dp_t) :: c2thres, dy, cminthres, cmaxthres

  real(kind=dp_t), allocatable :: c_fab(:,:,:), c_avg(:), c2_avg(:)
  type(plotfile) :: pf
  integer narg
  character(len=128) :: phome
  character(len=256) ::  fname

  integer unit
  integer ii, jj

  integer :: nx

  real(kind=dp_t) :: dx
  integer, dimension(2) :: tlo(2), thi(2)

  logical :: verbose
  integer :: ntime
  integer :: max_level

  integer :: cut

  c2thres = 6.25d-3

  cminthres = 0.05d0
  cmaxthres = 0.45d0

  narg = command_argument_count()
  call get_environment_variable("HOME", value = phome)

  verbose = .FALSE.

  max_level = 1
  cut = 1

  farg = 1

  ! parse the runtime parameters
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('--max_level')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) max_level

     case ('--cut')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) cut

     case default
        exit
     end select
     farg = farg + 1
  end do

  if ( farg > narg )      return

  print *, '#  <x*(0.5-x)> threshold = ', c2thres
  print *, '#  <x> thresholds (min, max) = ', cminthres, cmaxthres

  unit = unit_new()

  ! ntime is the number of files to loop over
  ntime = narg - farg  + 1

  call get_command_argument(farg, value = fname)
  call build(pf, fname, unit)

  ! assume that all the files have the same components stored.  Convert
  ! components that were specified via string names into their 
  ! corresponding integer indices in the comps() array.  Don't worry
  ! about double counting for now.

  ! the names of the pf%nvars components are stored in pf%names()
  compname(1) = "Y(C12)"

  do j=1,nc
     found = .false.
     do i = 1, pf%nvars
        if (compname(j) == pf%names(i)) then
           index_comp(j) = i
           found = .true.
           exit
        endif
     enddo

     if (.NOT. found) then
        print *, 'ERROR: component = ', compname(j), ' not found'
        print *, pf%nvars, ' components stored'
        print *, 'available components are: '
        do i = 1, pf%nvars
           print *, pf%names(i)
        enddo
        
        stop
     endif
  enddo

  ! allocate storage for the data on the finest mesh we are using
  call get_command_argument(farg, value = fname)
  call build(pf, fname, unit)

  plo = lwb(plotfile_get_pd_box(pf, max_level))
  phi = upb(plotfile_get_pd_box(pf, max_level))

  allocate(c_avg(plo(2):phi(2)))
  allocate(c2_avg(plo(2):phi(2)))

  ! loop over all the files, build a single FAB, an then loop over the 
  ! individual components we wish to store, and store them one by one
  dy = (pf%phi(2) - pf%plo(2))/(phi(2)-plo(2)+1)


  nx = phi(1) - plo(1) + 1
  dx = nx/cut
 
  print *, '#         <x> method    <x*(0.5-x)> method'
  print *, '# time     min  max          min  max'
  
  do f = 1, ntime

     c_avg(:) = 0.0
     c2_avg(:) = 0.0

     call get_command_argument(f + farg - 1, value = fname)
     if ( f /= 1 ) call build(pf, fname, unit)

     ! we are going to blow this out to a fab in pieces, so we can
     ! handle large domains.  We cut in the x direction.
     do n = 1, cut

        ! allocate the storage for the current subdomain
        tlo(1) = (n-1)*dx
        thi(1) = n*dx -1

        if (n == cut) then
           thi(1) = nx - 1
        endif

        tlo(2) = plo(2)
        thi(2) = phi(2)
        
        allocate(c_fab(tlo(1):thi(1), tlo(2):thi(2),nc))
        call blow_out_to_sub_fab(c_fab, tlo, thi, pf, index_comp, max_level)

        ! do the lateral averaging -- don't normalize here -- we'll save
        ! that to the end
        
        do j = tlo(2), thi(2)
           c_avg(j) = c_avg(j) + sum(c_fab(:,j,1))
           c2_avg(j) = c2_avg(j) + sum(c_fab(:,j,1)*(0.5 - c_fab(:,j,1)))
        enddo

        deallocate (c_fab)

     enddo

     ! normalize
     c_avg = c_avg/(phi(1) - plo(1) + 1)
     c2_avg = c2_avg/(phi(1) - plo(1) + 1)

     ! find the width for the <x> case

     do j = plo(2), phi(2)
        if (c_avg(j) < cmaxthres) then
           jmax = j
           exit
        endif
     enddo

     do j = plo(2), phi(2)
        if (c_avg(j) < cminthres) then
           jmin = j
           exit
        endif
     enddo


     ! find the width for the <x*(0.5-x)> case

     do j = phi(2), plo(2), -1
        if (c2_avg(j) > c2thres) then
           j2max = j
           exit
        endif
     enddo

     do j = plo(2), phi(2)
        if (c2_avg(j) > c2thres) then
           j2min = j
           exit
        endif
     enddo
     
     print *, real(pf%tm), real(dy*jmax + pf%plo(2)) , &
          real(dy*jmin + pf%plo(2)), real(dy*j2min + pf%plo(2)) , &
          real(dy*j2max + pf%plo(2))

     call destroy(pf)
  end do

end subroutine process_snapshots
