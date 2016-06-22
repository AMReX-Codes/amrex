! compute the amplitude of the LD cusps


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

  integer :: index_comp(1)
  character (len=64) :: compname
  logical :: found

  integer :: jmin, jmax
  real(kind=dp_t) :: cmaxthres, cminthres, dy

  real(kind=dp_t), allocatable :: c_fab(:,:,:), c_avg(:), c_min(:), c_max(:)
  type(plotfile) :: pf
  integer narg
  character(len=128) :: phome
  character(len=256) ::  fname

  integer unit
  integer ii, jj

  logical :: verbose
  integer :: ntime
  integer :: max_level

  cminthres = 0.05d0
  cmaxthres = 0.45d0

  narg = command_argument_count()
  call get_environment_variable("HOME", value = phome)

  verbose = .FALSE.

  max_level = 1

  farg = 1

  ! parse the runtime parameters
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('--max_level')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) max_level

     case default
        exit
     end select
     farg = farg + 1
  end do

  if ( farg > narg )      return

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
  compname = "Y(C12)"

  found = .false.
  do i = 1, pf%nvars
     if (compname == pf%names(i)) then
        index_comp = i
        found = .true.
        exit
     endif
  enddo

  if (.NOT. found) then
     print *, 'ERROR: component = ', compname, ' not found'
     print *, pf%nvars, ' components stored'
     print *, 'available components are: '
     do i = 1, pf%nvars
        print *, pf%names(i)
     enddo
        
     stop
  endif


  print *, ' Files  =', narg - farg + 1
  print *, ' '

  ! allocate storage for the data on the finest mesh we are using
  call get_command_argument(farg, value = fname)
  call build(pf, fname, unit)

  plo = lwb(plotfile_get_pd_box(pf, max_level))
  phi = upb(plotfile_get_pd_box(pf, max_level))
  allocate(c_fab(plo(1):phi(1), plo(2):phi(2),1))
  allocate(c_avg(plo(2):phi(2)))
  allocate(c_min(plo(2):phi(2)))
  allocate(c_max(plo(2):phi(2)))

  ! loop over all the files, build a single FAB, an then loop over the 
  ! individual components we wish to store, and store them one by one
  dy = (pf%phi(2) - pf%plo(2))/(phi(2)-plo(2)+1)
  print *, 'dy = ', dy

  do f = 1, ntime
     call get_command_argument(f + farg - 1, value = fname)
     !     print *, 'File: ', trim(fname)
     if ( f /= 1 ) call build(pf, fname, unit)

     call blow_out_to_fab(c_fab, plo, pf, index_comp, max_level)

     ! do the lateral averaging
     
     do j = plo(2), phi(2)
        c_avg(j) = sum(c_fab(:,j,1))/(phi(1) - plo(1) + 1)
        c_min(j) = minval(c_fab(:,j,1))
        c_max(j) = maxval(c_fab(:,j,1))
     enddo

     do j = plo(2), phi(2)
        if (c_min(j) < cmaxthres) then
           jmax = j
           exit
        endif
     enddo

     do j = plo(2), phi(2)
        if (c_max(j) < cmaxthres) then
           jmin = j
           exit
        endif
     enddo
     

     print *, pf%tm, dy*jmin + pf%plo(2) , dy*jmax + pf%plo(2)

     call destroy(pf)
  end do

end subroutine process_snapshots
