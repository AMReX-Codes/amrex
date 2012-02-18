! dump out a Fortran unformatted binary file containing plotfile data 
! mapped onto a single, uniform grid, with resolution set to the 
! resolution of the highest AMR level we specify.
! 
! in this version, we write out a single precision array, to save space.  
! Optionally, the domain can be divided up so it is read in multiple passes,
! cast into 4-byte reals, and stored in the main array
!
! 
! currently, the following is written out:
!
! xmin, ymin, zmin, xmax, ymax, zmax
! nx, ny, nz
! time
! data
!
! xmin, ymin, ... double precision reals.
! nx, ny, and nz are integers
! time is double precision
! data is SINGLE precision -- to conserve space.


program fIDLdump3d

  use plotfile_module
  use filler_module
  use bl_parmparse_module
  use f2kcli
  use bl_IO_module

  implicit none

  integer, parameter :: NCOLOR = 256
  integer :: plo(3), phi(3), tlo(3), thi(3), rlo(3), rhi(3)
  integer f, i, farg
  integer :: n, nc, ncc, ncs
  integer :: index_store

  integer :: zmin, zmax

  integer, allocatable :: comps(:)
  character (len=64), allocatable :: compnames(:)
  character (len=64) :: compname, tempname
  logical :: found

! declare a single byte array to hold the scaled data
  real(kind=dp_t), allocatable :: c_fab(:,:,:,:)
  real, allocatable :: rData (:,:,:,:)

  type(plotfile) :: pf
  integer narg
  character(len=256) ::  fname
  character(len=256) :: ofname, outbase
  integer :: iout

  integer unit

  integer :: nx
  real (kind=dp_t) :: dx

  integer :: divide

  logical :: verbose
  integer :: indslsh, ntime
  integer :: max_level

  narg = command_argument_count()

  ! components to store can be specified either by number '-c' or name 
  ! '-cname' -- count the total number of components to store, regardless
  ! of how we requested them
  nc = pp_arg_count("-c", "--component")
  nc = nc + pp_arg_count("-cname", "--component_name")

  ncc = 0
  ncs = 0
  allocate(comps(nc))
  allocate(compnames(nc))

  verbose = .FALSE.

  max_level = 1

  iout = 0
  farg = 1

  zmin = -1
  zmax = 10420000

  divide = 1

  ! parse the runtime parameters
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-v','--verbose')
        verbose = .TRUE.

     case ('-c','--component')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        ncc = ncc + 1
        read(fname, *) comps(ncc)

     case ('-cname','--component_name')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        ncs = ncs + 1
        compnames(ncs) = trim(fname)

     case ('--max_level')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) max_level

     case ('-o')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) outbase
        iout = 1
     
     case ('-zmin')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) zmin

     case ('-zmax')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) zmax

     case ('--divide')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) divide
        
     case default
        exit
     end select
     farg = farg + 1
  end do

  if ( nc == 0 ) stop
  if ( farg > narg ) stop
  if (nc /= 1) then
     print *, 'ERROR: scaling only implemented for one component at a time'
     stop
  endif

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

  index_store = ncc+1
  do n = 1, ncs
     found = .false.
     do i = 1, pf%nvars
        if (compnames(n) == pf%names(i)) then
           comps(index_store) = i
           index_store = index_store + 1
           found = .true.
           exit
        endif
     enddo
     if (.NOT. found) then
        print *, 'ERROR: component = ', compnames(n), ' not found'
        print *, pf%nvars, ' components stored'
        print *, 'available components are: '
        do i = 1, pf%nvars
           print *, pf%names(i)
        enddo
        
        stop
     endif
  enddo

  ! now all nc components are in the comps array

  
  print *, ' Files  =', narg - farg + 1
  do n = 1, nc
     print *, '  Comp  = ', trim(pf%names(comps(n))), ' index = ', comps(n)
  end do
  print *, ' '

  ! allocate storage for the data on the finest mesh we are using
  call get_command_argument(farg, value = fname)
  call build(pf, fname, unit)

  plo = lwb(plotfile_get_pd_box(pf, max_level))
  phi = upb(plotfile_get_pd_box(pf, max_level))

  if (zmin < plo(3)) zmin = plo(3)
  if (zmax > phi(3)) zmax = phi(3)
  if (zmax < zmin) zmax = zmin+1

  tlo(1) = plo(1)
  tlo(2) = plo(2)
  tlo(3) = zmin

  thi(1) = phi(1)
  thi(2) = phi(2)
  thi(3) = zmax

  ! rlo and rhi are the bounds of the 4-byte array
  rlo = tlo
  rhi = thi

  nx = phi(1) - plo(1) + 1
  dx = nx/divide

  ! allocate the reduced precision data array
  allocate(rData(rlo(1):rhi(1), rlo(2):rhi(2), rlo(3):rhi(3), nc))
  print *, 'main array allocated'

  ! loop over all the files, build a single FAB, an then loop over the 
  ! individual components we wish to store, and store them one by one
  do f = 1, ntime
     call get_command_argument(f + farg - 1, value = fname)
     print *, 'File: ', trim(fname)
     if ( f /= 1 ) call build(pf, fname, unit)
     
     do n = 1, divide

        ! adjust the `1' limits to represent the slice we are reading in
        tlo(1) = (n-1)*dx
        thi(1) = n*dx - 1

        if (n == divide) then
           thi(1) = nx - 1
        endif

        print *, 'allocating ', tlo(1),thi(1), tlo(2),thi(2), tlo(3),thi(3)
        allocate(c_fab(tlo(1):thi(1), tlo(2):thi(2), tlo(3):thi(3), nc))
        print *, '... done'

        call blow_out_to_sub_fab_3d(c_fab, tlo, thi, pf, comps, max_level)

        ! do the casting
        rData(tlo(1):thi(1), tlo(2):thi(2), tlo(3):thi(3), nc) = &
             real(c_fab(tlo(1):thi(1), tlo(2):thi(2), tlo(3):thi(3), nc))

        deallocate(c_fab)
     enddo

     ! get basename of file
     indslsh = index(fname, '/', back = .TRUE.)

     if (iout == 0) then
        if ( indslsh /= 0 ) then
           ofname = fname(indslsh+1:)
        else
           ofname = fname
        end if
     else
        if ( indslsh /= 0 ) then
           ofname = trim(outbase) // '.' // fname(indslsh+1:)
        else
           ofname = trim(outbase) // '.' // fname
        end if

     endif

! here he have the raw data as an x by y by number of components array.
! we want to store it as double precision raw data, one component per file.
     do n = 1, nc

! if the component name has a '/', remove it
        tempname = trim(pf%names(comps(n)))
        indslsh = index(tempname, '/', back = .TRUE.)
        if (indslsh > 0) then
           compname = tempname(1:indslsh-1) // tempname(indslsh+1:)
        else
           compname = tempname
        endif

        open (unit=10, file = trim(ofname) // '.' // trim(compname) // '.raw', &
             form="unformatted")

        ! store the array size (nx, ny), the time, and the data itself
        write (unit=10) pf%plo, pf%phi
        write (unit=10) rhi(1)-rlo(1)+1, rhi(2)-rlo(2)+1, rhi(3)-rlo(3)+1
        write (unit=10) pf%tm

        if (n == 1) print *, '  time = ', pf%tm
               
        write (unit=10) rData(:,:,:,n)

        close (10)
        print *, '  wrote: ', trim(ofname) // '.' // trim(compname) // '.raw'
     end do
     call destroy(pf)
  end do

end program fIDLdump3d
