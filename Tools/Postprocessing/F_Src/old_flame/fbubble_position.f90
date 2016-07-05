! compute the extrema of the bubble by computing the
! horizontal and vertical averages of the ash (Mg24)
! mass fractions, and looking where it exceeds the
! threshold at the top, bottom, left, and right extents.

program fdata
  implicit none

  call process_snapshots

end program fdata

subroutine process_snapshots
  use plotfile_module
  use filler_module
!  use utility_module
  use bl_parmparse_module
  use bl_IO_module

  implicit none

  integer :: plo(2), phi(2)
  integer f, i, j, farg
  integer :: n
  integer :: index_store

  integer, parameter :: nc = 2
  integer :: index_comp(nc)
  character (len=64) :: compname(nc)
  logical :: found

  integer :: jmin, jmax, imin, imax, jmass
  real(kind=dp_t) :: dx, dy, cthres

  real(kind=dp_t), allocatable :: c_fab(:,:,:), c_h_avg(:), c_v_avg(:)
  type(plotfile) :: pf
  integer narg
  character(len=128) :: phome
  character(len=256) ::  fname

  real(kind=dp_t) :: yy, ymass, mass, ycm
  
  real(kind=dp_t) :: sum0
  integer unit
  integer ii, jj
  integer nn
  integer :: nx

  real(kind=dp_t) :: dx_cut
  integer, dimension(2) :: tlo(2), thi(2)

  logical :: verbose
  integer :: ntime
  integer :: max_level

  integer :: cut


  cthres = 0.001d0

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

  print *, '#  <x> threshold = ', cthres

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
  compname(1) = "Y(Mg24)"
  compname(2) = "density"

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

  allocate(c_h_avg(plo(2):phi(2)))
  allocate(c_v_avg(plo(1):phi(1)))

  ! loop over all the files, build a single FAB, an then loop over the 
  ! individual components we wish to store, and store them one by one
  dx = (pf%phi(1) - pf%plo(1))/(phi(1)-plo(1)+1)
  dy = (pf%phi(2) - pf%plo(2))/(phi(2)-plo(2)+1)


  nx = phi(1) - plo(1) + 1
  dx_cut = nx/cut
 
  print *, '# time     ymin    ymax    xmin    xmax    Y CM    mass (g/cm)'
  
  do f = 1, ntime

     c_h_avg(:) = 0.0d0
     c_v_avg(:) = 0.0d0

     ymass = 0.d0
     mass = 0.d0
   
     call get_command_argument(f + farg - 1, value = fname)
     if ( f /= 1 ) call build(pf, fname, unit)

     ! we are going to blow this out to a fab in pieces, so we can
     ! handle large domains.  We cut in the x direction.
     do n = 1, cut

        ! allocate the storage for the current subdomain
        tlo(1) = (n-1)*dx_cut
        thi(1) = n*dx_cut -1

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
           c_h_avg(j) = c_h_avg(j) + sum(c_fab(:,j,1))
        enddo

        do i = tlo(1), thi(1)
           c_v_avg(i) = c_v_avg(i) + sum(c_fab(i,:,1))
        enddo

        ! to compute the center of mass (in y), we need to compute
        ! the sum of {y rho Mg} and the sum of {rho Mg}
        do j = tlo(2), thi(2)
           yy = dy*(j+0.5) + pf%plo(2)

           do i = tlo(1), thi(1)
              mass = mass + c_fab(i,j,1)*c_fab(i,j,2)
              ymass = ymass + yy*c_fab(i,j,1)*c_fab(i,j,2)
           enddo

        enddo

        deallocate (c_fab)

     enddo

     ! normalize
     c_h_avg = c_h_avg/(phi(1) - plo(1) + 1)
     c_v_avg = c_v_avg/(phi(2) - plo(2) + 1)

     ycm = ymass/mass
     mass = mass*dx*dy

! find the upper and lower extrema

! lower
     do j = plo(2), phi(2)
        if (c_h_avg(j) > cthres) then
           jmin = j
           exit
        endif
     enddo

! upper
     do j = phi(2), plo(2), -1
        if (c_h_avg(j) > cthres) then
           jmax = j
           exit
        endif
     enddo

     
! find the left and right extrema

! left
     do i = plo(1), phi(1)
        if (c_v_avg(i) > cthres) then
           imin = i
           exit
        endif
     enddo

! right
     do i = phi(1), plo(1), -1
        if (c_v_avg(i) > cthres) then
           imax = i
           exit
        endif
     enddo



! plo begins at index 0, so the first zone center has a coordinate of
! dx*(plo(1)+0.5) + xmin

100 format(1x,7(g14.8,2x))
     print 100, real(pf%tm), &
          (dy*(jmin+0.5) + pf%plo(2)), &
          (dy*(jmax+0.5) + pf%plo(2)), &
          (dx*(imin+0.5) + pf%plo(1)), &
          (dx*(imax+0.5) + pf%plo(1)), &
          (ycm), &
          (mass)


     call destroy(pf)
  end do

end subroutine process_snapshots
