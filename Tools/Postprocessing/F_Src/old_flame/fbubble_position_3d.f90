! compute the extrema of the bubble by computing the
! horizontal and vertical averages of the ash (Mg24)
! mass fractions, and looking where it exceeds the
! threshold at the top, bottom, left, and right extents.
!
! This is the three-dimensional version of the bubble position 
! routine

program fdata
  implicit none

  call process_snapshots

end program fdata

subroutine process_snapshots
  use plotfile_module
  use filler_module
  use bl_parmparse_module
  use bl_IO_module

  implicit none

  integer :: plo(3), phi(3)
  integer f, i, j, k, farg
  integer :: n
  integer :: index_store

  integer, parameter :: nc = 2
  integer :: index_comp(nc)
  character (len=64) :: compname(nc)
  logical :: found

  integer :: kmin, kmax, jmin, jmax, imin, imax, kmass
  real(kind=dp_t) :: dx, dy, dz, cthres

  real(kind=dp_t), allocatable :: c_fab(:,:,:,:)

  ! c_N_avg holds the data averaged over all dimensions BUT N
  real(kind=dp_t), allocatable :: c_x_avg(:), c_y_avg(:), c_z_avg(:)

  type(plotfile) :: pf
  integer narg
  character(len=128) :: phome
  character(len=256) ::  fname

  real(kind=dp_t) :: zz, zmass, mass, zcm
  
  real(kind=dp_t) :: sum0
  integer unit

  integer :: nx, ny, nz

  real(kind=dp_t) :: dx_cut
  integer :: tlo(3), thi(3)

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

  print *, '#  <X> threshold = ', cthres

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
  ! Note, that we don't always save Y(Mg24) in 3-d runs, so we will need to
  ! derive it from the Mg partial density and the total mass density.
  compname(1) = "rho.Y(Mg24)"
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

!  print *, 'plo = ', plo
!  print *, 'phi = ', phi


  allocate(c_x_avg(plo(1):phi(1)))
  allocate(c_y_avg(plo(2):phi(2)))
  allocate(c_z_avg(plo(3):phi(3)))


  ! loop over all the files, build a single FAB, an then loop over the 
  ! individual components we wish to store, and store them one by one
  dx = (pf%phi(1) - pf%plo(1))/(phi(1)-plo(1)+1)
  dy = (pf%phi(2) - pf%plo(2))/(phi(2)-plo(2)+1)
  dz = (pf%phi(3) - pf%plo(3))/(phi(3)-plo(3)+1)


  nx = phi(1) - plo(1) + 1
  ny = phi(2) - plo(2) + 1
  nz = phi(3) - plo(3) + 1

  dx_cut = nx/cut
 
  print *, '# time     xmin    xmax    ymin    ymax    zmin    zmax    Z CM    mass (g)'
  
  do f = 1, ntime

     c_x_avg(:) = 0.0d0
     c_y_avg(:) = 0.0d0
     c_z_avg(:) = 0.0d0

     zmass = 0.d0
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

        tlo(3) = plo(3)
        thi(3) = phi(3)
        
        allocate(c_fab(tlo(1):thi(1), tlo(2):thi(2), tlo(3):thi(3), nc))
        call blow_out_to_sub_fab_3d(c_fab, tlo, thi, pf, index_comp, max_level)

        ! do the lateral averaging -- don't normalize here -- we'll save
        ! that to the end
        
        do k = tlo(3), thi(3)
           c_z_avg(k) = c_z_avg(k) + sum(c_fab(:,:,k,1)/c_fab(:,:,k,2))
        enddo

        do j = tlo(2), thi(2)
           c_y_avg(j) = c_y_avg(j) + sum(c_fab(:,j,:,1)/c_fab(:,j,:,2))
        enddo

        do i = tlo(1), thi(1)
           c_x_avg(i) = c_x_avg(i) + sum(c_fab(i,:,:,1)/c_fab(i,:,:,2))
        enddo

        ! to compute the center of mass (in z), we need to compute
        ! the sum of {z rho Mg} and the sum of {rho Mg}
        do k = tlo(3), thi(3)
           zz = dz*(k+0.5) + pf%plo(3)

           do j = tlo(2), thi(2)
              do i = tlo(1), thi(1)
                 mass = mass + c_fab(i,j,k,1)
                 zmass = zmass + zz*c_fab(i,j,k,1)
              enddo
           enddo

        enddo

        deallocate (c_fab)

     enddo

     ! normalize
     c_x_avg = c_x_avg/(ny*nz)
     c_y_avg = c_y_avg/(nx*nz)
     c_z_avg = c_z_avg/(nx*ny)


     zcm = zmass/mass
     mass = mass*dx*dy*dz


! find the bubble height (zmax and zmin)

! lower
     do k = plo(3), phi(3)
        if (c_z_avg(k) > cthres) then
           kmin = k
           exit
        endif
     enddo

! upper
     do k = phi(3), plo(3), -1
        if (c_z_avg(k) > cthres) then
           kmax = k
           exit
        endif
     enddo

     
! find the lateral extrema

! xmin
     do i = plo(1), phi(1)
        if (c_x_avg(i) > cthres) then
           imin = i
           exit
        endif
     enddo

! xmax
     do i = phi(1), plo(1), -1
        if (c_x_avg(i) > cthres) then
           imax = i
           exit
        endif
     enddo


! ymin
     do j = plo(2), phi(2)
        if (c_y_avg(j) > cthres) then
           jmin = j
           exit
        endif
     enddo

! ymax
     do j = phi(2), plo(2), -1
        if (c_y_avg(j) > cthres) then
           jmax = j
           exit
        endif
     enddo


! plo begins at index 0, so the first zone center has a coordinate of
! dx*(plo(1)+0.5) + xmin

100 format(1x,9(g14.8,2x))
     print 100, real(pf%tm), &
          (dx*(imin+0.5) + pf%plo(1)), &
          (dx*(imax+0.5) + pf%plo(1)), &
          (dy*(jmin+0.5) + pf%plo(2)), &
          (dy*(jmax+0.5) + pf%plo(2)), &
          (dz*(kmin+0.5) + pf%plo(3)), &
          (dz*(kmax+0.5) + pf%plo(3)), &
          (zcm), &
          (mass)


     call destroy(pf)
  end do

end subroutine process_snapshots
