! Radially average a 3-d sedov problem to produce rho, u, and p as a
! function of r, for comparison to the analytic solution.

program fextract3d

  use f2kcli
  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module

  implicit none

  type(plotfile) pf
  integer :: unit
  integer :: i, j, ii, jj, kk
  real(kind=dp_t) :: xx, yy, zz
  integer :: rr, r1
  integer :: uno

  integer :: nbins
  real(kind=dp_t), allocatable :: r(:)
  real(kind=dp_t) :: maxdist, x_maxdist, y_maxdist, z_maxdist
  real(kind=dp_t) :: xctr, yctr, zctr

  real(kind=dp_t) :: dx(MAX_SPACEDIM)
  real(kind=dp_t) :: dx_fine

  real(kind=dp_t) :: r_zone
  integer :: index

  real(kind=dp_t), pointer :: p(:,:,:,:)

  integer, allocatable :: ncount(:)
  real(kind=dp_t), allocatable :: dens_bin(:), vel_bin(:), pres_bin(:), e_bin(:)

  integer :: dens_comp, xmom_comp, ymom_comp, zmom_comp, pres_comp, rhoe_comp

  logical, allocatable :: imask(:,:,:)
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
  zctr = HALF


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

     case ('--zctr')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) zctr

     case default
        exit

     end select
     farg = farg + 1
  end do

  if ( len_trim(pltfile) == 0 .OR. len_trim(slicefile) == 0 ) then
     print *, "usage: executable_name args"
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
  dens_comp = plotfile_var_index(pf, "density")

  ! x-momentum
  xmom_comp = plotfile_var_index(pf, "xmom") 

  ! y-momentum
  ymom_comp = plotfile_var_index(pf, "ymom") 

  ! z-momentum
  zmom_comp = plotfile_var_index(pf, "zmom")

  ! pressure
  pres_comp = plotfile_var_index(pf, "pressure") 

  ! rho*e
  rhoe_comp = plotfile_var_index(pf, "rho_e") 

  if ( dens_comp < 0 .or. xmom_comp < 0 .or. ymom_comp < 0 .or. &
       zmom_comp < 0 .or. pres_comp < 0 .or. rhoe_comp < 0) then
     call bl_error("ERROR: variable(s) not defined")
  endif


  ! get the index bounds and dx for the coarse level.  Note, lo and hi are
  ! ZERO based indicies
  lo = lwb(plotfile_get_pd_box(pf, 1))
  hi = upb(plotfile_get_pd_box(pf, 1))

  dx = plotfile_get_dx(pf, 1)


  ! get the index bounds for the finest level
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

  print *, flo
  print *, fhi

  ! compute the size of the radially-binned array -- we'll do it to
  ! the furtherest corner of the domain
  x_maxdist = max(abs(pf%phi(1) - xctr), abs(pf%plo(1) - xctr))
  y_maxdist = max(abs(pf%phi(2) - yctr), abs(pf%plo(2) - yctr))
  z_maxdist = max(abs(pf%phi(3) - zctr), abs(pf%plo(3) - zctr))
  
  maxdist = sqrt(x_maxdist**2 + y_maxdist**2 + z_maxdist**2)

  dx_fine = minval(plotfile_get_dx(pf, pf%flevel))
  nbins = int(maxdist/dx_fine)

  allocate(r(0:nbins-1))

  do i = 0, nbins-1
     r(i) = (dble(i) + HALF)*dx_fine
  enddo


  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here
  allocate(imask(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3)))
  imask(:,:,:) = .true.

  ! allocate storage for the data 
  allocate(dens_bin(0:nbins-1))
  allocate( vel_bin(0:nbins-1))
  allocate(pres_bin(0:nbins-1))
  allocate(   e_bin(0:nbins-1))
  allocate(  ncount(0:nbins-1))

  ncount(:) = 0
  dens_bin(:) = ZERO
  vel_bin(:) = ZERO
  pres_bin(:) = ZERO
     e_bin(:) = ZERO

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
        lo = lwb(get_box(pf, i, j))
        hi = upb(get_box(pf, i, j))

        ! get a pointer to the current patch
        p => dataptr(pf, i, j)

        
        ! loop over all of the zones in the patch.  Here, we convert
        ! the cell-centered indices at the current level into the
        ! corresponding RANGE on the finest level, and test if we've
        ! stored data in any of those locations.  If we haven't then
        ! we store this level's data and mark that range as filled.
        do kk = lbound(p,dim=3), ubound(p,dim=3)
           zz = (kk + HALF)*dx(3)/rr

           do jj = lbound(p,dim=2), ubound(p,dim=2)
              yy = (jj + HALF)*dx(2)/rr

              do ii = lbound(p,dim=1), ubound(p,dim=1)
                 xx = (ii + HALF)*dx(1)/rr

                 if ( any(imask(ii*r1:(ii+1)*r1-1, &
                                jj*r1:(jj+1)*r1-1, &
                                kk*r1:(kk+1)*r1-1) ) ) then

                    r_zone = sqrt((xx-xctr)**2 + (yy-yctr)**2 + (zz-zctr)**2)

                    index = r_zone/dx_fine

                    ! weight the zone's data by its size
                    dens_bin(index) = dens_bin(index) + &
                         p(ii,jj,kk,dens_comp)*r1**3

                    vel_bin(index) = vel_bin(index) + &
                         (sqrt(p(ii,jj,kk,xmom_comp)**2 + &
                               p(ii,jj,kk,ymom_comp)**2 + &
                               p(ii,jj,kk,zmom_comp)**2)/ &
                               p(ii,jj,kk,dens_comp))*r1**3

                    pres_bin(index) = pres_bin(index) + &
                         p(ii,jj,kk,pres_comp)*r1**3

                    e_bin(index) = e_bin(index) + &
                         (p(ii,jj,kk,rhoe_comp)/p(ii,jj,kk,dens_comp))*r1**3

                    ncount(index) = ncount(index) + r1**3

                    imask(ii*r1:(ii+1)*r1-1, &
                          jj*r1:(jj+1)*r1-1, &
                          kk*r1:(kk+1)*r1-1) = .false.
                 end if

              end do
           enddo
        enddo

     end do

     ! adjust r1 for the next lowest level
     if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
  end do


  ! normalize
  do i = 0, nbins-1
     if (ncount(i) /= 0) then
        dens_bin(i) = dens_bin(i)/ncount(i)
        vel_bin(i) = vel_bin(i)/ncount(i)
        pres_bin(i) = pres_bin(i)/ncount(i)
           e_bin(i) =    e_bin(i)/ncount(i)
     endif
  enddo

1000 format("#",100(a24,1x))
1001 format(1x, 100(g24.13,1x))

  ! slicefile
  open(unit=uno, file=slicefile, status = 'replace')

  ! write the header
  write(uno,1000) "x", "density", "velocity", "pressure", "int. energy"

  ! write the data in columns
  do i = 0, nbins-1
     ! Use this to protect against a number being xx.e-100
     !   which will print without the "e"
     if (abs(dens_bin(i)) .lt. 1.d-99)  dens_bin(i) = 0.d0
     if (abs( vel_bin(i)) .lt. 1.d-99)   vel_bin(i) = 0.d0
     if (abs(pres_bin(i)) .lt. 1.d-99)  pres_bin(i) = 0.d0
     if (abs(   e_bin(i)) .lt. 1.d-99)     e_bin(i) = 0.d0
     write(uno,1001) r(i), dens_bin(i), vel_bin(i), pres_bin(i), e_bin(i)
  end do

  close(unit=uno)

  do i = 1, pf%flevel
     call fab_unbind_level(pf, i)
  end do

  call destroy(pf)

end program fextract3d
