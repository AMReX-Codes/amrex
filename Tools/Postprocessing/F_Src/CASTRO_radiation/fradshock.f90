! extract a 1-d slice of the data (all variables or a single variable)
! along the specified coordinate direction from a plotfile.  The
! plotfile can be 1-, 2-, or 3-d.
!
! This routine is a generalized version is based on fextract3d, but geared
! toward the CASTRO radiating shock problem
!
! We read in all the variables, but only output a subset

program fradshock

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use sort_d_module
  use fundamental_constants_module, only: a_rad

  implicit none

  type(plotfile) pf
  integer :: unit
  integer :: i, j, ii, jj, kk
  integer :: rr, r1
  integer :: uno
  integer :: cnt, max_points, nvs

  real(kind=dp_t) :: dx(MAX_SPACEDIM)
  real(kind=dp_t), pointer :: p(:,:,:,:)
  real(kind=dp_t), allocatable :: sv(:,:)
  integer, allocatable :: isv(:)
  logical, allocatable :: imask(:)
  integer :: nspec, lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)
  integer :: dim

  character(len=256) :: slicefile
  character(len=256) :: pltfile
  integer :: indslsh
  integer :: iloc, jloc, kloc
  integer :: idir
  integer :: narg, farg
  character(len=256) :: fname

  character(len=1) :: dirstr
  integer :: dens_comp, xvel_comp, yvel_comp, zvel_comp, &
       pres_comp, eint_comp, rad_comp, temp_comp

  real(kind=dp_t) :: xmin, xmax, ymin, ymax, zmin, zmax

  unit = unit_new()
  uno =  unit_new()

  slicefile = ''
  pltfile  = ''
  idir = 1

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

     case ('-d', '--direction')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) idir

     case default
        exit

     end select
     farg = farg + 1
  end do

  if ( len_trim(pltfile) == 0 ) then
     print *, " "
     print *, "extract at 1D slice through a plotfile in any coordinate direction"
     print *, "usage: fradshock -p plotfile [-s outfile] [-d dir]"
     print *, "args [-p|--pltfile]   plotfile   : plot file directory (required)"
     print *, "     [-s|--slicefile] slice file : slice file          (optional)"
     print *, "     [-d|--direction] idir       : slice direction {1 (default), 2, or 3}"
     print *, " "
     stop
  end if

  ! slicefile not defined, default to plotfile.slice
  if ( len_trim(slicefile) == 0 ) then

     ! get basename of file
     indslsh = index(pltfile, '/', back = .TRUE.)

     if ( indslsh /= 0 ) then
        slicefile = trim(pltfile(:indslsh-1)) // ".slice"
     else
        slicefile = trim(pltfile) // ".slice"
     end if
     
  endif

  print *, 'pltfile   = "', trim(pltfile), '"'
  print *, 'slicefile = "', trim(slicefile), '"'
  print *, 'direction =  ', idir

  call build(pf, pltfile, unit)

  nvs = pf%nvars

  ! get the indices for the variables we need
  dens_comp = plotfile_var_index(pf, "density")
  xvel_comp = plotfile_var_index(pf, "x_velocity")
  if (dim > 1) yvel_comp = plotfile_var_index(pf, "y_velocity")
  if (dim > 2) zvel_comp = plotfile_var_index(pf, "z_velocity")
  eint_comp = plotfile_var_index(pf, "eint_E")  
  temp_comp = plotfile_var_index(pf, "Temp")  
  rad_comp  = plotfile_var_index(pf, "rad")  
  pres_comp = plotfile_var_index(pf, "pressure")  

  if (dens_comp < 0 .or. xvel_comp < 0 .or. (yvel_comp < 0 .and. dim > 1) .or. &
       (zvel_comp < 0 .and. dim > 2) .or. eint_comp < 0 .or. temp_comp < 0 .or. &
       rad_comp < 0 .or. pres_comp < 0) then
     call bl_error("ERROR: variable components undefined")
  endif


  do i = 1, pf%flevel
     call fab_bind_level(pf, i)
  end do

  ! get the index bounds and dx for the coarse level.  Note, lo and hi are
  ! ZERO based indicies
  lo(:) = ZERO
  hi(:) = ZERO

  lo = lwb(plotfile_get_pd_box(pf, 1))
  hi = upb(plotfile_get_pd_box(pf, 1))

  dx = plotfile_get_dx(pf, 1)

  dim = pf%dim

  ! get the physical domain extrema
  xmin = pf%plo(1)
  xmax = pf%phi(1)

  if (dim >= 2) then
     ymin = pf%plo(2)
     ymax = pf%phi(2)
  endif
  
  if (dim == 3) then
     zmin = pf%plo(3)
     zmax = pf%phi(3)
  endif


  ! get the index bounds for the finest level
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

  ! compute the index of the center of the domain on the coarse grid.
  ! These are used to set the position of the slice in the transverse
  ! direction.
  iloc = (hi(1)-lo(1)+1)/2 + lo(1)
  if (dim > 1) jloc = (hi(2)-lo(2)+1)/2 + lo(2)
  if (dim > 2) kloc = (hi(3)-lo(3)+1)/2 + lo(3)

  cnt = 0

  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here

  if (idir > dim) call bl_error("idir > dim")

  if (idir == 1) then
     allocate(imask(flo(1):fhi(1)))
     max_points = fhi(1) - flo(1) + 1

  else if (idir == 2 .and. dim >= 2) then
     allocate(imask(flo(2):fhi(2)))
     max_points = fhi(2) - flo(2) + 1

  else if (idir == 3 .and. dim == 3) then
     allocate(imask(flo(3):fhi(3)))
     max_points = fhi(3) - flo(3) + 1

  else 
     call bl_error("invalid direction")
  endif

  imask = .true.

  ! allocate storage for the data -- we will assume that we are
  ! completely refined, and allocate enough storage for that,
  ! although, we won't really need all of this
  !
  ! note: sv(:,1) will be the coordinate information.  
  ! the variables will be stored in sv(:,2:nvvs+1)
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
        lo(:) = 1
        hi(:) = 1
        lo = lwb(get_box(pf, i, j))
        hi = upb(get_box(pf, i, j))

        select case (idir)

        case (1)

           if (dim == 1) then

              p => dataptr(pf, i, j)
              jj = 1
              kk = 1

              do ii = lo(1), hi(1)
                 if ( any(imask(ii*r1:(ii+1)*r1-1) ) ) then
                    cnt = cnt + 1
                    
                    sv(cnt,1) = xmin + (ii + HALF)*dx(1)/rr
                    sv(cnt,2:) = p(ii,jj,kk,:)
                    
                    imask(ii*r1:(ii+1)*r1-1) = .false.
                 end if
              end do

           else

              ! if the current patch stradles our slice, then get a data
              ! pointer to it
              if ( rr*jloc >= lo(2) .and. rr*jloc <= hi(2) .and. &
                   ( (dim .eq. 2) .or. (rr*kloc >= lo(3) .and. rr*kloc <= hi(3)) ) ) then
                 p => dataptr(pf, i, j)
                 jj = jloc*rr
                 if (dim .eq. 3) then
                    kk = kloc*rr
                 else
                    kk = 1
                 end if

                 ! loop over all of the zones in the slice direction.
                 ! Here, we convert the cell-centered indices at the
                 ! current level into the corresponding RANGE on the
                 ! finest level, and test if we've stored data in any of
                 ! those locations.  If we haven't then we store this
                 ! level's data and mark that range as filled.
                 do ii = lo(1), hi(1)
                    if ( any(imask(ii*r1:(ii+1)*r1-1) ) ) then
                       cnt = cnt + 1

                       sv(cnt,1) = xmin + (ii + HALF)*dx(1)/rr
                       sv(cnt,2:) = p(ii,jj,kk,:)

                       imask(ii*r1:(ii+1)*r1-1) = .false.
                    end if
                 end do

              end if

           end if

        case (2)

           ! if the current patch stradles our slice, then get a data
           ! pointer to it
           if ( rr*iloc >= lo(1) .and. rr*iloc <= hi(1) .and. &
                ( (dim .eq. 2) .or. (rr*kloc >= lo(3) .and. rr*kloc <= hi(3)) ) ) then
              p => dataptr(pf, i, j)
              ii = iloc*rr
              if (dim .eq. 3) then
                kk = kloc*rr
              else
                kk = 1
              end if

           
              ! loop over all of the zones in the slice direction.
              ! Here, we convert the cell-centered indices at the
              ! current level into the corresponding RANGE on the
              ! finest level, and test if we've stored data in any of
              ! those locations.  If we haven't then we store this
              ! level's data and mark that range as filled.
              do jj = lo(2), hi(2)
                 if ( any(imask(jj*r1:(jj+1)*r1-1) ) ) then
                    cnt = cnt + 1

                    sv(cnt,1) = ymin + (jj + HALF)*dx(2)/rr
                    sv(cnt,2:) = p(ii,jj,kk,:)

                    imask(jj*r1:(jj+1)*r1-1) = .false.
                 end if
              end do

           end if

        case (3)

           ! if the current patch stradles our slice, then get a data
           ! pointer to it
           if ( rr*iloc >= lo(1) .and. rr*iloc <= hi(1) .and. &
                rr*jloc >= lo(2) .and. rr*jloc <= hi(2)) then
              p => dataptr(pf, i, j)
              ii = iloc*rr
              jj = jloc*rr

           
              ! loop over all of the zones in the slice direction.
              ! Here, we convert the cell-centered indices at the
              ! current level into the corresponding RANGE on the
              ! finest level, and test if we've stored data in any of
              ! those locations.  If we haven't then we store this
              ! level's data and mark that range as filled.
              do kk = lo(3), hi(3)
                 if ( any(imask(kk*r1:(kk+1)*r1-1) ) ) then
                    cnt = cnt + 1

                    sv(cnt,1) = zmin + (kk + HALF)*dx(3)/rr
                    sv(cnt,2:) = p(ii,jj,kk,:)

                    imask(kk*r1:(kk+1)*r1-1) = .false.
                 end if
              end do

           end if

        end select

        
     end do

     ! adjust r1 for the next lowest level
     if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
  end do


  ! sort the data based on the coordinates
  call sort(sv(1:cnt,1),isv(1:cnt))

  if (idir == 1) then
     dirstr = "x"
  else if (idir == 2) then
     dirstr = "y"
  else 
     dirstr = "z"
  endif


 998 format("# 1-d slice in ", a1 "-direction, file: ", a)
 999 format("# time = ", g24.12)
1000 format("#",100(a24,1x))
1001 format(1x, 100(g24.12,1x))

  ! slicefile
  open(unit=uno, file=slicefile, status = 'replace')
  write(uno,998) dirstr, trim(pltfile)
  write(uno,999) pf%tm

  ! output selected variables

  ! write the header
  if (dim == 1) then
     write(uno,1000) dirstr, "density", &
          "x-velocity", &
          "pressure", "int. energy", "temperature", &
          "rad energy", "rad temp"
  else if (dim == 2) then
     write(uno,1000) dirstr, "density", &
          "x-velocity", "y-velocity", &
          "pressure", "int. energy", "temperature", &
          "rad energy", "rad temp"
  else if (dim == 3) then
     write(uno,1000) dirstr, "density", &
          "x-velocity", "y-velocity", "z-velocity", &
          "pressure", "int. energy", "temperature", &
          "rad energy", "rad temp"
  endif

  ! Use this to protect against a number being xx.e-100 
  !   which will print without the "e"
  do i=1,cnt
     do j=2,nvs+1
        if (abs(sv(isv(i),j)) .lt. 1.d-99) then
           sv(isv(i),j) = 0.d0
        end if
     end do
  end do

  ! write the data in columns
  do i = 1, cnt
     if (dim == 1) then
        write(uno,1001) sv(isv(i),1), sv(isv(i),dens_comp+1), &
             sv(isv(i),xvel_comp+1), &
             sv(isv(i),pres_comp+1), sv(isv(i),eint_comp+1), sv(isv(i),temp_comp+1), &
             sv(isv(i),rad_comp+1), (sv(isv(i),rad_comp+1)/a_rad)**0.25
     else if (dim == 2) then
        write(uno,1001) sv(isv(i),1), sv(isv(i),dens_comp+1), &
             sv(isv(i),xvel_comp+1), sv(isv(i),yvel_comp+1), &  
             sv(isv(i),pres_comp+1), sv(isv(i),eint_comp+1), sv(isv(i),temp_comp+1), &
             sv(isv(i),rad_comp+1), (sv(isv(i),rad_comp+1)/a_rad)**0.25
     else if (dim == 3) then
        write(uno,1001) sv(isv(i),1), sv(isv(i),dens_comp+1), &
             sv(isv(i),xvel_comp+1), sv(isv(i),yvel_comp+1), sv(isv(i),zvel_comp+1), &
             sv(isv(i),pres_comp+1), sv(isv(i),eint_comp+1), sv(isv(i),temp_comp+1), &
             sv(isv(i),rad_comp+1), (sv(isv(i),rad_comp+1)/a_rad)**0.25
     endif
  end do
  
  close(unit=uno)

  do i = 1, pf%flevel
     call fab_unbind_level(pf, i)
  end do

  call destroy(pf)

end program fradshock
