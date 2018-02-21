! extract a 1-d slice of the data (all variables or a single variable)
! along the specified coordinate direction from a plotfile.  The
! plotfile can be 1-, 2-, or 3-d.
!
! Note: by default, this routine slices through the center of the
! domain in the transverse directions.  At the moment, there is no
! runtime method of overridding this.

program fextract3d

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use sort_d_module

  implicit none

  type(plotfile) pf
  integer :: unit, uno, un_ji
  integer :: i, j, ii, jj, kk
  integer :: rr, r1
  integer :: cnt, max_points, nvs

  real(kind=dp_t) :: dx(MAX_SPACEDIM)
  real(kind=dp_t), pointer :: p(:,:,:,:)
  real(kind=dp_t), allocatable :: sv(:,:)
  integer, allocatable :: isv(:)
  logical, allocatable :: imask(:)
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)
  integer :: dim
  integer :: fine_level, coarse_level

  character(len=256) :: slicefile
  character(len=256) :: pltfile
  integer :: indslsh
  integer :: iloc, jloc, kloc
  integer :: idir
  integer :: narg, farg
  character(len=256) :: fname, varnames, line, temp

  character(len=1) :: dirstr
  integer :: ivar, idx
  integer, allocatable :: var_indices(:)

  real(kind=dp_t) :: xmin, xmax, ymin, ymax, zmin, zmax

  logical :: ji_exist, valid, center
  integer :: io

  unit = unit_new()

  slicefile = ''
  pltfile  = ''
  idir = 1
  varnames = ''
  center = .true.
  coarse_level = 1
  fine_level = -1 ! This will be fixed later if the user doesn't select one

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

     case ('-v', '--variable')
        farg = farg + 1
        call get_command_argument(farg, value = varnames)

     case ('-l', '--lower_left')
        center = .false.

     case ('-c', '--coarse_level')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) coarse_level

     case ('-f', '--fine_level')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) fine_level

     case default
        exit

     end select
     farg = farg + 1
  end do

  if (len_trim(pltfile) == 0 .and. farg <= narg) then
     call get_command_argument(farg, value = pltfile)
  endif

  if ( len_trim(pltfile) == 0 ) then
     print *, " "
     print *, "Extract at 1D slice through a plotfile in any coordinate direction."
     print *, "Works with 1-, 2-, or 3-d datasets."
     print *, " "
     print *, "Usage:"
     print *, "   fextract [-p plotfile] [-s outfile] [-d dir] [-v variable] [-c coarse_level] [-f fine_level] plotfile"
     print *, " "
     print *, "args [-p|--pltfile]   plotfile        : plot file directory (depreciated, optional)"
     print *, "     [-s|--slicefile] slice file      : slice file          (optional)"
     print *, "     [-d|--direction] idir            : slice direction {1 (default), 2, or 3}"
     print *, "     [-v|--variable]  varname(s)      : only output the values of variable varname"
     print *, "                                        (space separated string for multiple variables)"
     print *, "     [-l|--lower_left]                : slice through lower left corner instead of center"
     print *, "     [-c|--coarse_level] coarse level : coarsest level to extract from"
     print *, "     [-f|--fine_level]   fine level   : finest level to extract from"
     print *, " "
     print *, "Note the plotfile at the end of the commandline is only required if you do"
     print *, "not use the deprecated '-p' option"
     print *, " "
     print *, "If a job_info file is present in the plotfile, that information is made"
     print *, "available at the end of the slice file (commented out), for reference."
     print *, " "
     stop
  end if

  ! slicefile not defined, default to plotfile.slice
  if ( len_trim(slicefile) == 0 ) then

     ! get basename of file
     indslsh = index(pltfile, '/', back = .TRUE.)

     if ( indslsh == len(trim(pltfile)) ) then
        slicefile = trim(pltfile(:indslsh-1)) // ".slice"
     else
        slicefile = trim(pltfile) // ".slice"
     end if

  endif

  select case (idir)
  case (1)
     print *, 'slicing along x-direction and outputting to ', trim(slicefile)
  case (2)
     print *, 'slicing along y-direction and outputting to ', trim(slicefile)
  case (3)
     print *, 'slicing along z-direction and outputting to ', trim(slicefile)
  end select


  call build(pf, pltfile, unit)

  nvs = pf%nvars

  ! if we are outputting only a single variable, make sure it exists
  ivar = -1
  if (varnames /= '') then
     ! flag that indicates we have variable
     ivar = 0
     allocate(var_indices(pf%nvars))
     var_indices(:) = 0

     do while (.not. trim(varnames) == "")
        ivar = ivar + 1

        idx = index(varnames, " ")
        temp = varnames(:idx)
        varnames = trim(adjustl(varnames(idx+1:)))

        var_indices(ivar) = plotfile_var_index(pf, trim(temp))
     enddo
  endif

  ! Get the index bounds and dx for the coarse level.  
  ! Note, lo and hi are 0-based indices
  lo(:) = 0
  hi(:) = 0

  lo(1:pf%dim) = lwb(plotfile_get_pd_box(pf, 1))
  hi(1:pf%dim) = upb(plotfile_get_pd_box(pf, 1))

  dx(1:pf%dim) = plotfile_get_dx(pf, 1)

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
  flo(1:pf%dim) = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi(1:pf%dim) = upb(plotfile_get_pd_box(pf, pf%flevel))

  ! compute the index of the center or lower left of the domain on the
  ! coarse grid.  These are used to set the position of the slice in
  ! the transverse direction.

  if (center) then
     iloc = (hi(1)-lo(1)+1)/2 + lo(1)
     if (dim > 1) jloc = (hi(2)-lo(2)+1)/2 + lo(2)
     if (dim > 2) kloc = (hi(3)-lo(3)+1)/2 + lo(3)
  else
     iloc = 0
     if (dim > 1) jloc = 0
     if (dim > 2) kloc = 0
  endif

  cnt = 0

  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here

  if (idir > dim) call bl_error("idir > dim")

  select case (idir)
  case (1)
     allocate(imask(flo(1):fhi(1)))
     max_points = fhi(1) - flo(1) + 1

  case (2)
     allocate(imask(flo(2):fhi(2)))
     max_points = fhi(2) - flo(2) + 1

  case (3)
     allocate(imask(flo(3):fhi(3)))
     max_points = fhi(3) - flo(3) + 1
  end select

  imask = .true.

  ! allocate storage for the data -- we will assume that we are
  ! completely refined, and allocate enough storage for that,
  ! although, we won't really need all of this
  !
  ! note: sv(:,1) will be the coordinate information.
  ! the variables will be stored in sv(:,2:nvvs+1)
  allocate(sv(max_points,nvs+1), isv(max_points))

  ! If the user didn't select a finest level, assume we want the finest level available

  if (fine_level < 0) then
     fine_level = pf%flevel
  end if

  ! sanity check on valid selected levels

  if (fine_level > pf%flevel .or. coarse_level < 1 .or. coarse_level > fine_level) then
     call bl_error("Invalid level selection")
  end if

  ! loop over the data, starting at the finest grid, and if we haven't
  ! already store data in that grid location (according to imask),
  ! store it.  We'll put it in the correct order later.

  cnt = 0

  ! r1 is the factor between the current level grid spacing and the
  ! FINEST level
  r1  = 1

  do i = fine_level, coarse_level, -1

     ! rr is the factor between the COARSEST level grid spacing and
     ! the current level
     rr = product(pf%refrat(1:i-1,1))

     do j = 1, nboxes(pf, i)

        ! read in the data 1 patch at a time
        call fab_bind(pf, i, j)

        lo(:) = 1
        hi(:) = 1
        lo(1:pf%dim) = lwb(get_box(pf, i, j))
        hi(1:pf%dim) = upb(get_box(pf, i, j))

        p => dataptr(pf, i, j)

        select case (idir)

        case (1)

           if (dim == 1) then

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

        call fab_unbind(pf, i, j)
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


 998 format("# 1-d slice in ", a1, "-direction, file: ", a)
 999 format("# time = ", g24.12)
1000 format("#",100(a24,1x))
1001 format(1x, 100(g24.12,1x))

  ! slicefile
  uno =  unit_new()
  open(unit=uno, file=slicefile, status = 'replace')

  write(uno,998) dirstr, trim(pltfile)
  write(uno,999) pf%tm

  if (ivar == -1) then

     ! output all variables

     ! write the header
     un_ji = unit_new()
     write(uno,1000) dirstr, (trim(pf%names(j)),j=1,nvs)

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
        write(uno,1001) sv(isv(i),1), (sv(isv(i),j),j=2,nvs+1)
     end do

  else

     ! output a specific variables

     ! write the header
     write(uno,1000, advance="no") dirstr
     do i = 1, ivar
        write(uno,fmt="(a24,1x)", advance="no") trim(pf%names(var_indices(i)))
     enddo
     write(uno, *) ""

     ! Use this to protect against a number being xx.e-100
     !   which will print without the "e"
     do j=1,cnt
        do i = 1, ivar
           if (abs(sv(isv(j),var_indices(i)+1)) .lt. 1.d-99) then
              sv(isv(j),var_indices(i)+1) = 0.d0
           end if
        enddo
     enddo

     ! write the data in columns
     do j = 1, cnt
        write(uno,1001, advance="no") sv(isv(j),1)
        do i = 1, ivar
           write(uno, "(g24.12,1x)", advance="no") sv(isv(j), var_indices(i)+1)
        enddo
        write(uno, *) ""
     enddo

  endif

  ! job_info? if so write it out to the slice file end
  inquire (file=trim(pltfile)//"/job_info", exist=ji_exist)
  if (ji_exist) then
     open(unit=un_ji, file=trim(pltfile)//"/job_info", status="old")
     write (uno, fmt="(a)") " "
     valid = .true.
     do while (valid)
        read  (un_ji, fmt="(a)", iostat=io) line
        if (io < 0) then
           valid = .false.
           exit
        endif
        write (uno, fmt="('#', a)") trim(line)
     enddo

  endif

  close(unit=uno)

  call destroy(pf)

end program fextract3d
