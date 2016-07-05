! produce an image of a 2-d slice through a 3-d dataset

program fsnapshot3d

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use ppm_util_module

  implicit none

  integer, parameter :: NCOLOR = 256
  integer r(NCOLOR), g(NCOLOR), b(NCOLOR), a(NCOLOR)

  real(kind=dp_t) :: gmx, gmn
  real(kind=dp_t) :: def_mx, def_mn
  logical :: ldef_mx, ldef_mn, do_log

  type(plotfile) :: pf
  integer :: unit

  integer :: i, j, ii, jj, kk
  integer :: iloc, jloc, kloc

  integer :: ndir, ndir_pass, ndir_start, ndir_end
  character (len=2) :: ndir_string

  integer :: comp
  character (len=64) :: compname

  logical :: origin

  integer :: rr, r1

  real(kind=dp_t) :: dx(MAX_SPACEDIM)
  real(kind=dp_t) :: dx_fine

  real(kind=dp_t), pointer :: p(:,:,:,:)

  logical, allocatable :: imask(:,:)
  real(kind=dp_t), allocatable :: slicedata(:,:)
  integer, allocatable :: intdata(:,:)

  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  character(len=256) :: pltfile
  character(len=256) :: ofname

  integer :: indslsh

  integer :: narg, farg
  character(len=256) :: fname

  character(len=128) :: phome
  character(len=256) :: pfname

  unit = unit_new()

  ! set the defaults
  pltfile  = ''
  ndir_pass = 1
  compname = "density"
  ldef_mx = .false.
  ldef_mn = .false.
  do_log = .false.
  origin = .false.

  ! process the runtime arguments
  narg = command_argument_count()

  call get_environment_variable("HOME", value = phome)
  pfname = trim(phome) // "/.amrvis.Palette"


  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-p', '--pltfile')
        farg = farg + 1
        call get_command_argument(farg, value = pltfile)

     case ('--palette')
              farg = farg + 1
        call get_command_argument(farg, value = pfname)

     case ('-n', '--normaldir')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) ndir_pass

     case ('-cname','--component_name')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        compname = trim(fname)

     case ('-M', '--max')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) def_mx
        ldef_mx = .true.

     case ('-m', '--min')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) def_mn
        ldef_mn = .true.

     case ('-l', '--log')
        do_log = .true.


     case ('--origin')
        origin = .true.

     case default
        exit

     end select
     farg = farg + 1
  end do

  if ( len_trim(pltfile) == 0 .OR. len_trim(compname) == 0 .OR. ndir_pass == -1) then
     print *, "produce an image of a 2-d slice through a 3-d dataset"
     print *, " "
     print *, "usage: fsnapshot3d args"
     print *, "args [-p    |--pltfile]   plotfile   : plot file directory (required)"
     print *, " "
     print *, "     [-n    |--normaldir] {0,1,2,3}  : direction normal to slice (required)"
     print *, "                                       (note 0 means do ALL directions)"
     print *, " "
     print *, "     [-cname|--compname]  name       : variable to plot (default: density)"
     print *, "     --palette pfname                : use the file pfname as the Palette"
     print *, "                                       (default ~/amrvis.Palette)"
     print *, " "
     print *, "     -m val                          : set the minimum value of the data to val"
     print *, "     -M val                          : set the maximum value of the data to val"
     print *, "     [-l|--log]                      : toggle log plot"
     stop
  end if


  ! get the palette
  if ( pfname /= '') then
     call load_palette(pfname, r, g, b, a)
  endif

  ! build the plotfile to get the level and component information
  call build(pf, pltfile, unit)


  ! figure out the variable indices
  comp = plotfile_var_index(pf, compname)
  if (comp < 0) then
     call bl_error("ERROR: variable not found in plotfile")
  endif

  ! get the extrema for scaling the plot
  gmx = plotfile_maxval(pf, comp, pf%flevel)
  gmn = plotfile_minval(pf, comp, pf%flevel)

  print *, "plotfile variable maximum = ", gmx
  print *, "plotfile variable minimum = ", gmn

  if (ldef_mx) then 
     print *, "resetting variable maximum to ", def_mx
     gmx = def_mx
  endif

  if (ldef_mn) then
     print *, "resetting variable minimum to ", def_mn
     gmn = def_mn
  endif

  if (do_log) then
     gmn = log10(gmn)
     gmx = log10(gmx)
  endif


  ! get the index bounds and dx for the coarse level.  Note, lo and hi are
  ! ZERO based indicies
  lo = lwb(plotfile_get_pd_box(pf, 1))
  hi = upb(plotfile_get_pd_box(pf, 1))

  dx = plotfile_get_dx(pf, 1)


  ! get the index bounds for the finest level
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

  dx_fine = minval(plotfile_get_dx(pf, pf%flevel))


  ! set the location of the slices to be the center of the grid
  ! as defined on the coarse grid
  if (origin) then
     iloc = lo(1)
     jloc = lo(2)
     kloc = lo(3)
  else
     iloc = (hi(1)-lo(1)+1)/2 + lo(1)
     jloc = (hi(2)-lo(2)+1)/2 + lo(2)
     kloc = (hi(3)-lo(3)+1)/2 + lo(3)
  endif


  ! produce the image perpendicuar to the desired normal direction
  if (ndir_pass == 0) then
     ndir_start = 1
     ndir_end   = 3
  else
     ndir_start = ndir_pass
     ndir_end   = ndir_pass
  endif


  do ndir = ndir_start, ndir_end


     ! imask will be set to false if we've already output the data.
     ! Note, imask is defined in terms of the finest level.  As we loop
     ! over levels, we will compare to the finest level index space to
     ! determine if we've already output here
     if (ndir == 1) then
        allocate(imask(flo(2):fhi(2),flo(3):fhi(3)))

     else if (ndir == 2) then
        allocate(imask(flo(1):fhi(1),flo(3):fhi(3)))

     else if (ndir == 3) then
        allocate(imask(flo(1):fhi(1),flo(2):fhi(2)))     

     else
        call bl_error("invalid ndir")
     endif

     imask(:,:) = .true.


     ! allocate storage for the data slice plane
     if (ndir == 1) then
        allocate(slicedata(flo(2):fhi(2),flo(3):fhi(3)))
        allocate(  intdata(flo(2):fhi(2),flo(3):fhi(3)))

     else if (ndir == 2) then
        allocate(slicedata(flo(1):fhi(1),flo(3):fhi(3)))
        allocate(  intdata(flo(1):fhi(1),flo(3):fhi(3)))
        
     else if (ndir == 3) then
        allocate(slicedata(flo(1):fhi(1),flo(2):fhi(2)))
        allocate(  intdata(flo(1):fhi(1),flo(2):fhi(2)))     
        
     else
        call bl_error("invalid ndir")
     endif
     
     slicedata(:,:) = ZERO
     intdata(:,:) = ZERO


     !-------------------------------------------------------------------------
     ! loop over the data, starting at the finest grid, and if we haven't
     ! already store data in that grid location (according to imask),
     ! store it.  
     !-------------------------------------------------------------------------


     ! r1 is the factor between the current level grid spacing and the
     ! FINEST level
     r1  = 1

     do i = pf%flevel, 1, -1

        ! rr is the factor between the COARSEST level grid spacing and
        ! the current level
        rr = product(pf%refrat(1:i-1,1))

        do j = 1, nboxes(pf, i)

           call fab_bind_comp_vec(pf, i, j, (/comp/) )

           lo = lwb(get_box(pf, i, j))
           hi = upb(get_box(pf, i, j))


           select case (ndir)

           case (1)

              ! produce a slice perpendicular to the x-axis

              ! if the current patch stradles our slice plane, then get a data
              ! pointer to it
              if (rr*iloc >= lo(1) .and. rr*iloc <= hi(1)) then
              
                 p => dataptr(pf, i, j)
                 lo(:) = lwb(get_box(pf, i, j))
                 hi(:) = upb(get_box(pf, i, j))

                 ii = iloc*rr

                 ! loop over all of the zones in the patch.  Here, we convert
                 ! the cell-centered indices at the current level into the
                 ! corresponding RANGE on the finest level, and test if we've
                 ! stored data in any of those locations.  If we haven't then
                 ! we store this level's data and mark that range as filled.
                 do kk = lo(3), hi(3)
                    do jj = lo(2), hi(2)

                       if ( any(imask(jj*r1:(jj+1)*r1-1, &
                                      kk*r1:(kk+1)*r1-1) ) ) then
                 
                          ! since we've only bound one component to pf, we 
                          ! index p with 1 for the component
                          if (do_log) then
                             slicedata(jj*r1:(jj+1)*r1-1, &
                                       kk*r1:(kk+1)*r1-1) = log10(p(ii,jj,kk,1))
                          else
                             slicedata(jj*r1:(jj+1)*r1-1, &
                                       kk*r1:(kk+1)*r1-1) = p(ii,jj,kk,1)
                          endif

                          imask(jj*r1:(jj+1)*r1-1, &
                                kk*r1:(kk+1)*r1-1) = .false.
                 
                       end if

                    end do
                 enddo

              endif


           case (2)

              ! produce a slice perpendicular to the y-axis

              ! if the current patch stradles our slice plane, then get a data
              ! pointer to it
              if (rr*jloc >= lo(2) .and. rr*jloc <= hi(2)) then
              
                 p => dataptr(pf, i, j)
                 lo(:) = lwb(get_box(pf, i, j))
                 hi(:) = upb(get_box(pf, i, j))

                 jj = jloc*rr
              
                 ! loop over all of the zones in the patch.  Here, we convert
                 ! the cell-centered indices at the current level into the
                 ! corresponding RANGE on the finest level, and test if we've
                 ! stored data in any of those locations.  If we haven't then
                 ! we store this level's data and mark that range as filled.
                 do kk = lo(3), hi(3)
                    do ii = lo(1), hi(1)

                       if ( any(imask(ii*r1:(ii+1)*r1-1, &
                                      kk*r1:(kk+1)*r1-1) ) ) then
                 
                          ! since we've only bound one component to pf, we 
                          ! index p with 1 for the component
                          if (do_log) then
                             slicedata(ii*r1:(ii+1)*r1-1, &
                                       kk*r1:(kk+1)*r1-1) = log10(p(ii,jj,kk,1))
                          else
                             slicedata(ii*r1:(ii+1)*r1-1, &
                                       kk*r1:(kk+1)*r1-1) = p(ii,jj,kk,1)
                          endif

                          imask(ii*r1:(ii+1)*r1-1, &
                                kk*r1:(kk+1)*r1-1) = .false.
                 
                       end if

                    end do
                 enddo

              endif


           case (3)

              ! produce a slice perpendicular to the z-axis

              ! if the current patch stradles our slice plane, then get a data
              ! pointer to it
              if (rr*kloc >= lo(3) .and. rr*kloc <= hi(3)) then
              
                 p => dataptr(pf, i, j)
                 lo(:) = lwb(get_box(pf, i, j))
                 hi(:) = upb(get_box(pf, i, j))

                 kk = kloc*rr

                 ! loop over all of the zones in the patch.  Here, we convert
                 ! the cell-centered indices at the current level into the
                 ! corresponding RANGE on the finest level, and test if we've
                 ! stored data in any of those locations.  If we haven't then
                 ! we store this level's data and mark that range as filled.
                 do jj = lo(2), hi(2)
                    do ii = lo(1), hi(1)
                 
                       if ( any(imask(ii*r1:(ii+1)*r1-1, &
                                      jj*r1:(jj+1)*r1-1) ) ) then
                 
                          ! since we've only bound one component to pf, we 
                          ! index p with 1 for the component
                          if (do_log) then
                             slicedata(ii*r1:(ii+1)*r1-1, &
                                       jj*r1:(jj+1)*r1-1) = log10(p(ii,jj,kk,1))
                          else
                             slicedata(ii*r1:(ii+1)*r1-1, &
                                       jj*r1:(jj+1)*r1-1) = p(ii,jj,kk,1)
                          endif

                          imask(ii*r1:(ii+1)*r1-1, &
                                jj*r1:(jj+1)*r1-1) = .false.
                 
                       end if

                    enddo
                 enddo

              endif

           end select

           call fab_unbind(pf, i, j)

        enddo

        ! adjust r1 for the next lowest level
        if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
     
     end do

     deallocate(imask)

     !-------------------------------------------------------------------------
     ! output the slice as an image
     !-------------------------------------------------------------------------

     ! get basename of file
     indslsh = index(pltfile, '/', back = .TRUE.)
     if ( indslsh /= 0 ) then
        ofname = pltfile(indslsh+1:)
     else
        ofname = pltfile
     end if

     ! scale the slice data to range from 0:255, taking into account the 
     ! data limits
     intdata(:,:) = max(min(int(255.999*(slicedata(:,:) - gmn)/(gmx-gmn)),255),0)

     select case (ndir)
     case (1) 
        ndir_string = "YZ"
     case (2)
        ndir_string = "XZ"
     case (3)
        ndir_string = "XY"
     end select

     call store_ppm(trim(ofname) // "." // trim(compname) // "." // ndir_string // ".ppm", intdata(:,:), r, g, b)
  
     deallocate(slicedata,intdata)

  enddo

  call destroy(pf)

end program fsnapshot3d
