program fwritecontents2d
  use BoxLib
  implicit none

  character (len=256) :: infile, outfile, input_field, component
  integer :: i, j, narg, farg, nx, ny, nz, ierr
  character (len=256) :: fname
  real*8, allocatable :: array(:,:)
  
  !---------------------------------------------------------------------------
  ! process the command line arguments

  narg = command_argument_count()

  ! defaults
  infile = ""
  outfile = ""

  farg = 1
  do while (farg <= narg)
     call get_command_argument(farg, value = fname)
     
     select case (fname)

     case ('-i', '--infile')
        farg = farg + 1
        call get_command_argument(farg, value = infile)

     case ('-o', '--outfile')
        farg = farg + 1
        call get_command_argument(farg, value = outfile)

     case ('-c', '--component')
        farg = farg + 1
        call get_command_argument(farg, value = component)

     case default
        exit

     end select
     farg = farg + 1
  enddo

  if ((len_trim(infile) == 0).or.(len_trim(outfile) == 0)) then
     stop "All arguments are needed"
  endif

  call boxlib_initialize()
  
  call fplotfile_get_size(infile, nx, ny, nz)
  write(*,*) "SIZE=", nx, ny, nz
  
  if(nz>1) stop "Only works for 2d plotfiles"
  
  allocate(array(nx,ny))
  
  call fplotfile_get_data_2d(infile, component, array, nx, ny, ierr)
  
  open(unit=11, file=trim(outfile), status="unknown", action="write")  
  do j=1,ny
     do i=1,nx
        write(11,*) real(array(i,j))
     end do
  end do  
  close(11)
  
  deallocate(array)

contains

!------------------------------------------------------------------------------
! fplotfile_get_size
!------------------------------------------------------------------------------
subroutine fplotfile_get_size(pltfile, nx, ny, nz)

  use bl_IO_module
  use plotfile_module

  implicit none

  character (len=*), intent(in) :: pltfile
  integer, intent(out) :: nx, ny, nz

!f2py intent(in) :: pltfile
!f2py intent(out) :: nx, ny, nz
!f2py note return the dimensions of the data at the finest level of refinement
  type(plotfile) :: pf
  integer :: unit

  integer, allocatable :: flo(:), fhi(:)

  integer :: max_level


  ! build the plotfile to get the level information
  unit = unit_new()
  call build(pf, pltfile, unit)

  max_level = pf%flevel

  if (pf%dim == 1) then
     
     ! 1-d -- return maximum size of finest level
     allocate(flo(1), fhi(1))

     flo = lwb(plotfile_get_pd_box(pf, max_level))
     fhi = upb(plotfile_get_pd_box(pf, max_level))

     nx = fhi(1) - flo(1) + 1
     ny = -1
     nz = -1

  else if (pf%dim == 2) then
     
     ! 2-d -- return maximum size of finest level
     allocate(flo(2), fhi(2))

     flo = lwb(plotfile_get_pd_box(pf, max_level))
     fhi = upb(plotfile_get_pd_box(pf, max_level))

     nx = fhi(1) - flo(1) + 1
     ny = fhi(2) - flo(2) + 1
     nz = -1

  else if (pf%dim == 3) then

     ! 3-d -- return maximum size of finest level
     allocate(flo(3), fhi(3))

     flo = lwb(plotfile_get_pd_box(pf, max_level))
     fhi = upb(plotfile_get_pd_box(pf, max_level))

     nx = fhi(1) - flo(1) + 1
     ny = fhi(2) - flo(2) + 1
     nz = fhi(3) - flo(3) + 1

  endif

  deallocate(flo, fhi)
  call destroy(pf)

end subroutine fplotfile_get_size


!------------------------------------------------------------------------------
! fplotfile_get_time
!------------------------------------------------------------------------------
subroutine fplotfile_get_time(pltfile, time)

  use bl_IO_module
  use plotfile_module

  implicit none

  character (len=*), intent(in) :: pltfile
  double precision, intent(out) :: time

!f2py intent(in) :: pltfile
!f2py intent(out) :: time

  type(plotfile) :: pf
  integer :: unit


  ! build the plotfile to get the level information
  unit = unit_new()
  call build(pf, pltfile, unit)

  time = pf%tm

  call destroy(pf)

end subroutine fplotfile_get_time


!------------------------------------------------------------------------------
! fplotfile_get_nvar
!------------------------------------------------------------------------------
subroutine fplotfile_get_nvar(pltfile, nvar)

  use bl_IO_module
  use plotfile_module

  implicit none

  character (len=*), intent(in) :: pltfile
  integer, intent(out) :: nvar

!f2py intent(in) :: pltfile
!f2py intent(out) :: nvar

  type(plotfile) :: pf
  integer :: unit


  ! build the plotfile to get the level information
  unit = unit_new()
  call build(pf, pltfile, unit)

  nvar = pf%nvars

  call destroy(pf)

end subroutine fplotfile_get_nvar


!------------------------------------------------------------------------------
! fplotfile_get_varnames
!------------------------------------------------------------------------------
subroutine fplotfile_get_varinfo(pltfile, ivar, varname, varmin, varmax, ierr)

  use bl_IO_module
  use plotfile_module

  implicit none

  character (len=*) , intent(in) :: pltfile
  integer           , intent(in ) :: ivar
  character (len=64), intent(out) :: varname
  real              , intent(out) :: varmin
  real              , intent(out) :: varmax
  integer           , intent(out) :: ierr

!f2py intent(in) :: pltfile
!f2py intent(in) :: ivar
!f2py intent(out) :: varname
!f2py intent(out) :: varmin
!f2py intent(out) :: varmax
!f2py intent(out) :: ierr

  type(plotfile) :: pf
  integer :: unit

  ierr = 0

  ! build the plotfile to get the level information
  unit = unit_new()
  call build(pf, pltfile, unit)

  if (ivar > pf%nvars .or. ivar < 1) then
     print *, "ERROR: component ", ivar, " outside of valid range."
     ierr = 1
     return
  endif

  varname = pf%names(ivar)

  varmin = plotfile_minval(pf,ivar,pf%flevel)
  varmax = plotfile_maxval(pf,ivar,pf%flevel)

  call destroy(pf)

end subroutine fplotfile_get_varinfo


!------------------------------------------------------------------------------
! fplotfile_get_limits
!------------------------------------------------------------------------------
subroutine fplotfile_get_limits(pltfile, xmin, xmax, ymin, ymax, zmin, zmax)

  use bl_IO_module
  use plotfile_module

  implicit none

  character (len=*), intent(in) :: pltfile
  double precision, intent(out) :: xmin, xmax, ymin, ymax, zmin, zmax

!f2py intent(in) :: pltfile
!f2py intent(out) :: xmin, xmax, ymin, ymax, zmin, zmax

  type(plotfile) :: pf
  integer :: unit

  integer :: max_level


  ! build the plotfile to get the level information
  unit = unit_new()
  call build(pf, pltfile, unit)

  xmin = pf%plo(1)
  xmax = pf%phi(1)

  if (pf%dim >= 2) then
     ymin = pf%plo(2)
     ymax = pf%phi(2)
  else
     ymin = -1.0
     ymax = -1.0
  endif

  if (pf%dim == 3) then
     zmin = pf%plo(3)
     zmax = pf%phi(3)
  else
     zmin = -1.0
     zmax = -1.0
  endif

  call destroy(pf)

end subroutine fplotfile_get_limits


!------------------------------------------------------------------------------
! fplotfile_get_data_1d
!------------------------------------------------------------------------------
subroutine fplotfile_get_data_1d(pltfile, component, mydata, x, nx_max, nx, ierr)

  ! return a single variable in mydata with coordinate locations x --
  ! it may not be evenly gridded due to AMR.  nx_max is the size of
  ! the mydata and x arrays as allocated (hidden in this interface)
  ! while nx is the number of points actually stored, since some may
  ! be at lower resolution.

  use bl_space
  use bl_constants_module, ONLY: ZERO, HALF
  use bl_IO_module
  use plotfile_module
  use filler_module
  use sort_d_module

  implicit none

  character (len=*), intent(in) :: pltfile
  character (len=*), intent(in) :: component
  integer          , intent(in) :: nx_max
  double precision , intent(inout) :: mydata(nx_max)
  double precision , intent(inout) :: x(nx_max)
  integer          , intent(out) :: nx
  integer, intent(out) :: ierr

!f2py intent(in) :: pltfile, component
!f2py intent(in,out) :: mydata
!f2py intent(in,out) :: x
!f2py intent(hide) :: nx_max
!f2py intent(out) :: nx
!f2py intent(out) :: ierr

  type(plotfile) :: pf
  integer :: unit

  integer :: comp

  integer :: flo(1), fhi(1)

  integer :: max_level
  
  integer :: i, j, ii

  integer :: cnt

  real(kind=dp_t), pointer :: p(:,:,:,:)
  real(kind=dp_t), allocatable :: x_tmp(:), mydata_tmp(:)
  logical, allocatable :: imask(:)
  integer, allocatable :: isv(:)
 
  integer :: r1, rr
  real(kind=dp_t) :: dx(MAX_SPACEDIM)
  real(kind=dp_t) :: xmin

  ierr = 0

  ! build the plotfile to get the level and component information
  unit = unit_new()
  call build(pf, pltfile, unit)

  ! figure out the variable indices
  comp = plotfile_var_index(pf, component)
  if (comp < 0) then
     print *, "ERROR: component ", trim(component), " not found in plotfile"
     ierr = 1
     return
  endif

  max_level = pf%flevel


  ! get dx for the coarse level
  dx = plotfile_get_dx(pf, 1)


  ! domain limit
  xmin = pf%plo(1)
    

  ! domain index limits on the fine level
  flo = lwb(plotfile_get_pd_box(pf, max_level))
  fhi = upb(plotfile_get_pd_box(pf, max_level))


  ! imask will be set to false if we've already output the data.                
  ! Note, imask is defined in terms of the finest level.  As we loop            
  ! over levels, we will compare to the finest level index space to             
  ! determine if we've already output here  
  allocate(imask(flo(1):fhi(1)))
  imask = .true.


  ! in general, we can get the data out of order.  Allocate an array
  ! of indices for sorting later
  allocate(isv(nx_max))


  ! temporary storage for the unsorted data
  allocate(x_tmp(nx_max))
  allocate(mydata_tmp(nx_max))

  cnt = 0

  ! r1 is the factor between the current level grid spacing and the             
  ! FINEST level                                                                
  r1  = 1

  do i = pf%flevel, 1, -1

     ! rr is the factor between the COARSEST level grid spacing and             
     ! the current level                                                        
     rr = product(pf%refrat(1:i-1,1))

     do j = 1, nboxes(pf, i)

        ! read in the data 1 patch at a time, single component only
        call fab_bind_comp_vec(pf, i, j, (/comp/) )

        p => dataptr(pf, i, j)

        do ii = lbound(p,dim=1), ubound(p,dim=1)
           if ( any(imask(ii*r1:(ii+1)*r1-1) ) ) then

              cnt = cnt + 1

              x_tmp(cnt) = xmin + (ii + HALF)*dx(1)/rr
              mydata_tmp(cnt) = p(ii,1,1,1)

              imask(ii*r1:(ii+1)*r1-1) = .false.
           endif
        enddo

        call fab_unbind(pf, i, j)

     enddo

     ! adjust r1 for the next lowest level                                      
     if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)

  enddo

  ! sort the data based on the coordinates
  call sort(x_tmp(1:cnt),isv(1:cnt))

  x(:) = ZERO
  mydata(:) = ZERO

  do i = 1, cnt
     x(i)      = x_tmp(isv(i))
     mydata(i) = mydata_tmp(isv(i))
  enddo

  nx = cnt

  deallocate(imask)
  deallocate(isv)
  deallocate(x_tmp,mydata_tmp)

  call destroy(pf)

end subroutine fplotfile_get_data_1d


!------------------------------------------------------------------------------
! fplotfile_get_data_2d
!------------------------------------------------------------------------------
subroutine fplotfile_get_data_2d(pltfile, component, mydata, nx, ny, ierr)

  use bl_IO_module
  use plotfile_module
  use filler_module

  implicit none

  character (len=*), intent(in) :: pltfile
  character (len=*), intent(in) :: component
  integer, intent(in) :: nx, ny
  double precision, intent(inout) :: mydata(nx, ny)
  integer, intent(out) :: ierr

!f2py intent(in) :: pltfile, component
!f2py intent(in,out) :: mydata
!f2py intent(hide) :: nx, ny
!f2py intent(out) :: ierr

  type(plotfile) :: pf
  integer ::unit

  integer :: comp

  integer :: flo(2), fhi(2)
  real(kind=dp_t), allocatable :: c_fab(:,:,:)

  integer :: max_level

  ierr = 0
  
  ! build the plotfile to get the level and component information
  unit = unit_new()
  call build(pf, pltfile, unit)

  ! figure out the variable indices
  comp = plotfile_var_index(pf, component)
  if (comp < 0) then
     print *, "ERROR: component ", trim(component), " not found in plotfile"
     ierr = 1
     return
  endif

  max_level = pf%flevel

  ! for 2-d, we will put the dataset onto a uniform grid 
  flo = lwb(plotfile_get_pd_box(pf, max_level))
  fhi = upb(plotfile_get_pd_box(pf, max_level))
  allocate(c_fab(flo(1):fhi(1), flo(2):fhi(2), 1))

  if ( (nx /= (fhi(1) - flo(1) + 1)) .or. &
       (ny /= (fhi(2) - flo(2) + 1)) ) then
     print *, "ERROR: input array wrong dimensions"
     deallocate(c_fab)
     ierr = 1
     return
  endif

  call blow_out_to_fab(c_fab, flo, pf, (/comp/), max_level)
  
  mydata(:,:) = c_fab(:,:,1)

  deallocate(c_fab)

  call destroy(pf)

end subroutine fplotfile_get_data_2d


!------------------------------------------------------------------------------
! fplotfile_get_data_3d
!------------------------------------------------------------------------------
subroutine fplotfile_get_data_3d(pltfile, component, &
                                 indir, origin, mydata, m, n, ierr)

  use bl_constants_module, ONLY: ZERO
  use bl_IO_module
  use plotfile_module
  use filler_module

  implicit none

  character (len=*), intent(in) :: pltfile
  character (len=*), intent(in) :: component
  integer, intent(in) :: indir     ! normal direction
  integer, intent(in) :: origin    ! pass through origin (1), or center
  integer, intent(in) :: m, n
  double precision, intent(inout) :: mydata(m, n)
  integer, intent(out) :: ierr

!f2py intent(in) :: pltfile, component
!f2py intent(in) :: indir
!f2py intent(in) :: origin
!f2py intent(in,out) :: mydata
!f2py intent(hide) :: m, n
!f2py intent(out) :: ierr

  real (kind=dp_t), allocatable :: slicedata(:,:)

  type(plotfile) :: pf
  integer ::unit

  integer :: comp

  integer ::  lo(3),  hi(3)
  integer :: flo(3), fhi(3)

  integer :: rr, r1

  integer :: i, j
  integer :: ii, jj, kk
  integer :: iloc, jloc, kloc

  real(kind=dp_t) :: dx(3), dx_fine(3)

  real(kind=dp_t), pointer :: p(:,:,:,:)

  logical, allocatable :: imask(:,:)

  integer :: max_level


  ierr = 0

  ! build the plotfile to get the level and component information
  unit = unit_new()
  call build(pf, pltfile, unit)

  ! figure out the variable indices
  comp = plotfile_var_index(pf, component)
  if (comp < 0) then
     print *, "ERROR: component ", trim(component), " not found in plotfile"
     ierr = 1
     return
  endif


  ! get the index bounds and dx for the coarse level.  Note, lo and
  ! hi are ZERO based indicies
  lo = lwb(plotfile_get_pd_box(pf, 1))
  hi = upb(plotfile_get_pd_box(pf, 1))

  dx = plotfile_get_dx(pf, 1)

  ! get the index bounds for the finest level
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

  dx_fine = minval(plotfile_get_dx(pf, pf%flevel))


  ! sanity checks
  select case (indir)

  case (1)
     if ( (m /= (fhi(2) - flo(2) + 1)) .or. &
          (n /= (fhi(3) - flo(3) + 1)) ) then
        print *, "ERROR: input array wrong dimensions"
        ierr = 1
        return
     endif

  case (2)
     if ( (m /= (fhi(1) - flo(1) + 1)) .or. &
          (n /= (fhi(3) - flo(3) + 1)) ) then
        print *, "ERROR: input array wrong dimensions"
        ierr = 1
        return
     endif

  case (3)
     if ( (m /= (fhi(1) - flo(1) + 1)) .or. &
          (n /= (fhi(2) - flo(2) + 1)) ) then
        print *, "ERROR: input array wrong dimensions"
        ierr = 1
        return
     endif

  case default

     print *, "ERROR: invalid indir"
     ierr = 1
     return

  end select



  ! determine the location of the slices -- either through the center
  ! or through the origin -- these are with respect to the coarse grid
  if (origin == 1) then
     iloc = lo(1)
     jloc = lo(2)
     kloc = lo(3)
  else
     iloc = (hi(1)-lo(1)+1)/2 + lo(1)
     jloc = (hi(2)-lo(2)+1)/2 + lo(2)
     kloc = (hi(3)-lo(3)+1)/2 + lo(3)
  endif


  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here
  select case (indir)

  case (1)
     allocate(imask(flo(2):fhi(2),flo(3):fhi(3)))
     allocate(slicedata(flo(2):fhi(2),flo(3):fhi(3)))

  case (2)
     allocate(imask(flo(1):fhi(1),flo(3):fhi(3)))
     allocate(slicedata(flo(1):fhi(1),flo(3):fhi(3)))

  case (3)
     allocate(imask(flo(1):fhi(1),flo(2):fhi(2)))
     allocate(slicedata(flo(1):fhi(1),flo(2):fhi(2)))

  end select

  imask(:,:) = .true.


  slicedata(:,:) = ZERO


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


        select case (indir)

        case (1)
           ! produce a slice perpendicular to the x-axis

           ! if the current patch stradles our slice plane, then get a data  
           ! pointer to it                                                   
           if (rr*iloc >= lo(1) .and. rr*iloc <= hi(1)) then

              p => dataptr(pf, i, j)
              ii = iloc*rr

              ! loop over all of the zones in the patch.  Here, we convert   
              ! the cell-centered indices at the current level into the      
              ! corresponding RANGE on the finest level, and test if we've   
              ! stored data in any of those locations.  If we haven't then   
              ! we store this level's data and mark that range as filled.    
              do kk = lbound(p,dim=3), ubound(p,dim=3)
                 do jj = lbound(p,dim=2), ubound(p,dim=2)

                    if ( any(imask(jj*r1:(jj+1)*r1-1, &
                                   kk*r1:(kk+1)*r1-1) ) ) then

                       ! since we've only bound one component to pf, we      
                       ! index p with 1 for the component                    
                       slicedata(jj*r1:(jj+1)*r1-1, &
                                 kk*r1:(kk+1)*r1-1) = p(ii,jj,kk,1)

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
              jj = jloc*rr

              ! loop over all of the zones in the patch.  Here, we convert   
              ! the cell-centered indices at the current level into the      
              ! corresponding RANGE on the finest level, and test if we've   
              ! stored data in any of those locations.  If we haven't then   
              ! we store this level's data and mark that range as filled.    
              do kk = lbound(p,dim=3), ubound(p,dim=3)
                 do ii = lbound(p,dim=1), ubound(p,dim=1)

                    if ( any(imask(ii*r1:(ii+1)*r1-1, &
                                   kk*r1:(kk+1)*r1-1) ) ) then

                       ! since we've only bound one component to pf, we      
                       ! index p with 1 for the component                    
                       slicedata(ii*r1:(ii+1)*r1-1, &
                                 kk*r1:(kk+1)*r1-1) = p(ii,jj,kk,1)

                       imask(ii*r1:(ii+1)*r1-1, &
                             kk*r1:(kk+1)*r1-1) = .false.

                    endif

                 end do
              enddo
              
           endif

        case (3)
           ! produce a slice perpendicular to the z-axis

           ! if the current patch stradles our slice plane, then get a data  
           ! pointer to it                                                   
           if (rr*kloc >= lo(3) .and. rr*kloc <= hi(3)) then

              p => dataptr(pf, i, j)
              kk = kloc*rr

              ! loop over all of the zones in the patch.  Here, we convert   
              ! the cell-centered indices at the current level into the      
              ! corresponding RANGE on the finest level, and test if we've   
              ! stored data in any of those locations.  If we haven't then   
              ! we store this level's data and mark that range as filled.    
              do jj = lbound(p,dim=2), ubound(p,dim=2)
                 do ii = lbound(p,dim=1), ubound(p,dim=1)

                    if ( any(imask(ii*r1:(ii+1)*r1-1, &
                                   jj*r1:(jj+1)*r1-1) ) ) then

                       ! since we've only bound one component to pf, we      
                       ! index p with 1 for the component                    
                       slicedata(ii*r1:(ii+1)*r1-1, &
                                 jj*r1:(jj+1)*r1-1) = p(ii,jj,kk,1)

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

  mydata(:,:) = slicedata(:,:)

  deallocate(imask)
  deallocate(slicedata)

  call destroy(pf)

end subroutine fplotfile_get_data_3d

end program
