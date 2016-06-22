! produce an image of a 2-d dataset

program fsnapshot2d

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use filler_module
  use ppm_util_module

  implicit none

  integer, parameter :: NCOLOR = 256
  integer r(NCOLOR), g(NCOLOR), b(NCOLOR), a(NCOLOR)

  real(kind=dp_t) :: gmx, gmn
  real(kind=dp_t) :: def_mx, def_mn
  logical :: ldef_mx, ldef_mn, do_log

  type(plotfile) :: pf
  integer :: unit

  integer :: i

  integer :: comp
  character(len=64) :: compname

  integer :: plo(2), phi(2)

  real(kind=dp_t), allocatable :: c_fab(:,:,:)
  integer, allocatable :: intdata(:,:)

  character(len=256) :: pltfile
  character(len=256) :: ofname

  integer :: indslsh

  integer :: narg, farg
  character(len=256) ::  fname

  character(len=128) :: phome
  character(len=256) :: pfname

  integer :: max_level

  unit = unit_new()

  ! set the defaults
  pltfile = ''
  compname = "density"
  max_level = -1
  ldef_mx = .false.
  ldef_mn = .false.
  do_log = .false.


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

     case ('--max_level')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) max_level

     case ('-l', '--log')
        do_log = .true.

     case default
        exit

     end select

     farg = farg + 1
  end do

  if ( len_trim(pltfile) == 0 .OR. len_trim(compname) == 0) then
     print *, "produce an image of a 2-d dataset"
     print *, " "
     print *, "usage: fsnapshot2d args"
     print *, "args [-p    |--pltfile]   plotfile   : plot file directory (required)"
     print *, "     [-cname|--compname]  name       : variable to plot (default: density)"
     print *, "     --palette pfname                : use the file pfname as the Palette"
     print *, "                                       (default ~/amrvis.Palette)"
     print *, " "
     print *, "     -m val                          : set the minimum value of the data to val"
     print *, "     -M val                          : set the maximum value of the data to val"
     print *, "     [-l|--log]                      : toggle log plot"
     stop
  end if
     

  ! make sure we have valid options set
  if (do_log) then
     if (ldef_mx .and. def_mx < ZERO) call bl_error("ERROR: log plot specified with negative maximum")
     if (ldef_mn .and. def_mn < ZERO) call bl_error("ERROR: log plot specified with negative minimum")
  endif

  ! get the palette
  if ( pfname /= '' ) then
     call load_palette(pfname, r, g, b, a)
  end if

  ! build the plotfile to get the level and component information
  call build(pf, pltfile, unit)


  ! figure out the variable indices
  comp = plotfile_var_index(pf, compname)
  if (comp < 0) then
     call bl_error("ERROR: variable not found in plotfile")
  endif

  ! set the level to visualize -- default to the finest level
  if (max_level == -1) then
     max_level = pf%flevel
  else
     if (max_level < 1 .OR. max_level > pf%flevel) then
        call bl_error("ERROR: specified level not allowed")
     endif
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

  ! for 2-d, we will put the dataset onto a uniform grid 
  plo = lwb(plotfile_get_pd_box(pf, max_level))
  phi = upb(plotfile_get_pd_box(pf, max_level))
  allocate(intdata(plo(1):phi(1), plo(2):phi(2)))
  allocate(c_fab(plo(1):phi(1), plo(2):phi(2), 1))

  call blow_out_to_fab(c_fab, plo, pf, (/comp/), max_level)

  ! change the scale if need be
  if (do_log) then
     c_fab(:,:,:) = log10(c_fab(:,:,:))
  endif


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

  ! scale the data to range from 0:255, taking into account the
  ! data limits
  intdata(:,:) = max(min(int(255.999*(c_fab(:,:,1) - gmn)/(gmx-gmn)),255),0)

  call store_ppm(trim(ofname) // "." // trim(compname) // '.ppm', intdata(:,:), r, g, b)

  deallocate(c_fab,intdata)

  call destroy(pf)
  
end program fsnapshot2d
