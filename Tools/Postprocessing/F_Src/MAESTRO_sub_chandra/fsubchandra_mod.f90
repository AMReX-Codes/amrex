! Module for analysis routines for the 3-d sub_chandra MAESTRO problem 
! This code borrows from fsedov3d_sph.f90 and fwdconvect.f90

module subchandra

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module

  implicit none
  private 

  !Type definitions 
  !A collection of state variable component indexes
  type state_comps
    integer :: dens_comp = -1, temp_comp = -1, pi_comp = -1, p0_comp = -1
    integer :: magv_comp = -1, h_comp = -1
    integer :: XC12_comp = -1, XO16_comp = -1, XHe4_comp = -1, rhopert_comp = -1
    integer :: XdotC12_comp = -1, XdotO16_comp = -1, XdotHe4_comp = -1
    integer :: tpert_comp = -1, enuc_comp = -1, xvel_comp = -1
    integer :: yvel_comp = -1, zvel_comp = -1, s_comp = -1
    integer :: spert_comp = -1
  end type

  !Geometric attributes
  type geometry
    real(kind=dp_t) :: dx(MAX_SPACEDIM)
    integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)
    real(kind=dp_t) :: xctr, yctr, zctr
    !Note this may not be portable.  Technically f90/95 doesn't allow
    !allocatables in derived types, but in practice most compilers implemented
    !the f2003 fix that allows this.
    real(kind=dp_t), allocatable :: r(:) 
    real(kind=dp_t) :: dx_fine
  end type geometry

  !Collection of radial averages
  type radial_averages
    !For convenience and readability we store nbins
    integer :: nbins = -1
    ! ncount keeps track of how many zones were written into
    ! each bin
    integer, allocatable :: ncount(:)
    !binned averaging data for various physical quantities
    real(kind=dp_t), allocatable :: dens_avg_bin(:), dens_rms_bin(:)
    real(kind=dp_t), allocatable :: temp_avg_bin(:), temp_rms_bin(:)
    real(kind=dp_t), allocatable :: pres_avg_bin(:), pres_rms_bin(:)
    real(kind=dp_t), allocatable :: magv_avg_bin(:), enuc_avg_bin(:)
    real(kind=dp_t), allocatable :: XC12_avg_bin(:), XO16_avg_bin(:), XHe4_avg_bin(:)
    real(kind=dp_t), allocatable :: XdotC12_avg_bin(:), XdotO16_avg_bin(:), XdotHe4_avg_bin(:)
    real(kind=dp_t), allocatable :: entropy_avg_bin(:), entropy_rms_bin(:)
  end type radial_averages

  !Collection of global results based on evaluation of
  !the entire computational domain
  type globals
    !Peak temperature attributes
    real(kind=dp_t) :: T_peak = 0.0_dp_t
    real(kind=dp_t) :: xloc_Tpeak = 0.0_dp_t, yloc_Tpeak = 0.0_dp_t
    real(kind=dp_t) :: zloc_Tpeak = 0.0_dp_t, R_Tpeak = 0.0_dp_t
    real(kind=dp_t) :: vx_Tpeak = 0.0_dp_t, vy_Tpeak = 0.0_dp_t
    real(kind=dp_t) :: vz_Tpeak = 0.0_dp_t, vr_Tpeak = 0.0_dp_t
    
    !Peak energy generation attributes
    real(kind=dp_t) :: enuc_peak = 0.0_dp_t
    real(kind=dp_t) :: xloc_enucpeak = 0.0_dp_t, yloc_enucpeak = 0.0_dp_t
    real(kind=dp_t) :: zloc_enucpeak = 0.0_dp_t, R_enucpeak = 0.0_dp_t
  end type globals
  
  !Store the salient properties of a single hotspot, defined here as a
  !hot cell.
  type hotspot
    !Temperature and density
    real(kind=dp_t) :: temp = -1.0_dp_t
    real(kind=dp_t) :: rho = -1.0_dp_t

    !Location info
    real(kind=dp_t) :: x = -1.0_dp_t, y = -1.0_dp_t, z = -1.0_dp_t

    !Cell info
    integer :: level = -1
  end type hotspot

  !Histogram data for average temperatures at various lengthscales
  type temp_hist
     !Lengthscales for each level (nlevs)
     real(kind=dp_t), allocatable :: ell(:)

     !Temperature limits and delta
     real(kind=dp_t) :: min_temp = 9.0d6
     real(kind=dp_t) :: max_temp = 9.0d8
     real(kind=dp_t) :: dT = 1.0d5
     
     !Histogram temperature data (nlevs, nbins)
     integer, allocatable :: tavg(:,:)   !Use the temperature of the cell
     integer, allocatable :: tfeavg(:,:) !Derive temp from EoS
     integer :: nbins = nint((9.0d8 - 9.0d6)/1.0d5)
  end type temp_hist


  !A min-heap of the hottest hotspots with metadata about the heap.
  !The min-heap implementation allows for quick (log(n)) removal of the
  !lowest temperature hotspot.
  !  Heap implementation:
  !    lptrs acts as the heap with pseudopointers (integer values representing
  !    indices in list).  The heap is implemented as a binary tree with the root
  !    at lptrs(1).  Any element i has a left child at 2i, a right child at
  !    2i+1, and a parent at floor(i/2).  All elements in a min-heap binary tree
  !    satisfy the heap property that they're less than or equal to their
  !    children.  The tree must also be complete (all levels filled) except for
  !    the last level, which is filled from left to right.
  !
  !    Schematic where a<=b,c ; b<=d,e ; c<=f:
  !    |a   |b |c |d |e |f |  |
  !     root la ra lb rb lc rc
  !
  !           a
  !           |
  !         b---c
  !         |   |
  !        d-e f- 
  !
  !TODO: The use of pseudopointers prevents us from having to move around
  !elements in the list of type(hotspot)'s, but I don't think this minor
  !boost in performance is worth the much less readable code.
  type hheap
    !Max, min value
    real(kind=dp_t) :: max_temp = -1.0_dp_t, min_temp = -1.0_dp_t

    !Index of max and min value, index of available element
    integer :: maxdex = -1, mindex = -1, adex = 1

    !List and integer array of pseudopointers for heap management
    type(hotspot), allocatable :: list(:)
    integer, allocatable :: lptrs(:)
    integer :: length = -1
  end type hheap

  !Public subroutines
  public :: parse_args, init_comps, init_geometry, init_averages
  public :: analyze, writeout
  !Public types
  public :: globals, radial_averages, state_comps, geometry, hheap, hotspot
  public :: temp_hist

contains
  !Parse command line args, initialize variables
  subroutine parse_args(slicefile, pltfile, globals_only, fullstar, hpcnt, pf)

    implicit none

    !Arguments
    real(kind=dp_t),    intent(out) :: hpcnt
    logical,            intent(out) :: globals_only, fullstar
    character(len=256), intent(out) :: slicefile, pltfile
    type(plotfile),     intent(out) :: pf

    !Local data
    integer :: pffd, narg, iarg, indslsh, ierr
    character(len=256) :: fname, hpstring

    ! set the defaults
    slicefile = ''
    globals_only = .false.
    fullstar = .false.
    hpcnt = -1.0_dp_t !hpcnt < 0 indicates we don't want to do hotspot statistics

    ! read user's command line options 
    narg = command_argument_count()
    iarg = 1
    do while ( iarg <= narg )
      call get_command_argument(iarg, value = fname)

      select case (fname)
      case ('-s', '--slicefile')
        iarg = iarg + 1
        call get_command_argument(iarg, value = slicefile)

      case ('-h', '--do-hotspots')
        iarg = iarg + 1
        call get_command_argument(iarg, value = hpstring)
        read(unit=hpstring, fmt=*, iostat=ierr) hpcnt
        if(ierr > 0) then
          call usage()
          stop 'Invalid argument for --do-hotspots, must be a number!'
        endif

      case ('--globals_only', '--globals-only')
        globals_only = .true.

      case ('-f', '--full-star')
        fullstar = .true.

      case default
        !if it's not the last argument (plotfile) and we don't recognize
        !from above cases, then we have incorrect usage
        if(iarg .ne. narg) then
          call usage()
        end if
        exit
      end select

      iarg = iarg + 1
    end do

    ! print usage if --slicefile misused (i.e. it's case was able to
    ! increment iarg beyond narg)
    if (iarg > narg) then
      call usage()
    end if

    ! plotfile is the last argument, grab it
    call get_command_argument(iarg, value = pltfile)

    ! if slicefile not defined, default to <pltfile>.slice
    if ( len_trim(slicefile) == 0 ) then

      ! get basename of file
      indslsh = index(pltfile(:len_trim(pltfile)-1), '/', back = .TRUE.)  !The slicing of pltfile is to cut off any trailing '/'

      if ( indslsh /= 0 ) then
        slicefile = trim(pltfile(indslsh+1:)) // ".slice"
      else
        slicefile = trim(pltfile) // ".slice"
      end if

    endif

    !Build plotfile object 
    pffd = unit_new()
    call build(pf, pltfile, pffd)
  end subroutine parse_args

  !Print proper usage and stop program
  subroutine usage
    implicit none

    print *, "usage: fsubchandra [args] plotfile"
    print *, "args: [-s|--slicefile]   slice file : specify output file"
    print *, "      [-h|--do-hotspots] pcnt       : calculate hotspot statistics"
    print *, "                                      for hottest pcnt of cells"
    print *, "                                      WARNING: can generate lots of data."
    print *, "      [-f|--full-star]              : the plotfile comes from a full star simulation (not octant)"
    print *, "      --globals-only                : only output global quantities "
    stop
  end subroutine usage

  !Add a new hotspot to hheap and update metadata
  !If the list is full, then replace the lowest temperature
  !hotspot with this one
  !Assumptions:
  ! +hh is properly initialized before first call to hheap_add
  subroutine hheap_add(hh, newhs)
    !Modules
    implicit none

    !Args
    type(hheap), intent(inout) :: hh
    type(hotspot), intent(in) :: newhs

    !Local
    integer, save :: it=0
    integer :: i, pdex, left, right, smallest, mydex, avail, swap
    type(hotspot) :: cur_hs
    real(kind=dp_t) :: parent_temp, my_temp

    !This is the currently available index
    avail = min(hh%adex, hh%length)
    
    !If the heap isn't full update adex
    if(hh%adex <= hh%length) then
      !Update adex, set pointer
      hh%adex = hh%adex + 1
      hh%lptrs(avail) = avail
    !Otherwise we need to remove the heap root/minimum temperature 
    !(maintaining all heap properties) before adding the new hotspot
    else
      !Sanity check
      if( (hh%maxdex .eq. -1) .or. (hh%mindex .eq. -1) ) then
        call bl_error("ERROR, hheap_add: this shouldn't happen")
      end if

      !Replace root with rightmost element of bottom level in tree.
      !  This insures we maintain the property that
      !  the binary tree is full except possibly for the last level
      !  which must always be filled from left to right.
      !hh%list(hh%lptrs(1)) = hh%list(hh%lptrs(hh%length))
      swap = hh%lptrs(1)
      hh%lptrs(1) = hh%lptrs(hh%length)

      !Set rightmost element to point to the old root list element.
      !We'll overwrite the list element this points to 
      !with the new hotspot later
      hh%lptrs(avail) = swap

      !Percolate the new root down 
      !until min-heap property is satisfied
      !Note that for this loop we're short one element so the highest
      !valid element is hh%length - 1.
      i = 1
      do while (.true.)
        left = 2*i
        right = 2*i+1

        !Find the index of the lowest temperature member of the set
        !(parent, left child, right child)
        smallest = i
        if(left < hh%length .and. &
            hh%list(hh%lptrs(left))%temp < hh%list(hh%lptrs(smallest))%temp) then
          smallest = left
        end if
        if(right < hh%length .and. &
            hh%list(hh%lptrs(right))%temp < hh%list(hh%lptrs(smallest))%temp) then
          smallest = right
        end if

        if(smallest == i) then
          !Heap property satisfied, exit loop
          exit
        else
          !Parent wasn't smallest, so swap parent with child and
          !update i for next loop cycle
          swap = hh%lptrs(i) 
          hh%lptrs(i) = hh%lptrs(smallest)
          hh%lptrs(smallest) = swap
          i = smallest
        end if
      end do
    end if

    !Add to list (overwriting the element previously pointed to by the root
    !if heap was full). 
    !Avail is either the last element if heap was full 
    !or the next available index is heap wasn't full.
    hh%list(hh%lptrs(avail)) = newhs

    !"Percolate" newly added value up
    mydex = avail
    do while(mydex > 0 .and. mydex <= avail) !within bounds
      !At root, stop
      if(mydex == 1) then
        exit
      end if 

      !Calculate info
      pdex = int(mydex/2)
      parent_temp = hh%list(hh%lptrs(pdex))%temp
      my_temp = hh%list(hh%lptrs(mydex))%temp

      !If parent is hotter, swap
      if(parent_temp > my_temp) then
        swap = hh%lptrs(mydex)
        hh%lptrs(mydex) = hh%lptrs(pdex)
        hh%lptrs(pdex) = swap
        mydex = pdex
      !Otherwise the heap property is satisfied, break out of loop
      else
        exit
      end if
    end do
    
    !Update extrema metadata
    if( (hh%maxdex .eq. -1) .or. (newhs%temp > hh%max_temp) ) then
      hh%max_temp = newhs%temp
      hh%maxdex = hh%lptrs(avail)
    end if
    hh%min_temp = hh%list(hh%lptrs(1))%temp
    hh%mindex = hh%lptrs(1)
  end subroutine hheap_add

  !Initialize heap of hotspots
  subroutine init_hheap(hh, len)
    !Modules
    implicit none

    !Args
    type(hheap), intent(out) :: hh
    integer, intent(in) :: len

    !NOTE! Here we do NOT follow the typical Maestro practice of starting
    !indexing from 0.  We do this because it makes the heap operations much
    !clearer to implement
    allocate(hh%list(1:len))
    allocate(hh%lptrs(1:len))
    !TODO: set all initial hotspots to some value?
    hh%length = len
    hh%adex = 1
  end subroutine init_hheap

  !Initialize component indices
  subroutine init_comps(pf, sc)
    implicit none

    type(plotfile), intent(in) :: pf
    type(state_comps), intent(out) :: sc

    sc%dens_comp    = plotfile_var_index(pf, "density")
    sc%h_comp       = plotfile_var_index(pf, "h")
    sc%temp_comp    = plotfile_var_index(pf, "tfromp") 
    sc%pi_comp      = plotfile_var_index(pf, "pi") 
    sc%p0_comp      = plotfile_var_index(pf, "p0") 
    sc%magv_comp    = plotfile_var_index(pf, "magvel") 
    sc%XHe4_comp    = plotfile_var_index(pf, "X(He4)")
    sc%XC12_comp    = plotfile_var_index(pf, "X(C12)")
    sc%XO16_comp    = plotfile_var_index(pf, "X(O16)")
    sc%XdotHe4_comp = plotfile_var_index(pf, "omegadot(He4)")
    sc%XdotC12_comp = plotfile_var_index(pf, "omegadot(C12)")
    sc%XdotO16_comp = plotfile_var_index(pf, "omegadot(O16)")
    sc%rhopert_comp = plotfile_var_index(pf, "rhopert")
    sc%tpert_comp   = plotfile_var_index(pf, "tpert") 
    sc%enuc_comp    = plotfile_var_index(pf, "enucdot") 
    sc%xvel_comp    = plotfile_var_index(pf, "x_vel")
    sc%yvel_comp    = plotfile_var_index(pf, "y_vel") 
    sc%zvel_comp    = plotfile_var_index(pf, "z_vel") 
    sc%s_comp       = plotfile_var_index(pf, "entropy")
    sc%spert_comp   = plotfile_var_index(pf, "entropypert")

    if ( sc%dens_comp < 0 .or. sc%temp_comp < 0 .or. sc%h_comp < 0 .or. &
         sc%pi_comp < 0 .or. sc%p0_comp < 0 .or. sc%magv_comp < 0 .or. &
         sc%XHe4_comp < 0 .or. sc%XC12_comp < 0 .or. sc%XO16_comp < 0 .or. &
         sc%XdotHe4_comp < 0 .or. sc%XdotC12_comp < 0 .or. sc%XdotO16_comp < 0 .or. &
         sc%rhopert_comp < 0 .or. sc%tpert_comp < 0 .or. &
         sc%enuc_comp < 0 .or. &
         sc%xvel_comp < 0 .or. sc%yvel_comp < 0 .or. sc%zvel_comp < 0 .or. &
         sc%s_comp < 0 .or. sc%spert_comp < 0) then
       call bl_error("ERROR: plotfile varaible(s) not defined")
    endif
  end subroutine init_comps

  !Initialize geometry
  subroutine init_geometry(pf, globals_only, fullstar, geo)
    implicit none

    !Args
    logical, intent(in) :: globals_only, fullstar
    type(plotfile), intent(in) :: pf

    type(geometry), intent(out) :: geo

    !Local
    real(kind=dp_t) :: x_maxdist, y_maxdist, z_maxdist, maxdist
    integer :: i, nbins

    ! get dx for the coarsest level.  
    geo%dx = plotfile_get_dx(pf, 1)

    ! get the index bounds for the finest level.  Note, lo and hi are
    ! ZERO based indicies
    geo%flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
    geo%fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

    if (.not. globals_only) then
       print *, 'Size of domain (zones): ', geo%fhi(1)-geo%flo(1)+1, geo%fhi(2)-geo%flo(2)+1, geo%fhi(3)-geo%flo(3)+1
    endif

    ! determine location of star's center
    if(fullstar) then
      ! full star simulation, center is at 
      ! the center of the computational domain
      geo%xctr = HALF*(pf%phi(1) - pf%plo(1))
      geo%yctr = HALF*(pf%phi(2) - pf%plo(2))
      geo%zctr = HALF*(pf%phi(3) - pf%plo(3))
    else
      ! we have an octant simulation, center is origin
      geo%xctr = ZERO
      geo%yctr = ZERO
      geo%zctr = ZERO
    endif

    if (.not. globals_only) then
       print *, 'Center of the star: ', geo%xctr, geo%yctr, geo%zctr
    endif

    ! compute the size of the radially-binned array -- we'll do it to
    ! the furtherest corner of the domain
    x_maxdist = max(abs(pf%phi(1) - geo%xctr), abs(pf%plo(1) - geo%xctr))
    y_maxdist = max(abs(pf%phi(2) - geo%yctr), abs(pf%plo(2) - geo%yctr))
    z_maxdist = max(abs(pf%phi(3) - geo%zctr), abs(pf%plo(3) - geo%zctr))
    maxdist = sqrt(x_maxdist**2 + y_maxdist**2 + z_maxdist**2)

    geo%dx_fine = minval(plotfile_get_dx(pf, pf%flevel))
    nbins = int(maxdist/geo%dx_fine)

    allocate(geo%r(0:nbins-1))
    do i = 0, nbins-1
       geo%r(i) = (dble(i) + HALF)*geo%dx_fine
    enddo
  end subroutine init_geometry

  subroutine init_averages(nbins, radav)
    implicit none

    !Args
    integer, intent(in) :: nbins
    type(radial_averages), intent(out) :: radav

    !store nbins
    radav%nbins = nbins

    ! ncount keeps track of how many zones were written into
    ! each bin
    allocate(radav%ncount(0:nbins-1))

    ! allocate storage for the data 
    allocate(   radav%dens_avg_bin(0:nbins-1))
    allocate(   radav%temp_avg_bin(0:nbins-1))
    allocate(   radav%pres_avg_bin(0:nbins-1))
    allocate(   radav%magv_avg_bin(0:nbins-1))
    allocate(   radav%enuc_avg_bin(0:nbins-1))
    allocate(   radav%XHe4_avg_bin(0:nbins-1))
    allocate(   radav%XC12_avg_bin(0:nbins-1))
    allocate(   radav%XO16_avg_bin(0:nbins-1))
    allocate(   radav%XdotHe4_avg_bin(0:nbins-1))
    allocate(   radav%XdotC12_avg_bin(0:nbins-1))
    allocate(   radav%XdotO16_avg_bin(0:nbins-1))
    allocate(   radav%dens_rms_bin(0:nbins-1))
    allocate(   radav%temp_rms_bin(0:nbins-1))
    allocate(   radav%pres_rms_bin(0:nbins-1))
    allocate(radav%entropy_avg_bin(0:nbins-1))
    allocate(radav%entropy_rms_bin(0:nbins-1))

    radav%ncount(:) = 0
    radav%dens_avg_bin(:) = ZERO
    radav%temp_avg_bin(:) = ZERO
    radav%pres_avg_bin(:) = ZERO
    radav%magv_avg_bin(:) = ZERO
    radav%enuc_avg_bin(:) = ZERO
    radav%XHe4_avg_bin(:) = ZERO
    radav%XC12_avg_bin(:) = ZERO
    radav%XO16_avg_bin(:) = ZERO
    radav%XdotHe4_avg_bin(:) = ZERO
    radav%XdotC12_avg_bin(:) = ZERO
    radav%XdotO16_avg_bin(:) = ZERO
    radav%dens_rms_bin(:) = ZERO
    radav%temp_rms_bin(:) = ZERO
    radav%pres_rms_bin(:) = ZERO
    radav%entropy_avg_bin(:) = ZERO
    radav%entropy_rms_bin(:) = ZERO
  end subroutine init_averages
  
  subroutine analyze(pf, geo, sc, thist, glb, radav, hh, hheap_frac_input)
    !Modules
    use eos_module, only: eos_input_rh, eos_input_re, eos
    use eos_type_module, only: eos_t
    use network, only: nspec, network_species_index
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    implicit none

    !Args
    type(plotfile), intent(inout) :: pf  !Must be inout because fab binding writes to pf
    type(geometry), intent(in) :: geo
    type(state_comps), intent(in) :: sc
    type(temp_hist), intent(out) :: thist

    type(globals), intent(out), optional :: glb
    type(radial_averages), intent(out), optional :: radav
    type(hheap), intent(out), optional :: hh
    real(kind=dp_t), intent(in), optional :: hheap_frac_input

    !Local data
    integer, parameter :: HHEAP_LEN = 100000
    real(kind=dp_t), parameter :: HHEAP_FRAC = 0.005
    type(box) :: cur_box
    type(eos_t) :: eos_state
    type(globals) :: glb_local
    type(hheap)   :: hh_local
    logical, allocatable :: imask(:,:,:), fmask(:,:,:)
    logical :: do_globals, do_averages, do_hotspots
    integer :: r1, rr, n, i, j, k, ii, jj, kk
    integer :: il, ir, jl, jr, kl, kr
    integer :: indx, b, hs_count
    integer, allocatable :: cell_count(:)
    real(kind=dp_t) :: xx, yy, zz, r_zone
    real(kind=dp_t) :: xlo, xhi, cur_xlo, cur_xhi
    real(kind=dp_t) :: ylo, yhi, cur_ylo, cur_yhi
    real(kind=dp_t) :: tfe, xmass(nspec)
    real(kind=dp_t) :: zlo, zhi, cur_zlo, cur_zhi, eps, lhheap_frac
    real(kind=dp_t), pointer :: p(:,:,:,:)

    !See which optional arguments are available
    do_globals = present(glb)
    do_averages = present(radav)
    do_hotspots = present(hh)

    !Set hottest percent of cells to track
    if( present(hheap_frac_input) ) then
      lhheap_frac = hheap_frac_input
    else
      lhheap_frac = HHEAP_FRAC
    endif

    !If no optionals given, nothing to do
    if(.not. any( (/do_globals, do_averages, do_hotspots/) )) then
      return
    end if

    !Initialize averaging data structures, size(geo%r) is number of radial bins
    if(do_averages) then
      call init_averages(size(geo%r), radav)
    end if

    !Initialize temperature histogram
    allocate(thist%ell(pf%flevel))
    allocate(thist%tavg(pf%flevel,thist%nbins))
    allocate(thist%tfeavg(pf%flevel,thist%nbins))
    thist%tavg = 0
    thist%tfeavg = 0

    ! imask will be set to false if we've already output the data.
    ! Note, imask is defined in terms of the finest level.  As we loop
    ! over levels, we will compare to the finest level index space to
    ! determine if we've already output here
    !
    ! NOTE: For large enough grids this mask can take up quite a bit of RAM.
    ! A 4 level 512^3 run's mask took up 10s of GB.
    !
    allocate(imask(geo%flo(1):geo%fhi(1),geo%flo(2):geo%fhi(2),geo%flo(3):geo%fhi(3)))
    imask(:,:,:) = .true.
    allocate(cell_count(pf%flevel))
    cell_count(:) = 0

    !Initialize hotspot list
    if(do_hotspots) then
      !Track the hottest (lhheap_frac)*(total number of valid cells) cells
      r1  = 1 !Refinement ratio of current:finest
      do n = pf%flevel, 1, -1 !Loop over all levels, finest to coarsest
        !rr = product(pf%refrat(1:n-1,1)) !Refinement ratio of current:coarsest

        print *, 'level ', n, ' of ', pf%flevel
        print *, '========================'
        !For each level, loop over all boxes in that level
        do i = 1, nboxes(pf, n)
          !Get bounds of current box
          cur_box = get_box(pf, n, i)
          il = lwb(cur_box,1)
          ir = upb(cur_box,1)
          jl = lwb(cur_box,2)
          jr = upb(cur_box,2)
          kl = lwb(cur_box,3)
          kr = upb(cur_box,3)

          !Count cells based on imask.  imask has a resolution of the finest
          !level, so we must divide by the cube (it's 3D) of the current ratio 
          !between this level and the finest to properly count cells.
          cell_count(n) = cell_count(n) +                 &
                          count(imask(il*r1:(ir+1)*r1-1,  &
                                      jl*r1:(jr+1)*r1-1,  &
                                      kl*r1:(kr+1)*r1-1)) &
                          /r1**3

          !We've counted these, so set to false
          imask(il*r1:(ir+1)*r1-1,  &
                jl*r1:(jr+1)*r1-1,  &
                kl*r1:(kr+1)*r1-1) = .false.

        enddo
        ! adjust r1 for the next level in the loop
        ! WARNING: As of now Maestro only has uniform grids in a given level, so
        !   we just pick dimension 1 for getting refrat, but this won't work for
        !   codes that can have different ratios along different dimensions
        if ( n /= 1 ) r1 = r1*pf%refrat(n-1,1)
        
        print *, 'Unique cells on level ', n , ': ', cell_count(n)
      enddo
     
      hs_count = int(lhheap_frac*sum(cell_count))
      call init_hheap(hh, hs_count)
      print *, 'Tracking ', hs_count, ' hotspots'

      !Reset imask
      imask(:,:,:) = .true.
    end if

    ! loop over the data, starting at the finest grid, and if we haven't
    ! already stored data in that grid location (according to imask),
    ! store it. 
    ! Calculation kernels are executed based on which arguments
    ! are passed in.

    ! r1 is the factor between the current level grid spacing and the
    ! FINEST level
    r1  = 1
    do i = pf%flevel, 1, -1
      ! rr is the factor between the COARSEST level grid spacing and
      ! the current level
      rr = product(pf%refrat(1:i-1,1))

      !If we're doing more than just globals to stdout, update the user on progress
      if ( any((/do_averages, do_hotspots/)) ) then
         print *, 'processing level ', i, ' rr = ', rr
      endif
      !$omp parallel default(none) firstprivate(i, r1, rr)                     &
      !$omp firstprivate(glb_local, hs_count)                                  &
      !$omp private(j,ii,jj,kk, cur_box, hh_local)                             &
      !$omp private(il, ir, jl, jr, kl, kr, p, xx, yy, zz)                     &
      !$omp shared(pf, geo, sc, imask, do_averages, do_globals, do_hotspots)   &
      !$omp shared(radav, glb, hh)   
      
      if (do_hotspots) then
         call init_hheap(hh_local, hs_count)
      endif

      !$omp do
      do j = 1, nboxes(pf, i)
         ! read in the data 1 patch at a time -- read in all the variables
         call fab_bind(pf, i, j)

         !Get bounds of current box
         cur_box = get_box(pf, i, j)
         il = lwb(cur_box,1)
         ir = upb(cur_box,1)
         jl = lwb(cur_box,2)
         jr = upb(cur_box,2)
         kl = lwb(cur_box,3)
         kr = upb(cur_box,3)

         ! get a pointer to the current patch
         p => dataptr(pf, i, j)

         ! loop over all of the zones in the patch.  Here, we convert
         ! the cell-centered indices at the current level into the
         ! corresponding RANGE on the finest level, and test if we've
         ! stored data in any of those locations.  If we haven't then
         ! we store this level's data and mark that range as filled.
         do kk = lbound(p,dim=3), ubound(p,dim=3)
            zz = (kk + HALF)*geo%dx(3)/rr + pf%plo(3)

            do jj = lbound(p,dim=2), ubound(p,dim=2)
               yy = (jj + HALF)*geo%dx(2)/rr + pf%plo(2)

               do ii = lbound(p,dim=1), ubound(p,dim=1)
                  xx = (ii + HALF)*geo%dx(1)/rr + pf%plo(1)

                  ! since a coarse zone will completely encompass this group
                  ! of fine zones, we can use any() here or all().  If any()
                  ! is true, then all() will be too.
                  if ( any(imask(ii*r1:(ii+1)*r1-1, &
                                 jj*r1:(jj+1)*r1-1, &
                                 kk*r1:(kk+1)*r1-1) ) ) then

                     if(do_averages) then
                        call averages_kernel(xx, yy, zz, ii, jj, kk, r1, p, geo, sc, radav)
                     end if

                     if(do_globals) then
                        !call globals_kernel(xx, yy, zz, ii, jj, kk, p, geo, sc, glb)
                        call globals_kernel(xx, yy, zz, ii, jj, kk, p, geo, sc, glb_local)
                     endif

                     if(do_hotspots) then
                        call hotspots_kernel(xx, yy, zz, ii, jj, kk, p, geo, sc, i, hh_local)
                     endif

                     imask(ii*r1:(ii+1)*r1-1, &
                           jj*r1:(jj+1)*r1-1, &
                           kk*r1:(kk+1)*r1-1) = .false.
                  end if
               end do
            enddo
         enddo

         call fab_unbind(pf, i, j)
      end do
      !$omp end do nowait

      if(do_globals) then
         !Set the shared, global glb object's variables
         !to be that of the thread with the highest global variables
         !$omp critical (glb_reduction)
         if(glb_local%T_peak > glb%T_peak) then
            glb%T_peak     = glb_local%T_peak     
            glb%xloc_Tpeak = glb_local%xloc_Tpeak 
            glb%yloc_Tpeak = glb_local%yloc_Tpeak 
            glb%zloc_Tpeak = glb_local%zloc_Tpeak 
            glb%R_Tpeak    = glb_local%R_Tpeak    
            glb%vx_Tpeak   = glb_local%vx_Tpeak   
            glb%vy_Tpeak   = glb_local%vy_Tpeak   
            glb%vz_Tpeak   = glb_local%vz_Tpeak   
            glb%vr_Tpeak   = glb_local%vr_Tpeak   
         endif
         if (glb_local%enuc_peak > glb%enuc_peak) then
            glb%enuc_peak     = glb_local%enuc_peak     
            glb%xloc_enucpeak = glb_local%xloc_enucpeak 
            glb%yloc_enucpeak = glb_local%yloc_enucpeak 
            glb%zloc_enucpeak = glb_local%zloc_enucpeak 
            glb%R_enucpeak    = glb_local%R_enucpeak    
         endif
         !$omp end critical (glb_reduction)
      endif

      if(do_hotspots) then
         !Merge the thread's local heap into the global, shared heap
         !$omp critical (hh_reduction)
         do k=1, size(hh_local%list)
            if(hh_local%list(k)%temp > hh%min_temp) then
               call hheap_add(hh, hh_local%list(k))
            endif
         enddo
         !$omp end critical (hh_reduction)
      endif

      !$omp end parallel

      ! adjust r1 for the next lowest level
      if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
    end do

    r1  = 1
    ! fmask is used to mark cells covered by the finest level of refinement.
    deallocate(imask)
    allocate(fmask(geo%flo(1):geo%fhi(1),geo%flo(2):geo%fhi(2),geo%flo(3):geo%fhi(3)))
    fmask(:,:,:) = .false.
    do i = pf%flevel, 1, -1
      ! rr is the factor between the COARSEST level grid spacing and
      ! the current level
      rr = product(pf%refrat(1:i-1,1))

      ! Set the lengthscale for this level
      !  ASSUMPTION: dx is the same for all dimensions
      thist%ell(i) = geo%dx(1)/rr
      
      !update the user on progress
      print *, 'thist processing level ', i, ' rr = ', rr
      
      !$omp parallel do default(none) firstprivate(i, r1, rr)                     &
      !$omp private(j,ii,jj,kk, b, cur_box, xmass, eos_state)                  &
      !$omp private(il, ir, jl, jr, kl, kr, p, xx, yy, zz)                     &
      !$omp shared(pf, geo, sc, fmask, thist)                     
      do j = 1, nboxes(pf, i)
         ! read in the data 1 patch at a time -- read in all the variables
         call fab_bind(pf, i, j)

         !Get bounds of current box
         cur_box = get_box(pf, i, j)
         il = lwb(cur_box,1)
         ir = upb(cur_box,1)
         jl = lwb(cur_box,2)
         jr = upb(cur_box,2)
         kl = lwb(cur_box,3)
         kr = upb(cur_box,3)

         !Mark the finest level
         if(i == pf%flevel) then
            fmask(il*r1:(ir+1)*r1-1,  &
                  jl*r1:(jr+1)*r1-1,  &
                  kl*r1:(kr+1)*r1-1) = .true.
         endif

         ! get a pointer to the current patch
         p => dataptr(pf, i, j)

         ! loop over all of the zones in the patch.  Here, we convert
         ! the cell-centered indices at the current level into the
         ! corresponding RANGE on the finest level, and test if we've
         ! stored data in any of those locations.  If we haven't then
         ! we store this level's data and mark that range as filled.
         do kk = lbound(p,dim=3), ubound(p,dim=3)
            zz = (kk + HALF)*geo%dx(3)/rr + pf%plo(3)

            do jj = lbound(p,dim=2), ubound(p,dim=2)
               yy = (jj + HALF)*geo%dx(2)/rr + pf%plo(2)

               do ii = lbound(p,dim=1), ubound(p,dim=1)
                  xx = (ii + HALF)*geo%dx(1)/rr + pf%plo(1)

                  !The finest level tracks the regions of vigorous burning, so
                  !we use it to mark where we want to calculate turbulent
                  !fluctuations at varying lengthscales.
                  if ( any(fmask(ii*r1:(ii+1)*r1-1, &
                                 jj*r1:(jj+1)*r1-1, &
                                 kk*r1:(kk+1)*r1-1) ) ) then
                     b = min( &
                           nint((p(ii,jj,kk,sc%temp_comp) - thist%min_temp)/thist%dT), &
                           thist%nbins - 1 &
                         )
                     !$omp atomic
                     thist%tavg(i,b) = thist%tavg(i,b) + 1

                     xmass(network_species_index('He4')) = p(ii,jj,kk,sc%XHe4_comp)
                     xmass(network_species_index('C12')) = p(ii,jj,kk,sc%XC12_comp)
                     xmass(network_species_index('O16')) = p(ii,jj,kk,sc%XO16_comp)
                     xmass(network_species_index('Fe56')) = 0.0_dp_t

                     eos_state%rho   = p(ii,jj,kk,sc%dens_comp)
                     eos_state%h     = p(ii,jj,kk,sc%h_comp)
                     eos_state%T     = p(ii,jj,kk,sc%temp_comp) !initial guess
                     eos_state%xn(:) = xmass(:)
                       
                     call eos(eos_input_rh, eos_state, .false.)
                     
                     b = min( &
                           nint((eos_state%T - thist%min_temp)/thist%dT), &
                           thist%nbins - 1 &
                         )
                     b = max(1,b)
                     !$omp atomic
                     thist%tfeavg(i,b) = thist%tfeavg(i,b) + 1
                  endif
               end do
            enddo
         enddo

         call fab_unbind(pf, i, j)
      end do
      !$omp end parallel do 

      ! adjust r1 for the next lowest level
      if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
    end do

    ! normalize
    if(do_averages) then
      do i = 0, radav%nbins-1
         if (radav%ncount(i) /= 0) then
            ! simple averages
            radav%dens_avg_bin(i) = radav%dens_avg_bin(i)/radav%ncount(i)
            radav%temp_avg_bin(i) = radav%temp_avg_bin(i)/radav%ncount(i)
            radav%pres_avg_bin(i) = radav%pres_avg_bin(i)/radav%ncount(i)
            radav%magv_avg_bin(i) = radav%magv_avg_bin(i)/radav%ncount(i)
            radav%enuc_avg_bin(i) = radav%enuc_avg_bin(i)/radav%ncount(i)
            radav%entropy_avg_bin(i) = radav%entropy_avg_bin(i)/radav%ncount(i)

            ! Favre averaged composition
            radav%XHe4_avg_bin(i) = (radav%XHe4_avg_bin(i)/radav%ncount(i)) / radav%dens_avg_bin(i)
            radav%XC12_avg_bin(i) = (radav%XC12_avg_bin(i)/radav%ncount(i)) / radav%dens_avg_bin(i)
            radav%XO16_avg_bin(i) = (radav%XO16_avg_bin(i)/radav%ncount(i)) / radav%dens_avg_bin(i)
            radav%XdotHe4_avg_bin(i) = (radav%XdotHe4_avg_bin(i)/radav%ncount(i)) / radav%dens_avg_bin(i)
            radav%XdotC12_avg_bin(i) = (radav%XdotC12_avg_bin(i)/radav%ncount(i)) / radav%dens_avg_bin(i)
            radav%XdotO16_avg_bin(i) = (radav%XdotO16_avg_bin(i)/radav%ncount(i)) / radav%dens_avg_bin(i)

            ! RMS quantities
            radav%dens_rms_bin(i) = sqrt(radav%dens_rms_bin(i)/radav%ncount(i))        
            radav%temp_rms_bin(i) = sqrt(radav%temp_rms_bin(i)/radav%ncount(i))
            radav%pres_rms_bin(i) = sqrt(radav%pres_rms_bin(i)/radav%ncount(i))
            radav%entropy_rms_bin(i) = sqrt(radav%entropy_rms_bin(i)/radav%ncount(i))
         endif
      enddo
    end if

  end subroutine analyze

  !An innermost loop kernel for calculating averages
  !
  ! Here we compute angle-averaged quantities, <q>, and RMS quantities,
  ! q' = sqrt { < (q - <q>)**2 > }, where the averaging is done at constant
  ! radius.
  !
  ! For density, rho_0 = <rho>, and in the plotfiles we store rhopert =
  ! rho - rho0, so we can compute the RMS density directly from this.
  !
  ! Similarly, for temperature, and entropy (using tpert and spert from
  ! the plotfiles).!
  !
  !Assumptions: 
  !  +the radial_averages object is initialized approriately
  !  before the first call to this routine
  !  +pointer p is properly initialized, bound to a plotfile patch
  subroutine averages_kernel(xx, yy, zz, ii, jj, kk, r1, p, geo, sc, radav)
    !Modules
    implicit none

    !Args
    real(kind=dp_t), intent(in) :: xx, yy, zz
    integer, intent(in) :: ii, jj, kk, r1
    real(kind=dp_t), pointer, intent(inout) :: p(:,:,:,:)
    type(geometry), intent(in) :: geo
    type(state_comps), intent(in) :: sc
    type(radial_averages), intent(inout) :: radav

    !Local
    real(kind=dp_t) :: r_zone
    integer :: indx

    !Determine bin index
    r_zone = sqrt((xx-geo%xctr)**2 + (yy-geo%yctr)**2 + (zz-geo%zctr)**2)
    indx = r_zone/geo%dx_fine

    ! weight the zone's data by its size
    ! note, for p(:,:,:,n), n refers to index of the
    ! variable as found via plotfile_var_index
    !$omp atomic
    radav%dens_avg_bin(indx) = radav%dens_avg_bin(indx) + &
         p(ii,jj,kk,sc%dens_comp) * r1**3

    !$omp atomic
    radav%temp_avg_bin(indx) = radav%temp_avg_bin(indx) + &
         p(ii,jj,kk,sc%temp_comp) * r1**3
    
    !$omp atomic
    radav%pres_avg_bin(indx) = radav%pres_avg_bin(indx) + & !P = P_0 + \pi
         (p(ii,jj,kk,sc%p0_comp) + p(ii,jj,kk,sc%pi_comp)) * r1**3
    
    !$omp atomic
    radav%magv_avg_bin(indx) = radav%magv_avg_bin(indx) + &
         p(ii,jj,kk,sc%magv_comp) * r1**3
    
    !$omp atomic
    radav%enuc_avg_bin(indx) = radav%enuc_avg_bin(indx) + &
         p(ii,jj,kk,sc%enuc_comp) * r1**3
    
    !$omp atomic
    radav%entropy_avg_bin(indx) = radav%entropy_avg_bin(indx) + &
         p(ii,jj,kk,sc%s_comp) * r1**3

    ! do the Favre-average here, < rho * X(He4) > / < rho >
    !$omp atomic
    radav%XHe4_avg_bin(indx) = radav%XHe4_avg_bin(indx) + &
         p(ii,jj,kk,sc%dens_comp)*p(ii,jj,kk,sc%XHe4_comp) * r1**3

    ! do the Favre-average here, < rho * X(C12) > / < rho >
    !$omp atomic
    radav%XC12_avg_bin(indx) = radav%XC12_avg_bin(indx) + &
         p(ii,jj,kk,sc%dens_comp)*p(ii,jj,kk,sc%XC12_comp) * r1**3

    ! do the Favre-average here, < rho * X(O16) > / < rho >
    !$omp atomic
    radav%XO16_avg_bin(indx) = radav%XO16_avg_bin(indx) + &
         p(ii,jj,kk,sc%dens_comp)*p(ii,jj,kk,sc%XO16_comp) * r1**3
    
    ! do the Favre-average here, < rho * Xdot(He4) > / < rho >
    !$omp atomic
    radav%XdotHe4_avg_bin(indx) = radav%XdotHe4_avg_bin(indx) + &
         p(ii,jj,kk,sc%dens_comp)*p(ii,jj,kk,sc%XdotHe4_comp) * r1**3

    ! do the Favre-average here, < rho * Xdot(C12) > / < rho >
    !$omp atomic
    radav%XdotC12_avg_bin(indx) = radav%XdotC12_avg_bin(indx) + &
         p(ii,jj,kk,sc%dens_comp)*p(ii,jj,kk,sc%XdotC12_comp) * r1**3

    ! do the Favre-average here, < rho * Xdot(O16) > / < rho >
    !$omp atomic
    radav%XdotO16_avg_bin(indx) = radav%XdotO16_avg_bin(indx) + &
         p(ii,jj,kk,sc%dens_comp)*p(ii,jj,kk,sc%XdotO16_comp) * r1**3

    ! for the RMS quantities, we use the perturbational quantities
    ! already stored in the plotfile
    !$omp atomic
    radav%dens_rms_bin(indx) = radav%dens_rms_bin(indx) + &
         p(ii,jj,kk,sc%rhopert_comp)*p(ii,jj,kk,sc%rhopert_comp) * r1**3

    !$omp atomic
    radav%temp_rms_bin(indx) = radav%temp_rms_bin(indx) + &
         p(ii,jj,kk,sc%tpert_comp)*p(ii,jj,kk,sc%tpert_comp) * r1**3

    !$omp atomic
    radav%pres_rms_bin(indx) = radav%pres_rms_bin(indx) + &  !\pi = (p - <p>)
         p(ii,jj,kk,sc%pi_comp)*p(ii,jj,kk,sc%pi_comp) * r1**3

    !$omp atomic
    radav%entropy_rms_bin(indx) = radav%entropy_rms_bin(indx) + &
         p(ii,jj,kk,sc%spert_comp)*p(ii,jj,kk,sc%spert_comp) * r1**3

    !$omp atomic
    radav%ncount(indx) = radav%ncount(indx) + r1**3
  end subroutine averages_kernel

  !An innermost loop kernel for calculating global values
  !Assumptions: 
  !   +the globals object is initialized appropriately before the first call to this routine
  !   +p is properly initialized, bound to a plotfile patch
  subroutine globals_kernel(xx, yy, zz, ii, jj, kk, p, geo, sc, glb)
    !Modules
    implicit none

    !Args
    real(kind=dp_t), intent(in) :: xx, yy, zz
    integer, intent(in) :: ii, jj, kk
    real(kind=dp_t), pointer, intent(in) :: p(:,:,:,:)
    type(geometry), intent(in) :: geo
    type(state_comps), intent(in) :: sc
    type(globals), intent(inout) :: glb

    !Local

    ! store the location and value of the peak temperature
    !!$omp critical (crit_tpeak)
    if (p(ii,jj,kk,sc%temp_comp) > glb%T_peak) then
       glb%T_peak = p(ii,jj,kk,sc%temp_comp)
       glb%xloc_Tpeak = xx
       glb%yloc_Tpeak = yy
       glb%zloc_Tpeak = zz

       glb%R_Tpeak = sqrt( (geo%xctr - glb%xloc_Tpeak)**2 + &
                       (geo%yctr - glb%yloc_Tpeak)**2 + &
                       (geo%zctr - glb%zloc_Tpeak)**2 )

       glb%vx_Tpeak = p(ii,jj,kk,sc%xvel_comp)
       glb%vy_Tpeak = p(ii,jj,kk,sc%yvel_comp)
       glb%vz_Tpeak = p(ii,jj,kk,sc%zvel_comp)

       ! this is (v . e_r)
       glb%vr_Tpeak = glb%vx_Tpeak*(xx - geo%xctr)/glb%R_Tpeak + &
                  glb%vy_Tpeak*(yy - geo%yctr)/glb%R_Tpeak + &
                  glb%vz_Tpeak*(zz - geo%zctr)/glb%R_Tpeak
    endif
    !!$omp end critical (crit_tpeak)

    ! store the location and value of the peak enucdot
    !!$omp critical (crit_enuc)
    if (p(ii,jj,kk,sc%enuc_comp) > glb%enuc_peak) then
       glb%enuc_peak = p(ii,jj,kk,sc%enuc_comp)

       glb%xloc_enucpeak = xx
       glb%yloc_enucpeak = yy
       glb%zloc_enucpeak = zz

       glb%R_enucpeak = sqrt( (geo%xctr - glb%xloc_enucpeak)**2 + &
                          (geo%yctr - glb%yloc_enucpeak)**2 + &
                          (geo%zctr - glb%zloc_enucpeak)**2 )
    endif
    !!$omp end critical (crit_enuc)
  end subroutine globals_kernel

  !An innermost loop kernel for calculating the hottest cells
  !Assumptions: 
  !   +the hheap object is initialized appropriately before the first call to this routine
  !   +p is properly initialized, bound to a plotfile patch
  subroutine hotspots_kernel(xx, yy, zz, ii, jj, kk, p, geo, sc, lev, hh)
    !Modules
    implicit none

    !Args
    real(kind=dp_t), intent(in) :: xx, yy, zz
    integer, intent(in) :: ii, jj, kk
    real(kind=dp_t), pointer, intent(in) :: p(:,:,:,:)
    type(geometry), intent(in) :: geo
    type(state_comps), intent(in) :: sc
    integer, intent(in) :: lev
    type(hheap), intent(inout) :: hh

    !Local
    type(hotspot) :: newhs
    integer :: i

    !If temp is more than hh's minimum, add this cell to the heap
    !!$omp critical (crit_hheap)
    if(p(ii,jj,kk,sc%temp_comp) > hh%min_temp) then
      newhs%temp = p(ii,jj,kk,sc%temp_comp)
      newhs%rho = p(ii,jj,kk,sc%dens_comp)
      newhs%x = xx
      newhs%y = yy
      newhs%z = zz
      newhs%level = lev
    
      call hheap_add(hh, newhs)
    end if
    !!$omp end critical (crit_hheap)
  end subroutine hotspots_kernel

  !An innermost loop kernel for calculating the temperatures at
  !the lenghthscale of each level to compare with fluctuation expectations 
  !from Kolmogorov theory.
  !Assumptions: 
  !   +
  subroutine fluct_kernel(xx, yy, zz, ii, jj, kk, p, geo, sc, lev, hh)
    !Modules
    implicit none

    !Args
    real(kind=dp_t), intent(in) :: xx, yy, zz
    integer, intent(in) :: ii, jj, kk
    real(kind=dp_t), pointer, intent(in) :: p(:,:,:,:)
    type(geometry), intent(in) :: geo
    type(state_comps), intent(in) :: sc
    integer, intent(in) :: lev
    type(hheap), intent(inout) :: hh

    !Local
    type(hotspot) :: newhs
    integer :: i

    !!If temp is more than hh's minimum, add this cell to the list
    !if(p(ii,jj,kk,sc%temp_comp) > hh%min_temp) then
    !  newhs%temp = p(ii,jj,kk,sc%temp_comp)
    !  newhs%rho = p(ii,jj,kk,sc%dens_comp)
    !  newhs%x = xx
    !  newhs%y = yy
    !  newhs%z = zz
    !  newhs%level = lev
    !
    !  call hheap_add(hh, newhs)
    !end if
    !!Generate file with (l, <T>)
    !ell = 
    !T_avg = 
  end subroutine fluct_kernel


  !Write output for all arguments passed
  subroutine writeout(pf, slicefile, geo, sc, thist, glb, radav, hh)
    !Modules
    implicit none

    !Args
    type(plotfile), intent(in) :: pf 
    character(len=256), intent(in) :: slicefile
    type(geometry), intent(in) :: geo
    type(state_comps), intent(in) :: sc
    type(temp_hist), intent(in) :: thist
    type(radial_averages), intent(inout), optional :: radav
    type(globals), intent(in), optional :: glb
    type(hheap), intent(inout), optional :: hh

    !Local data
    !Note: not sure if the below string continuation technique is portable.  It's
    !in the f90 standard so I imagine it's mostly implemented.
    character(len=*), parameter :: fmt_header = '("# time = ", g24.12)' 
    character(len=*), parameter :: fmt_section = &
      '("# ---------------------------------------------------------------------------")'
    character(len=*), parameter :: fmt_globals = &
      '("# peak temperature - ", g24.12, /,&
      &"# peak temp loc (x,y,z) = ", 3(g24.12,1x), /,&
      &"# peak temp radius = ", g24.12, /,&
      &"# velocity @ peak T loc (vx, vy, vz) = ", 3(g24.12,1x), /,&
      &"# radial velocity @ peak T loc = ", g24.12,1x, /,&
      &"# peak enucdot = ", g24.12, /,&
      &"# peak enucdot loc (x,y,z) = ", 3(g24.12,1x), /,&
      &"# peak enucdot radius = ", g24.12)'
    character(len=*), parameter :: fmt_labels = '("#",100(a24,1x))'
    character(len=*), parameter :: fmt_data = '(1x,100(g24.12,1x))'
    character(len=*), parameter :: fmt_hdata = '(1x,5(g24.12,1x),I24)'
    integer, parameter :: EXT_LEN = 6 !Length of '.slice'

    character(len=256) :: hotfile, histfile
    logical :: write_globals, write_averages, write_hotspots
    integer :: uno, i, j
    real(kind=dp_t) :: cur_temp

    !See which args passed
    write_globals = present(glb)
    write_averages = present(radav)
    write_hotspots = present(hh)

    ! slicefile
    if (any((/ write_averages/))) then
      uno = unit_new()
      open(unit=uno, file=slicefile, status = 'replace')

      ! write the header
      write(uno,fmt_header) pf%tm
      write(uno,fmt_section)

      ! write globals
      write(uno,fmt_globals) glb%T_peak, glb%xloc_Tpeak, glb%yloc_Tpeak, glb%zloc_Tpeak, &
        glb%R_Tpeak, glb%vx_Tpeak, glb%vy_Tpeak, glb%vz_Tpeak, glb%vr_Tpeak, glb%enuc_peak, &
        glb%xloc_enucpeak, glb%yloc_enucpeak, glb%zloc_enucpeak, glb%R_enucpeak
      write(uno,fmt_section)

      ! write averages labels
      write(uno,fmt_labels) "r", "density", "temperature", "pressure",  &
           "vel. magnitude", "Hnuc", "entropy", "X(He4)", "X(C12)", "X(O16)", &
           "Xdot(He4)", "Xdot(C12)", "Xdot(O16)", "RMS density", "RMS temperature", &
           "RMS pressure", "RMS entropy"

      ! write the data in columns
      do i = 0, radav%nbins-1
         ! Use this to protect against a number being xx.e-100
         !   which will print without the "e"
         if (abs(   radav%dens_avg_bin(i)) .lt. 1.d-99)    radav%dens_avg_bin(i) = 0.d0
         if (abs(   radav%temp_avg_bin(i)) .lt. 1.d-99)    radav%temp_avg_bin(i) = 0.d0
         if (abs(   radav%pres_avg_bin(i)) .lt. 1.d-99)    radav%pres_avg_bin(i) = 0.d0
         if (abs(   radav%magv_avg_bin(i)) .lt. 1.d-99)    radav%magv_avg_bin(i) = 0.d0
         if (abs(   radav%enuc_avg_bin(i)) .lt. 1.d-99)    radav%enuc_avg_bin(i) = 0.d0
         if (abs(radav%entropy_avg_bin(i)) .lt. 1.d-99) radav%entropy_avg_bin(i) = 0.d0
         if (abs(   radav%XHe4_avg_bin(i)) .lt. 1.d-99)    radav%XHe4_avg_bin(i) = 0.d0
         if (abs(   radav%XC12_avg_bin(i)) .lt. 1.d-99)    radav%XC12_avg_bin(i) = 0.d0
         if (abs(   radav%XO16_avg_bin(i)) .lt. 1.d-99)    radav%XO16_avg_bin(i) = 0.d0
         if (abs(   radav%XdotHe4_avg_bin(i)) .lt. 1.d-99) radav%XdotHe4_avg_bin(i) = 0.d0
         if (abs(   radav%XdotC12_avg_bin(i)) .lt. 1.d-99) radav%XdotC12_avg_bin(i) = 0.d0
         if (abs(   radav%XdotO16_avg_bin(i)) .lt. 1.d-99) radav%XdotO16_avg_bin(i) = 0.d0
         if (abs(   radav%dens_rms_bin(i)) .lt. 1.d-99)    radav%dens_rms_bin(i) = 0.d0
         if (abs(   radav%temp_rms_bin(i)) .lt. 1.d-99)    radav%temp_rms_bin(i) = 0.d0
         if (abs(   radav%pres_rms_bin(i)) .lt. 1.d-99)    radav%pres_rms_bin(i) = 0.d0
         if (abs(radav%entropy_rms_bin(i)) .lt. 1.d-99) radav%entropy_rms_bin(i) = 0.d0
         
         write(uno,fmt_data) geo%r(i), radav%dens_avg_bin(i), radav%temp_avg_bin(i),        &
              radav%pres_avg_bin(i), radav%magv_avg_bin(i), radav%enuc_avg_bin(i),          &
              radav%entropy_avg_bin(i),                                                     &
              radav%XHe4_avg_bin(i), radav%XC12_avg_bin(i), radav%XO16_avg_bin(i),          &
              radav%XdotHe4_avg_bin(i), radav%XdotC12_avg_bin(i), radav%XdotO16_avg_bin(i), &
              radav%dens_rms_bin(i), radav%temp_rms_bin(i),                                 &
              radav%pres_rms_bin(i), radav%entropy_rms_bin(i)
      end do

      close(unit=uno)
    endif

    ! hotspots file
    if (write_hotspots) then
      !build filename
      hotfile = trim(slicefile)
      hotfile = hotfile(1:len(trim(hotfile)) - EXT_LEN) // ".hotspots"
     
      !Get unit, open file 
      uno = unit_new()
      open(unit=uno, file=hotfile, status = 'replace')

      !Write header
      write(uno,fmt_header) pf%tm
      write(uno,fmt_labels) "temperature", "density", "x", "y", "z", "level"
      write(uno,fmt_section)

      !Write data in columns
      do i = 1, hh%length
         ! Use this to protect against a number being xx.e-100
         !   which will print without the "e"
         if (abs(hh%list(i)%temp) .lt. 1.d-99) hh%list(i)%temp = 0.d0
         if (abs(hh%list(i)%rho)  .lt. 1.d-99) hh%list(i)%rho  = 0.d0
         if (abs(hh%list(i)%x)    .lt. 1.d-99) hh%list(i)%x    = 0.d0
         if (abs(hh%list(i)%y)    .lt. 1.d-99) hh%list(i)%y    = 0.d0
         if (abs(hh%list(i)%z)    .lt. 1.d-99) hh%list(i)%z    = 0.d0
         
         write(uno,fmt_hdata) hh%list(i)%temp, hh%list(i)%rho, & 
              hh%list(i)%x, hh%list(i)%y, hh%list(i)%z, hh%list(i)%level
      end do


      close(unit=uno)
    endif

    ! temp histogram file
    !build filename
    histfile = trim(slicefile)
    histfile = histfile(1:len(trim(histfile)) - EXT_LEN) // ".temphist"

    !Get unit, open file
    uno = unit_new()
    open(unit=uno, file=histfile, status = 'replace')

    !Write header
    write(uno,fmt_header) pf%tm
    write(uno,fmt_labels) "lengthscale"
    write(uno,fmt_labels) "  temperature", "cell T counts", "T from EoS counts"
    write(uno,fmt_section)

    !Write data
    do i = pf%flevel, 1, -1
       write(uno,fmt_data) thist%ell(i)
       do j = 1, thist%nbins
          cur_temp = thist%min_temp + (j-1) * thist%dT
          write(uno,fmt_data) cur_temp, thist%tavg(i,j), thist%tfeavg(i,j)
       enddo
    enddo
    close(unit=uno)

    !stdout output
    if (any((/ write_averages/))) then
      print *, 'Peak temperature = ', glb%T_peak
      print *, 'Peak temperature location (x,y,z) = ', glb%xloc_Tpeak, glb%yloc_Tpeak, glb%zloc_Tpeak
      print *, 'Peak temperature radius = ', glb%R_Tpeak
      print *, 'Peak enucdot = ', glb%enuc_peak
      print *, 'Peak enucdot location (x,y,z) = ', glb%xloc_enucpeak, glb%yloc_enucpeak, glb%zloc_enucpeak
      print *, 'Peak enucdot radius = ', glb%R_enucpeak
    !Else, user specified --globals-only, so use that format
    else
      write (*,fmt_labels) "time", "T_peak", "x(T_peak)", "y(T_peak)", "z(T_peak)", "R(T_peak)", &
           "enuc_peak", "x(enuc_peak)", "y(enuc_peak)", "z(enuc_peak)", "R(enuc_peak)"
      write (*,fmt_data) pf%tm, glb%T_peak, glb%xloc_Tpeak, glb%yloc_Tpeak, glb%zloc_Tpeak, glb%R_Tpeak, &
           glb%enuc_peak, glb%xloc_enucpeak, glb%yloc_enucpeak, glb%zloc_enucpeak, glb%R_enucpeak
    endif
  end subroutine writeout
end module subchandra
