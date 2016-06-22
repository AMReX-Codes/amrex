!
! this program generates average and rms quantities as a function of "height"
! along a specified direction for a specified list of variables
! (defaults to all variables) in a specified plotfile.
!
! Calling sequence:
!    faverage -p <pltfile> [-o <outputfile> -d <dir> -v <num_vars> <var_list> -q]
!
! Flags:
!     -p|--pltfile <filename>    : specifies the input plotfile 
!
!     -o|--outputfile <filename> : specifies the output file
!                                  default is to stdout
!     
!     -d|--dir <integer>         : specifies the "height" direction
!                                  default is 2 (y) for 2-d plotfiles and 
!                                  3 (z) for 3-d plotfiles
!
!     -v|--vars <integer> <list> : specifies a specific <list> of variables to
!                                  analyze.  <integer> is the number of 
!                                  variables in the whitespace separated <list>
!                                  of variable names.
!                                  default is to analyze all variables
!
!     -f|--favre                 : toggles favre (density-weighted) averaging
!                                  default is to do normal averaging
!
!     -q|--quiet                 : perform averaging operations quietly
!

program faverage

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module

  implicit none

  ! argument variables
  character(len=256) :: pltfile, outputfile
  integer :: idir, nvars_in
  logical :: idir_specified, do_favre
  character(len=80), allocatable :: vars_in(:)
  logical :: quiet
  
  ! auxilliary variables for arguments
  logical :: do_specific_vars
  integer, allocatable :: icomps_in(:)
  logical, allocatable :: skip(:)
  character :: dir

  ! fk2cli variables
  integer :: narg, farg
  character(len=256) :: fname

  ! local variables
  type(plotfile) :: pf
  integer :: unit, uout
  integer :: dim, irho_comp, irho_comp_pass
  integer :: i, j
  real(kind=dp_t) :: dx(MAX_SPACEDIM)
  integer, dimension(MAX_SPACEDIM) :: flo, fhi, lo, hi
  logical, allocatable :: imask(:,:,:)
  real(kind=dp_t), allocatable :: avg(:,:), rms(:,:), cnt(:)
  real(kind=dp_t), pointer :: p(:,:,:,:)

  ! variables to be passed into averaging routines
  ! nvars_pass contains the total number of valid variables to be analyzed
  ! icomps_pass contains the indices (within the pf) of the valid variables
  !     being analyzed
  integer :: nvars_pass
  integer, allocatable :: icomps_pass(:)

  integer :: r1, max_points, ii, jj, kk, index, index_lo, index_hi
  real(kind=dp_t) :: weight

  ! make formatting flexible
  character(len=128) :: column_format, data_format, columns


  unit = unit_new()

  ! defaults
  pltfile = ''
  outputfile = 'stdout'
  idir_specified = .false.
  idir = 1
  dir = "x"
  do_favre = .false.
  do_specific_vars = .false.
  irho_comp = -10
  quiet = .false.

  !============================================================
  ! parse command line arguments                              |
  !============================================================
  narg = command_argument_count()

  farg = 1
  do while (farg <= narg)
     call get_command_argument(farg, value = fname)

     select case(fname)

     case('-p', '--pltfile')
        farg = farg + 1
        call get_command_argument(farg, value = pltfile)

     case('-o', '--outputfile')
        farg = farg + 1
        call get_command_argument(farg, value = outputfile)

     case('-d', '--dir')
        idir_specified = .true.
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) idir
        if (idir < 1 .or. idir > MAX_SPACEDIM) &
             call bl_error("-d option set with invalid direction")

     case('-v', '--vars')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) nvars_in

        ! allocate storage for the variables that are specified
        if (nvars_in > 0) then
           allocate(   vars_in(nvars_in), &
                     icomps_in(nvars_in), &
                          skip(nvars_in))
           do_specific_vars = .true.

           do i = 1, nvars_in
              farg = farg + 1
              call get_command_argument(farg, value = fname)
              vars_in(i) = trim(fname)
           enddo
        else
           call bl_error("-v option set with negative number of variables")
        endif

     case('-f', '--favre')
        do_favre = .true.

     case('-q', '--quiet')
        quiet = .true.

     case default
        exit

     end select
     farg = farg + 1
  enddo

  !============================================================
  ! sanity checks on specified options                        |
  !============================================================
  if (pltfile == '') then
     call print_usage()
     stop
  endif

  ! build the plotfile datatype
  call build(pf, pltfile, unit)

  dim = pf%dim

  ! make sure we don't specify z-direction for a 2d pltfile
  if (idir > dim) call bl_error("error: idir too large for pltfile")

  ! if the idir was not specified with the -d flag, then use the direction
  ! perpendicular to gravity
  if (.not. idir_specified) then
     if (dim == 2) idir = 2
     if (dim == 3) idir = 3
  endif

  if (idir == 2) dir = "y"
  if (idir == 3) dir = "z"

  ! make sure there weren't too many variables specified
  if (do_specific_vars) then
     if (nvars_in > pf%nvars) &
          call bl_error("error: -v option set with too many vars for pltfile")
  endif

  ! grab the density index in the plotfile
  ! if it is not in the plotfile and we are doing favre averaging then we abort
  do i = 1, pf%nvars
     if (var_name(pf, i) == "density") then
        irho_comp = i
        exit
     endif
  enddo

  if (do_favre .and. irho_comp < 0) &
       call bl_error("-f option set but density is not in pltfile")

  !============================================================
  ! build list of variables to analyze                        |
  !============================================================
  if (do_specific_vars) then

     ! get the icomps of the specified variables and find which ones are not
     ! in the pltfile - we skip those
     icomps_in = -1
     skip = .false.

     do i = 1, nvars_in
        do j = 1, pf%nvars

           if (vars_in(i) == var_name(pf, j)) then
              icomps_in(i) = j
              exit
           endif
        enddo

        if (icomps_in(i) < 0) then
           skip(i) = .true.
           print *, 'WARNING: skipping variable not found in pltfile:', &
                    vars_in(i)
        endif
     enddo

     ! if none of the specified variables are found in the pltfile then we 
     ! abort
     if (all(skip)) &
          call bl_error("none of the specified vars were found in the pltfile")

     ! if we are favre averaging and density was not specified then we need
     ! to add it to the list
     if (do_favre .and. (all(icomps_in(:) /= irho_comp))) then

        nvars_pass = count(.not. skip) + 1

        allocate(icomps_pass(nvars_pass))

        ! keep track of the index of density in the pass
        irho_comp_pass = nvars_pass
        
        ! fill the icomps_pass array - this is a map between the indices of 
        !    the specified <list> and the indices within the pf of the 
        !    valid variables
        ! first set it equal to the index for density
        icomps_pass = irho_comp
        ! now use the (.not. skip) array as a mask for icomps to fill the 
        !    icomps_pass array with valid indices.
        icomps_pass = pack(icomps_in, (.not. skip), icomps_pass)

      ! otherwise just set up the variables we are not skipping
     else
        
        nvars_pass = count(.not. skip)

        allocate(icomps_pass(nvars_pass))

        icomps_pass = pack(icomps_in, (.not. skip))

      ! keep track of the index of density in the pass if favre averaging
        if (do_favre) then

           do i = 1, nvars_pass
              if (icomps_pass(i) == irho_comp) then
                 irho_comp_pass = i
                 exit
              endif
           enddo

        endif

     endif

  ! a list of variables was not specified
  ! use all the variables in the plotfile
  else

     nvars_pass = pf%nvars

     allocate(icomps_pass(nvars_pass))

     do i = 1, nvars_pass
        icomps_pass(i) = i
     enddo

     ! keep track of the density in the pass if favre averaging
     if (do_favre) irho_comp_pass = irho_comp

  endif

  !============================================================
  ! Allocate storage for the averaging                        |
  !============================================================
  ! grab the finest level's dx
  dx = plotfile_get_dx(pf, pf%flevel)

  ! get the index bounds for the finest level
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop 
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here
  if (dim == 2) then
     allocate(imask(flo(1):fhi(1), &
                    flo(2):fhi(2), &
                    1              ))
  else if (dim == 3) then
     allocate(imask(flo(1):fhi(1), &
                    flo(2):fhi(2), &
                    flo(3):fhi(3) ))
  endif

  ! this is the maximum number of points along the specified direction if
  ! the entire domain were at the finest level resolution
  max_points = fhi(idir) - flo(idir) + 1

  ! allocate storage for the data
  allocate(avg(0:max_points-1,nvars_pass), &
           rms(0:max_points-1,nvars_pass), &
           cnt(0:max_points-1))

  avg = ZERO
  rms = ZERO

  !============================================================
  ! Average the quantities                                    |
  !============================================================
  ! dump some info for the user
  if (.not. quiet) then
     print *, 'averaging...'
     print *, 'pltfile = "', trim(pltfile), '"'
     print *, 'outputfile = "', trim(outputfile), '"'
  endif

  ! r1 is the refinement factor between the current level and the FINEST level
  r1 = 1

  imask = .true.
  cnt = 0

  ! loop over the data starting at the finest grid, and if we havn't already
  ! stored the data in that grid location (according to imask), store it.
  ! Since we are doing averaging and rms calculations, we don't care about the 
  ! order in which data is stored.
  do i = pf%flevel, 1, -1

     do j = 1, nboxes(pf, i)

        ! bind to a comp vector
        call fab_bind_comp_vec(pf, i, j, icomps_pass)

        ! get the data
        p => dataptr(pf, i, j)
        lo(:) = 1
        hi(:) = 1
        lo(1:dim) = lwb(get_box(pf, i, j))
        hi(1:dim) = upb(get_box(pf, i, j))

        ! loop over all the zones in the current box
        ! Here, we convert the cell-centered indices at the current level
        ! into the corresponding RANGE on the finest level, and test if 
        ! we've stored data in any of those locations.  If we havn't then
        ! we store this level's data and mark that range as filled.
        do kk = lo(3), hi(3)

           if (idir == 3) then
              index_lo = kk*r1
              index_hi = (kk+1)*r1 -1
           endif

           do jj = lo(2), hi(2)

              if (idir == 2) then
                 index_lo = jj*r1
                 index_hi = (jj+1)*r1 -1
              endif

              do ii = lo(1), hi(1)

                 if (idir == 1) then
                    index_lo = ii*r1
                    index_hi = (ii+1)*r1 -1
                 endif

                 ! calculate the average
                 if (dim == 2) then
                    if (any(imask(ii*r1:(ii+1)*r1-1, &
                                  jj*r1:(jj+1)*r1-1, &
                                  1))) then
                       ! apply the proper weighting such that data from coarse
                       ! and fine cells are given the same weight.  We have to
                       ! count data from coarse cells r1**(dim-1) times.  In
                       ! doing this we are treating coarse cells as if they are
                       ! broken up into the cells at the finest resolution all
                       ! containing the same value as the original coarse cell.
                       weight = r1**(dim-1)
                       if (do_favre) weight = p(ii,jj,kk,irho_comp_pass)*weight

                       do index = index_lo, index_hi
                          avg(index,:) = avg(index,:) + p(ii,jj,kk,:) * weight

                          cnt(index)   = cnt(index) + weight
                       enddo

                       imask(ii*r1:(ii+1)*r1-1, &
                             jj*r1:(jj+1)*r1-1, &
                             1) = .false.
                    endif

                 else if (dim == 3) then
                    if (any(imask(ii*r1:(ii+1)*r1-1, &
                                  jj*r1:(jj+1)*r1-1, &
                                  kk*r1:(kk+1)*r1-1))) then
                    
                       ! apply the proper weighting such that data from coarse
                       ! and fine cells are given the same weight.  We have to
                       ! count data from coarse cells r1**(dim-1) times.  In
                       ! doing this we are treating coarse cells as if they are
                       ! broken up into the cells at the finest resolution all
                       ! containing the same value as the original coarse cell.
                       weight = r1**(dim-1)
                       if (do_favre) weight = p(ii,jj,kk,irho_comp_pass)*weight

                       do index = index_lo, index_hi
                          avg(index,:) = avg(index,:) + p(ii,jj,kk,:) * weight
                          cnt(index)   = cnt(index) + weight
                       enddo

                       imask(ii*r1:(ii+1)*r1-1, &
                             jj*r1:(jj+1)*r1-1, &
                             kk*r1:(kk+1)*r1-1) = .false.
                    endif

                 else 
                    call bl_error("faverage does not support 1d!")
                 endif
              enddo
           enddo
        enddo

        call fab_unbind(pf, i, j)

     enddo

     ! adjust r1 for the next level
     if (i /= 1) r1 = r1*pf%refrat(i-1,1)
     
  enddo

  ! error checking; this should never happen with non-corrupted data
  if (any(cnt == 0)) call bl_error("pltfile contains zones with empty data!")

  ! normalize
   do i = 0, max_points-1
         avg(i,:) = avg(i,:) / cnt(i)
   enddo

  !============================================================
  ! compute the RMS quantities                                |
  !============================================================

  imask = .true.
  r1 = 1
  cnt = 0

  do i = pf%flevel, 1, -1

     do j = 1, nboxes(pf, i)

        call fab_bind_comp_vec(pf, i, j, icomps_pass)

        p => dataptr(pf, i, j)
        lo(:) = 1
        hi(:) = 1
        lo(1:dim) = lwb(get_box(pf, i, j))
        hi(1:dim) = upb(get_box(pf, i, j))

        do kk = lo(3), hi(3)

           if (idir == 3) then
              index_lo = kk*r1
              index_hi = (kk+1)*r1 - 1
           endif

           do jj = lo(2), hi(2)

              if (idir == 2) then
                 index_lo = jj*r1
                 index_hi = (jj+1)*r1 - 1
              endif

              do ii = lo(1), hi(1)

                 if (idir == 1) then
                    index_lo = ii*r1
                    index_hi = (ii+1)*r1 - 1
                 endif

                 if (dim == 2) then

                    if(any(imask(ii*r1:(ii+1)*r1-1, &
                                 jj*r1:(jj+1)*r1-1, &
                                 1))) then

                       weight = r1**(dim-1)
                    
                       do index = index_lo, index_hi
                          rms(index,:) = rms(index,:) + &
                               (p(ii,jj,kk,:) - avg(index,:))**2 * weight

                          cnt(index) = cnt(index) + weight
                       enddo

                       imask(ii*r1:(ii+1)*r1-1, &
                             jj*r1:(jj+1)*r1-1, &
                             1) = .false.

                    endif

                 else if (dim == 3) then

                    if(any(imask(ii*r1:(ii+1)*r1-1, &
                                 jj*r1:(jj+1)*r1-1, &
                                 kk*r1:(kk+1)*r1-1))) then
                       weight = r1**(dim-1)

                       do index = index_lo, index_hi
                          rms(index,:) = rms(index,:) + &
                               (p(ii,jj,kk,:) - avg(index,:))**2 * weight
                    
                          cnt(index) = cnt(index) + weight
                       enddo

                       imask(ii*r1:(ii+1)*r1-1, &
                             jj*r1:(jj+1)*r1-1, &
                             kk*r1:(kk+1)*r1-1) = .false.

                    endif

                 else 
                    call bl_error("faverage does not support 1d!")
                 endif

              enddo
           enddo
        enddo

        call fab_unbind(pf, i,j)

     enddo

     if (i/=1) r1 = r1*pf%refrat(i-1,1)

  enddo

  ! error checking
  if (any(cnt == 0)) call bl_error("pltfile contains zones with empty data!")

  ! normalize
  do i = 0, max_points-1
     rms(i,:) = sqrt(rms(i,:)/cnt(i))
  enddo

  !============================================================
  ! output                                                    |
  !============================================================
  
100  format("# time:", 1x, g24.12)

  ! flexible formats based on the number of variables
  ! number of columns = 2*nvars_pass + 1
  write(columns,'(i10)') 2*nvars_pass+1
  column_format = '("#",5x,a,5x,' // trim(columns) &
       // '(a,"_avg",5x,a,"_rms",5x))'
  data_format = '(1x,' // trim(columns) // '(g24.12e3,1x))'


  if (trim(outputfile) == "stdout") then
     uout = 6
  else
     uout = unit_new()
     open(unit=uout, file=outputfile, status='replace')
  endif

  write(uout, 100) pf%tm

  write(uout,trim(column_format)) dir, &
       (trim(pf%names(icomps_pass(j))), &
       trim(pf%names(icomps_pass(j))), j = 1, nvars_pass)

  do i = flo(idir), fhi(idir)
     write(uout,trim(data_format)) pf%plo(idir) + (i+HALF)*dx(idir), &
          (avg(i,j), rms(i,j), j=1, nvars_pass)
  enddo

  if (uout .ne. 6) close(unit=uout)

  call destroy(pf)
        
  
contains
  
  subroutine print_usage()

    implicit none

    print *, ''
    print *, 'Description: '
    print *, '  This program takes a specified pltfile and calculates the '
    print *, '  average and rms quantities as a function of "height" along a '
    print *, '  specified direction for a specified set of variables. '
    print *, ''
    print *, 'Usage: '
    print *, '  faverage <args> '
    print *, ''
    print *, 'Arguments: '
    print *, '  [-p|--pltfile]    <filename> : '
    print *, '      specifies the input pltfile; (required) '
    print *, ''
    print *, '  [-o|--outputfile] <filename> : '
    print *, '      specifies the output file; default is "output.dat"'
    print *, ''
    print *, '  [-d|--dir]        <integer>  : '
    print *, '      specifies "height" direction; default is 2 (y) for 2-d'
    print *, '      and 3 (z) for 3-d'
    print *, ''
    print *, '  [-v|--vars] <integer> <list> : '
    print *, '      specifies a specific <list> of variables to analyze. '
    print *, '      <integer> is the number of variables in the whitespace '
    print *, '      separated <list> of variable names.  any variables that '
    print *, '      are specified but are not in the pltfile will be skipped; '
    print *, '      default behaviour is to analyze all variables. '
    print *, ''
    print *, '  [-f|--favre]                 : '
    print *, '      toggles favre (density-weighted) averaging; default is to '
    print *, '      do normal averaging '
    print *, ''
    print *, '  [-q|--quiet]                 : '
    print *, '      perform averaging operations quietly'
    print *, ''

  end subroutine print_usage


end program faverage
