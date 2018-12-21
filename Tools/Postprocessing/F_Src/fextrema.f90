! print out the extrema of each variable in the plotfiles


subroutine get_extrema(pf, nvar_max, nvar, var_indices, var_min, var_max)

  use plotfile_module
  use filler_module
  use bl_IO_module

  implicit none

  type(plotfile), intent(inout) :: pf

  integer, intent(in) :: nvar, nvar_max
  integer, intent(in) :: var_indices(nvar_max)
  real(kind=dp_t), intent(out) :: var_min(nvar_max), var_max(nvar_max)

  logical, allocatable :: imask(:,:,:)

  integer, dimension(MAX_SPACEDIM) :: flo, fhi, lo, hi
  real(kind=dp_t), pointer :: p(:,:,:,:)
  type(box) :: pd

  integer :: i, j, n, ii, jj, kk
  integer :: r1
  integer :: dim
  integer :: ivar

  dim = pf%dim

  var_min(:) = 1.e30
  var_max(:) = -1.e30

  ! loop over levels (finest to coarest) and find extrema, marking the
  ! mask when we've visited a spatial location.

  ! get the index bounds for the finest level
  pd = plotfile_get_pd_box(pf, pf%flevel)
  flo(1:dim) = lwb(pd)
  fhi(1:dim) = upb(pd)

  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here
  if (dim == 1) then
     allocate(imask(flo(1):fhi(1), 1, 1))
  else if (dim == 2) then
     allocate(imask(flo(1):fhi(1), flo(2):fhi(2), 1))
  else if (dim == 3) then
     allocate(imask(flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3)))
  endif

  ! r1 is the refinement factor between the current level and the FINEST level
  r1 = 1

  imask(:,:,:) = .true.

  do i = pf%flevel, 1, -1

     do j = 1, nboxes(pf, i)

        ! bind to a comp vector
        call fab_bind_comp_vec(pf, i, j, var_indices(1:nvar))

        ! get the data and box info
        p => dataptr(pf, i, j)
        lo(:) = 1
        hi(:) = 1
        lo(1:dim) = lwb(get_box(pf, i, j))
        hi(1:dim) = upb(get_box(pf, i, j))

        ! loop over all the zones in the current box

        ! Here, we convert the cell-centered indices at the current
        ! level into the corresponding RANGE on the finest level, and
        ! test if we've stored data in any of those locations.  If we
        ! havn't then we store this level's data and mark that range
        ! as filled.
        do kk = lo(3), hi(3)
           do jj = lo(2), hi(2)
              do ii = lo(1), hi(1)

                 ! calculate the extrema
                 if (dim == 1) then
                    if (any(imask(ii*r1:(ii+1)*r1-1, 1, 1))) then
                       do n = 1, nvar
                          ivar = var_indices(n)
                          var_min(ivar) = min(var_min(ivar), p(ii, 1, 1, n))
                          var_max(ivar) = max(var_max(ivar), p(ii, 1, 1, n))
                       end do

                       imask(ii*r1:(ii+1)*r1-1, 1, 1) = .false.
                    endif

                 else if (dim == 2) then
                    if (any(imask(ii*r1:(ii+1)*r1-1, jj*r1:(jj+1)*r1-1, 1))) then
                       do n = 1, nvar
                          ivar = var_indices(n)
                          var_min(ivar) = min(var_min(ivar), p(ii, jj, 1, n))
                          var_max(ivar) = max(var_max(ivar), p(ii, jj, 1, n))
                       end do

                       imask(ii*r1:(ii+1)*r1-1, jj*r1:(jj+1)*r1-1, 1) = .false.
                    endif

                 else if (dim == 3) then
                    if (any(imask(ii*r1:(ii+1)*r1-1, &
                                  jj*r1:(jj+1)*r1-1, &
                                  kk*r1:(kk+1)*r1-1))) then
                       do n = 1, nvar
                          ivar = var_indices(n)
                          var_min(n) = min(var_min(ivar), p(ii, jj, kk, n))
                          var_max(n) = max(var_max(ivar), p(ii, jj, kk, n))
                       end do

                       imask(ii*r1:(ii+1)*r1-1, &
                             jj*r1:(jj+1)*r1-1, &
                             kk*r1:(kk+1)*r1-1) = .false.
                    endif

                 endif
              end do
           end do
        end do

        call fab_unbind(pf, i, j)

     end do

     ! adjust r1 for the next level
     if (i /= 1) r1 = r1*pf%refrat(i-1,1)

  end do

end subroutine get_extrema


program fextrema

  use plotfile_module
  use filler_module
  use bl_IO_module

  implicit none

  integer f, farg
  integer :: n

  integer :: unit

  type(plotfile) :: pf
  integer narg
  character(len=256) ::  fname, varnames, temp

  integer :: ntime, numvars
  integer :: max_level

  character (len = 3) :: minname(50), maxname(50)

  logical :: single

  integer :: i, idx
  integer, allocatable :: var_indices(:)

  real(kind=dp_t), allocatable :: vvmin(:), vvmax(:)

  single = .FALSE.

  data minname / 50*"min"/
  data maxname / 50*"max"/

  narg = command_argument_count()


  farg = 1

  varnames = ''

  do while (farg <= narg)

     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-s','--single')
        single = .TRUE.

     case ('-v','--variable')
        farg = farg + 1
        call get_command_argument(farg, value = varnames)

     case default
        exit

     end select

     farg = farg + 1

  enddo

  if ( farg > narg ) then
     print *, " "
     print *, "Report the extrema (min/max) for each variable in a plotfile"
     print *, "usage: "
     print *, "   fextrema [-s|--single] {[-v|--variable] name} plotfiles"
     print *, " "
     print *, "By default, the variable information is specified in columns, one line"
     print *, "per file."
     print *, " "
     print *, "  -s or --single : for each plotfile, each variable's information is"
     print *, "                   printed on a separate line.  This is the behavior"
     print *, "                   when only 1 plotfile is specified"
     print *, " "
     print *, "  -v names : output information only for specified variables, given"
     print *, "             as a space-spearated string"
     print *
     stop
  endif

  unit = unit_new()

  ! ntime is the number of files to loop over
  ntime = narg - farg  + 1

  if (ntime == 1) then
     single = .true.
  endif

  do f = 1, ntime

     call get_command_argument(f + farg - 1, value = fname)
     call build(pf, fname, unit)

     ! if we are outputting specific variables, find their indices
     if (f == 1) then
        allocate(var_indices(pf%nvars))
        var_indices(:) = 0

        if (varnames /= '') then
           ! we are doing only some of the variables

           numvars = 0
           do while (.not. trim(varnames) == "")
              numvars = numvars + 1

              idx = index(varnames, " ")
              temp = varnames(:idx)
              varnames = trim(adjustl(varnames(idx+1:)))

              var_indices(numvars) = plotfile_var_index(pf, trim(temp))
           end do
        else

           ! we are doing all variables
           numvars = pf%nvars
           do n = 1, numvars
              var_indices(n) = n
           end do
        end if

        ! We allocate this to be the full size for all variables in the plotfile
        allocate(vvmin(pf%nvars))
        allocate(vvmax(pf%nvars))

     end if

     ! get the extrema.
     call get_extrema(pf, pf%nvars, numvars, var_indices, vvmin, vvmax)

     if (single) then

        write (*,*) 'plotfile = ', trim(fname)
        write (*,*) 'time = ', pf%tm
        write (*,200) "variable", "minimum value", "maximum value"


        do n = 1, numvars
           write (*,201) pf%names(var_indices(n)), vvmin(var_indices(n)), vvmax(var_indices(n))
        enddo
        write (*,*) " "

     else

        ! we are doing a single plotfile per line with extrema in columns

        if (f == 1) then
           ! print the header
           write (*,100) "time", (pf%names(var_indices(n)), n=1,numvars)
           write (*,101) (minname(var_indices(n)), maxname(var_indices(n)), n=1,numvars)
        endif

        write (*,102) pf%tm, (vvmin(var_indices(n)), vvmax(var_indices(n)), n=1,numvars)

     end if

     call destroy(pf)

  enddo

100 format ("#", 1x, a22, 1x, 50("|",a44))
101 format ("#", 1x, 22x, 1x, 50("|",8x,a3,9x,1x,8x,a3,9x,3x))
102 format (1x, g22.11, 1x, 50(g20.10, 1x, g20.10, 3x))

200 format (1x, a22, 1x, a22, 1x, a22)
201 format (1x, a22, 1x, g22.11, 1x, g22.11)

end program fextrema
