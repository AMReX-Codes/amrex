program fdata
  implicit none

  call fturbkin

end program fdata


! Various statistical properties of the flow

subroutine fturbkin

  use bl_space
  use bl_IO_module
  use filler_module
  use plotfile_module

  implicit none

  integer :: f, level
  integer :: plo(2), phi(2), Nx(2)
  type(plotfile) pf
  integer :: unit

  character(len=64) :: fname
  character(len=128) :: times_out_filename
  character(len=128) :: avges_out_filename
  integer :: narg, farg
  integer, parameter :: rho_ind = 3
  integer, parameter :: u_ind = 1
  integer, parameter :: v_ind = 2
  integer :: max_level
  real(kind=dp_t) :: fmx

  integer :: nc, ntime, j, nxs, nn, ii, jj

  integer :: xstrt, xstop, xincr

  integer :: comps(3)
  
  real(kind=dp_t), allocatable :: f_fab(:,:,:)
  real(kind=dp_t), allocatable :: a_cor(:,:,:)
  real(kind=dp_t), allocatable :: v_cor(:,:)

  real(kind=dp_t), allocatable :: a_fab(:,:)
  real(kind=dp_t), allocatable :: f_times(:)

  narg = command_argument_count()
  farg = 1
  nc = 3

  xincr = 10
  xstrt = 0
  xstop = huge(xstop)

  max_level = huge(max_level)
  times_out_filename = ''
  avges_out_filename = ''
  comps = (/ rho_ind, u_ind, v_ind /)

  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)
     select case (fname)
     case ('--max_level')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) max_level
     case ('--times_out_filename')
        farg = farg + 1
        call get_command_argument(farg, times_out_filename)
     case ('--avges_out_filename')
        farg = farg + 1
        call get_command_argument(farg, avges_out_filename)
     case ('--xincr')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) xincr
     case ('--xstrt')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) xstrt
     case ('--xstop')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) xstop
     case default
        exit
     end select
     farg = farg + 1
  end do

  if ( farg > narg ) return

  unit = unit_new()

  call get_command_argument(farg, value = fname)
  call build(pf, fname, unit)

  level = min(pf%flevel, max_level)
  plo = lwb(plotfile_get_pd_box(pf, level))
  phi = upb(plotfile_get_pd_box(pf, level))
  Nx  = phi - plo + 1

  ntime = narg - farg + 1

  ! xstop = min(Nx(1)/2, xstop)
  xstop = min(Nx(1)-1, xstop)
  nxs   = max((xstop-xstrt+xincr)/xincr, 0)

  print *, '{NTIME} = ', ntime
  print *, '{XSTRT,XSTOP,XINCR,NXS} = ', xstrt, xstop, xincr, nxs

  allocate(f_fab  (plo(1):phi(1), plo(2):phi(2), nc))
  allocate(a_cor  (plo(1):phi(1), plo(2):phi(2), nc))

  allocate(a_fab  (plo(2):phi(2), nc))
  allocate(v_cor  (plo(2):phi(2), 2))
  allocate(f_times(1:ntime))

  if ( xstrt /= 0 ) then
     call bl_error('xstrt really ought to be zero')
  end if

  if ( times_out_filename /= '' ) then
     open(unit = 10, file = times_out_filename, status = 'replace')
  end if
  if ( avges_out_filename /= '' ) then
     open(unit = 11, file = avges_out_filename, status = 'replace')
  end if

  write(unit=11, fmt=*) ntime
  write(unit=11, fmt=*) 5
  write(unit=11, fmt=*) Nx(2)

  do f = 1, ntime
     call get_command_argument(f + farg - 1, value = fname)
     if ( f /= 1 ) call build(pf, fname, unit)
     print '(i5,a,a)', f, ': ', trim(fname)

     f_times(f) = pf%tm

     print '(7x,a)', 'filling'
     call blow_out_to_fab(f_fab, plo, pf, comps(1:nc), level)

     do nn = 2, 3
        f_fab(:,:,nn) = f_fab(:,:,nn)*f_fab(:,:,1)
     end do

     print '(7x,a)', 'avg'
     a_fab(:,:) = sum(f_fab, dim=1)/Nx(1)

     print '(7x,a)', 'ccor'
     do ii = plo(1), phi(1)
        a_cor(ii,:,:) = f_fab(ii,:,:) - a_fab(:,:)
     end do

     do ii = plo(1), phi(1)
        a_cor(ii,:,2) = a_cor(ii,:,2)/a_fab(:,1)
        a_cor(ii,:,3) = a_cor(ii,:,3)/a_fab(:,1)
     end do

     v_cor(:,1) = sum(a_cor(:,:,2)*a_cor(:,:,2), DIM=1)/Nx(1)
     v_cor(:,2) = sum(a_cor(:,:,3)*a_cor(:,:,3), DIM=1)/Nx(1)

     write(unit=11, fmt=*) a_fab(:,:)
     write(unit=11, fmt=*) v_cor

     call destroy(pf)
  end do

  write(unit=10, fmt=*) f_times
  close(unit=10)
  close(unit=11)

end subroutine fturbkin
