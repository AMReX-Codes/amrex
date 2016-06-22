program fdata
  implicit none

  call process_rate

end program fdata

subroutine process_rate

  use bl_space
  use bl_IO_module
  use plotfile_module

  implicit none

  real(kind=dp_t), allocatable :: tslice(:)
  real(kind=dp_t), allocatable :: bslice(:)
  real(kind=dp_t), allocatable :: allda(:,:)
  integer :: i, j, n, ii, jj, f
  real(kind=dp_t), pointer :: fb(:,:,:,:)
  real(kind=dp_t) :: dt, dy
  integer :: lo(2), hi(2), sz(2), plo(2), phi(2)
  real(kind=dp_t) :: psz(2)
  type(plotfile) pf
  integer :: unit
  character(len=64) :: fname
  integer :: narg, farg
  real(kind=dp_t) :: t0, y0, v0
  real(kind=dp_t) :: t1, y1, v1
  real(kind=dp_t) :: v_inflo
  real(kind=dp_t) :: r_inflo
  real(kind=dp_t) :: xc_inflo
  integer, parameter :: yvel_ind = 2
  integer, parameter :: rxc_ind  = 4
  real(kind=dp_t) :: a1, b1, b2, b3
  
  narg = command_argument_count()
  farg = 1

  v_inflo  = 2.934E3_dp_t
  xc_inflo = 0.5_dp_t
  r_inflo   = 1.0E7_dp_t

  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)
     select case (fname)
     case ('-v','--v_inflo')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) v_inflo
     case ('-x','--xc_inflo')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) v_inflo
     case ('-r','--r_inflo')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) r_inflo
     case default
        exit
     end select
     farg = farg + 1
  end do

  if ( farg > narg ) return

  unit = unit_new()
  t0 = 0.0
  y0 = 0.0
  v0 = 0.0

  ! we only need to consider the coarsest grid since the scheme is 
  ! conservative
  i = 1

  call get_command_argument(farg, value = fname)

  call build(pf, fname, unit)
  plo = lwb(plotfile_get_pd_box(pf, i))
  phi = upb(plotfile_get_pd_box(pf, i))
  psz = pf%phi - pf%plo
  sz = phi - plo + 1
  
  allocate(tslice(plo(1):phi(1)))
  allocate(bslice(plo(1):phi(1)))
  allocate(allda(plo(1):phi(1), plo(2):phi(2)))

  do f = farg, narg
     call get_command_argument(f, value = fname)
     if ( f /= farg ) call build(pf, fname, unit)
     tslice = 0.0_dp_t
     bslice = 0.0_dp_t
     allda = 0.0_dp_t
     do j = 1, nboxes(pf, i)
        call fab_bind_comp_vec(pf, i, j, (/ yvel_ind, rxc_ind/) )
        fb => dataptr(pf, i, j)
        do n = 1, 2
           lo(n) = lbound(fb,dim = n)
           hi(n) = ubound(fb,dim = n)
        end do
        if ( lo(2) == plo(2) ) then
           do ii = lo(1), hi(1)
              bslice(ii) = fb(ii, lo(2), 1, 1)
           end do
        end if
        if ( hi(2) == phi(2) ) then
           do ii = lo(1), hi(1)
              tslice(ii) = fb(ii, hi(2), 1, 1)
           end do
        end if
        do jj = lo(2), hi(2)
           do ii = lo(1), hi(1)
              allda(ii,jj) = fb(ii,jj, 1, 2)
           end do
        end do
        call fab_unbind(pf, 1, j)
     end do
     y1 = sum(allda)*product(psz)/product(sz)
     t1 = pf%tm
     v1 = sum(tslice*allda(:,phi(2)))*psz(1)/sz(1)
     ! print *, 'avg(tslice) = ', sum(tslice)/sz(1)
     ! print *, 'avt(alldaT) = ', sum(allda(:,phi(2)))/sz(1)
     if ( f /= farg ) then
        a1 =  -(y1 - y0)
        b1 =  +(r_inflo*xc_inflo)*v_inflo*psz(1)*(t1-t0)
        b2 =  -(v1 + v0)*(t1-t0)/2
        b3 = (t1-t0)*psz(1)*(r_inflo*xc_inflo)
        print '(5(g15.8, 2x))', (t0+t1)/2, a1/b3, b1/b3, b2/b3, (a1 + b1 + b2)/b3
     end if
     y0 = y1
     t0 = t1
     v0 = v1
     call destroy(pf)
  end do
end subroutine process_rate
