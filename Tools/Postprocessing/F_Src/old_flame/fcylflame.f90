program fdata
  implicit none

  call process

end program fdata

subroutine process

  use bl_space
  use bl_IO_module
  use plotfile_module
  use filler_module

  implicit none

  integer :: comps(2), f
  real(kind=dp_t), pointer :: fb(:,:,:,:)
  integer :: lo(2), hi(2), plo(2), phi(2), tlo(2), thi(2)
  type(plotfile) pf
  integer :: unit
  character(len=64) :: fname
  integer :: narg, farg
  real(kind=dp_t), allocatable :: f_fab(:,:,:)
  real(kind=dp_t) :: Dx(2), r0, r1
  integer :: nxc, nx, n
  integer :: ntime, max_level, level
  real(kind=dp_t) :: t0, t1, vol, vel, m0, m1, vel_jbb
  real(kind=dp_t) :: PI
  integer :: carbon
  real(kind=dp_t) :: rho

  integer :: cut

  PI = 4.0_dp_t*atan(1.0_dp_t)
  cut = 1

  narg = command_argument_count()
  farg = 1
  carbon = -1
  rho  = 5e7_dp_t

  max_level = huge(max_level)
  comps = (/20,5/)                     ! Mg for the C flame

  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)
     select case (fname)
     case ('--cut')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) cut
     case ('-M') 
        carbon = 0
     case ('-m','--max_level')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) max_level
     case ('-r','--rho')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) rho
     case default
        exit
     end select
     farg = farg + 1
  end do

  if ( farg > narg ) return

  unit = unit_new()

  comps = comps + carbon

  call get_command_argument(farg, value = fname)
  call build(pf, fname, unit)

  level = min(pf%flevel, max_level)
  plo = lwb(plotfile_get_pd_box(pf, level))
  phi = upb(plotfile_get_pd_box(pf, level))

  Dx = (pf%phi - pf%plo)/real(phi-plo+1)

  nx = phi(1) - plo(1) + 1
  nxc = nx/cut
  ntime = narg - farg + 1

  t0 = pf%tm
  m0 = 0.0_dp_t
  do f = 1, ntime
     call get_command_argument(f + farg - 1, value = fname)
     if ( f /= 1 ) call build(pf, fname, unit)
     t1 = pf%tm
     vol = 0.0_dp_t
     m1  = 0.0_dp_t
     do n = 1, cut
        tlo(1) = (n-1)*nxc
        thi(1) = n*nxc - 1
        if ( n == cut ) thi(1) = nx - 1
        tlo(2) = plo(2)
        thi(2) = phi(2)
        allocate(f_fab(tlo(1):thi(1), tlo(2):thi(2), 2))
        call blow_out_to_sub_fab(f_fab, tlo, thi, pf, comps, level)
        vol = vol + sum(f_fab(:,:,1))*product(Dx)
        m1 =  m1  + sum(f_fab(:,:,2))*product(Dx)
        deallocate(f_fab)
     end do
     r1 = sqrt(vol/PI)
     vel = 0.0_dp_t; vel_jbb = 0.0_dp_t
     if ( f > 1 ) then 
        vel = (r1-r0)/(t1-t0)
        vel_jbb = (m1-m0)/(PI*(t1-t0)*(r1+r0)*rho)
     end if
     print *, vol, r1, t1, m1, vel, vel_jbb
     t0 = t1
     r0 = r1
     m0 = m1
     call destroy(pf)
  end do

end subroutine process
