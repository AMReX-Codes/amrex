! Process a 2-d gaussian entropy pulse and compare it to the analytic solution for calculating the l1 norm
! this is based off of Castro_radiation/fgaussianpulse.f90
!
! this is used to analyze data from the test_diffusion test problem
!
! NOTE: THIS IS NOT MULTILEVEL AWARE
!

program fnorm

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module

  implicit none

  type(plotfile) pf
  integer :: unit
  integer :: i, j, ii, jj
  real(kind=dp_t) :: xx, yy

  real(kind=dp_t) :: dx(MAX_SPACEDIM)
  real(kind=dp_t) :: npts

  real(kind=dp_t) :: xctr, yctr, time, dist

  ! data for the analytic solution
  real(kind=dp_t) :: t0, h0, hp, D, analytic_solution

  real(kind=dp_t), pointer :: p(:,:,:,:)

  integer :: rhoh_comp, rho_comp

  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)

  character(len=256) :: pltfile

  integer :: narg, farg
  character(len=256) :: fname

  real(kind=dp_t) :: l1

  unit = unit_new()


  ! set the defaults
  pltfile  = ''
  t0 = 0.1d0
  D = 0.32009079924505196d0
  hp = TEN
  h0 = ONE

  l1 = ZERO
  npts = ZERO

  narg = command_argument_count()

  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-p', '--pltfile')
        farg = farg + 1
        call get_command_argument(farg, value = pltfile)

     case ('-h0')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) h0

     case ('-t0')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) t0

     case ('-hp')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) hp

     case ('-D')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) D

     case default
        exit

     end select
     farg = farg + 1
  end do

  if ( len_trim(pltfile) == 0 ) then
     print *, "usage: fgaussianpulse args"
     print *, "args [-p|--pltfile]   plotfile   : plot file directory (required)"
     print *, "     [-h0] background_enthalpy; default is ZERO"
     print *, "     [-t0] initial time; default is 0.1"
     print *, "     [-hp] peak_enthalpy; default is TEN"
     print *, "     [-D] diffusion_coefficient; default is 0.32009079924505196d0 "
     print *, '     based on the dump_gnuplot_analysis in test_diffusion/varden.f90'
     stop
  end if

  call build(pf, pltfile, unit)

  do i = 1, pf%flevel
     call fab_bind_level(pf, i)
  end do

  time = pf%tm

  dx = plotfile_get_dx(pf,pf%flevel)

  ! find the center coordinate
  xctr = HALF * (pf%phi(1) - pf%plo(1))
  yctr = HALF * (pf%phi(2) - pf%plo(2))

  ! figure out the variable indices
  rhoh_comp = -1
  do i = 1, pf%nvars
     if (pf%names(i) == "rhoh") then
        rhoh_comp = i
        exit
     endif
  enddo
  if (rhoh_comp < 0) call bl_error("rhoh not found")

  rho_comp = -1
  do i = 1, pf%nvars
     if (pf%names(i) == "density") then
        rho_comp = i
        exit
     endif
  enddo
  if (rho_comp < 0) call bl_error("density not found")

  do i = pf%flevel, 1, -1

     do j = 1, nboxes(pf, i)
        lo = lwb(get_box(pf, i, j))
        hi = upb(get_box(pf, i, j))

        ! get a pointer to the current patch
        p => dataptr(pf, i, j)
        
        ! loop over all of the zones in the patch.  
        do jj = lbound(p,dim=2), ubound(p,dim=2)
           yy = (jj + HALF)*dx(2)

           do ii = lbound(p,dim=1), ubound(p,dim=1)
              npts = npts + ONE

              xx = (ii + HALF)*dx(1)

              dist = sqrt((xx-xctr)**2 + (yy-yctr)**2)
              analytic_solution = f(time,dist)

                 l1 = l1 + abs(analytic_solution - &
                      p(ii,jj,1,rhoh_comp)/p(ii,jj,1,rho_comp))
                 
           end do
        enddo

     end do

  end do

  ! normalize
  l1 = l1 / npts

  print *, time, l1

  do i = 1, pf%flevel
     call fab_unbind_level(pf, i)
  end do

  call destroy(pf)

contains
  
  function f(t,x) result(r)
    real(kind=dp_t) :: t, x
    real(kind=dp_t) :: r

    r = (hp-h0)*(t0/(t+t0))*dexp(-x*x/(FOUR*D*(t+t0))) + h0
  end function f
  

end program fnorm
