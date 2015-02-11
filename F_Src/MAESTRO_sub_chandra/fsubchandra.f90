! Main driver for analysis routines for the 3-d sub_chandra MAESTRO problem 
! This code borrows from fsedov3d_sph.f90 and fwdconvect.f90
! See fsubchandra_mod.f90 for more.
!
! usage: fsubchandra [args] plotfile"
!          args: [-s|--slicefile] <slice file> : specify output file"
!                --globals-only                : only output global quantities "
!
! --globals_only will just print T_peak, it's location, and enuc_peak
! and its location to stdout.
!
! For now you have the choice between globals only or all analysis 
! (for now this is just averages and hotspot calculations).  If this proves
! inconvenient then see TODO list.
!
!

!TODO:
! 1) XXAdd ability to handle full star simulations
! 2) Add ability to select specific analysis to carry out

program fsubchandra
  !Modules
  use plotfile_module
  use subchandra
  use eos_module, only: eos_init, eos_finalize
  use network, only: network_init, network_finalize 
  implicit none

  !Data, variables
  type(plotfile) :: pf
  type(state_comps) :: sc
  type(geometry) :: geo
  type(radial_averages) :: radav
  type(globals) :: glb
  type(hheap) :: hspots
  type(temp_hist) :: thist

  real(kind=dp_t) :: hpcnt
  integer :: i
  character(len=256) :: slicefile, pltfile
  logical :: globals_only, fullstar

  !Read command line args, build pf
  call parse_args(slicefile, pltfile, globals_only, fullstar, hpcnt, pf)

  !Initialize the component indices
  call init_comps(pf, sc)

  !Initialize geometry attributes 
  call init_geometry(pf, globals_only, fullstar, geo)

  !Initialize microphysics
  call eos_init()
  call network_init()

  !This is the main analysis subroutine. It loops over the entire computational
  !domain applying calculations (averages, globals, etc...) based on the arguments passed.
  if(globals_only) then
    print *, 'analyze globals only...'
    call analyze(pf, geo, sc, thist, glb=glb)
  else if(hpcnt > 0.0) then
    !Only calculate hotspot statistics if hpcnt > 0.0
    print *, 'analyze with hotspots, hpcnt = ', hpcnt
    call analyze(pf, geo, sc, thist, glb=glb, radav=radav, hh=hspots, hheap_frac_input=hpcnt)
  else
    print *, 'analyze without hotspots...'
    call analyze(pf, geo, sc, thist, glb=glb, radav=radav)
  end if

  !Write output
  if(globals_only) then
    call writeout(pf, slicefile, geo, sc, thist, glb=glb)
  else
    call writeout(pf, slicefile, geo, sc, thist, glb=glb, radav=radav, hh=hspots)
  end if

  !Clean up
  do i = 1, pf%flevel
     call fab_unbind_level(pf, i)
  end do

  call network_finalize()
  call eos_finalize()
  call destroy(pf)
end program fsubchandra
