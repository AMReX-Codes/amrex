subroutine PROBINIT (init,name,namlen,problo,probhi)
  implicit none
  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  double precision, intent(in) :: problo(3), probhi(3)
end subroutine PROBINIT

! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  
! :::
! ::: NOTE:  cell size and locations are in velocity space 
! :::
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: plo,phi   => physical locations of the entire domain (not just this grid)
! ::: -----------------------------------------------------------
subroutine smc_initdata(level,time,lo,hi, &
     state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
     delta,plo,phi)

  implicit none

  integer, intent(in) :: level, lo(3), hi(3)
  integer, intent(in) :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision, intent(inout) :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,1)
  double precision, intent(in) :: time, delta(3), plo(3), phi(3)

        
end subroutine smc_initdata
