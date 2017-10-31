
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_BC_TYPES.H"

#define SDIM 1

! ::: -----------------------------------------------------------
! ::: This routine is intended to be a generic fill function
! ::: for cell-centered data.  It knows how to extrapolate
! ::: and reflect data and is used to supplement the problem-specific
! ::: fill functions which call it.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: q           <=  array to fill
! ::: lo,hi        => index extent of q array
! ::: domlo,domhi  => index extent of problem domain
! ::: dx           => cell spacing
! ::: xlo          => physical location of lower left hand
! :::	              corner of q array
! ::: bc	   => array of boundary flags bc(SPACEDIM,lo:hi)
! ::: 
! ::: NOTE: all corner as well as edge data is filled if not EXT_DIR
! ::: -----------------------------------------------------------

AMREX_DEVICE subroutine filcc(q,q_l1,q_h1,domlo,domhi,dx,xlo,bc)

  use amrex_filcc_module, only: filccn

  implicit none

  integer    q_l1, q_h1
  integer    domlo(SDIM), domhi(SDIM)
  integer    bc(SDIM,2)
  REAL_T     xlo(SDIM), dx(SDIM)
  REAL_T     q(q_l1:q_h1)

  integer :: q_lo(3), q_hi(3)

  q_lo = [q_l1, 0, 0]
  q_hi = [q_h1, 0, 0]

  call filccn(q_lo, q_hi, q, q_lo, q_hi, 1, domlo, domhi, dx, xlo, bc)

end subroutine filcc
