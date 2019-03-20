
! -----------------------------------------------------------
!> This routine is intended to be a generic fill function
!! for cell-centered data.  It knows how to extrapolate
!! and reflect data and is used to supplement the problem-specific
!! fill functions which call it.
!!
!! \param q           <=  array to fill
!! \param lo,hi        => index extent of q array
!! \param domlo,domhi  => index extent of problem domain
!! \param dx           => cell spacing
!! \param xlo          => physical location of lower left hand
!!	              corner of q array
!! \param bc	   => array of boundary flags bc(SPACEDIM,lo:hi)
!!
!! NOTE: all corner as well as edge data is filled if not EXT_DIR
! -----------------------------------------------------------

#ifndef AMREX_XSDK

subroutine filcc(q,q_l1,q_h1,domlo,domhi,dx,xlo,bc)

  use amrex_fort_module
  use amrex_filcc_module, only: filccn

  implicit none

  integer    q_l1, q_h1
  integer    domlo(1), domhi(1)
  integer    bc(1,2)
  real(amrex_real)     xlo(1), dx(1)
  real(amrex_real)     q(q_l1:q_h1)

  integer :: q_lo(3), q_hi(3)

  q_lo = [q_l1, 0, 0]
  q_hi = [q_h1, 0, 0]

  call filccn(q_lo, q_hi, q, q_lo, q_hi, 1, domlo, domhi, dx, xlo, bc)

end subroutine filcc

#endif
