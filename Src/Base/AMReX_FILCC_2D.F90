
#ifndef AMREX_XSDK

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

subroutine filcc(q,q_l1,q_l2,q_h1,q_h2,domlo,domhi,dx,xlo,bc)

  use amrex_fort_module
  use amrex_filcc_module, only: filccn

  implicit none

  integer    q_l1, q_l2, q_h1, q_h2
  integer    domlo(2), domhi(2)
  integer    bc(2,2)
  real(amrex_real)     xlo(2), dx(2)
  real(amrex_real)     q(q_l1:q_h1,q_l2:q_h2)

  integer :: q_lo(3), q_hi(3)

  q_lo = [q_l1, q_l2, 0]
  q_hi = [q_h1, q_h2, 0]

  call filccn(q_lo, q_hi, q, q_lo, q_hi, 1, domlo, domhi, dx, xlo, bc)

end subroutine filcc

subroutine hoextraptocc(q,q_l1,q_l2,q_h1,q_h2,domlo,domhi,dx,xlo)

  use amrex_fort_module
  use amrex_filcc_module, only : amrex_hoextraptocc_2d

  implicit none

  integer    q_l1, q_l2, q_h1, q_h2
  integer    domlo(2), domhi(2)
  real(amrex_real)     xlo(2), dx(2)
  real(amrex_real)     q(q_l1:q_h1,q_l2:q_h2)

  call amrex_hoextraptocc_2d(q,q_l1,q_l2,q_h1,q_h2,domlo,domhi,dx,xlo)

end subroutine hoextraptocc

#endif
