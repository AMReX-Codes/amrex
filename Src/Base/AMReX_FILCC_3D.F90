
#ifndef AMREX_XSDK

! -----------------------------------------------------------
!> This routine is intended to be a generic fill function
!! for cell centered data.  It knows how to exrapolate,
!! and reflect data and can be used to suppliment problem
!! specific fill functions (ie. EXT_DIR).
!!
!! \param q        <=  array to fill
!! \param q_l1,q_l2,q_l3,q_h1,q_h2,q_h3   => index extent of q array
!! \param domlo,hi  => index extent of problem domain
!! \param dx        => cell spacing
!! \param xlo       => physical location of lower left hand
!!	           corner of q array
!! \param bc	=> array of boundary flags bc(SPACEDIM,lo:hi)
!!
!! NOTE: corner data not used in computing soln but must have
!!       reasonable values for arithmetic to live
! -----------------------------------------------------------

subroutine filcc(q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,domlo,domhi,dx,xlo,bc)

  use amrex_fort_module, only: rt => amrex_real
  use amrex_filcc_module, only: filccn

  implicit none

  integer,  intent(in   ) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  integer,  intent(in   ) :: domlo(3), domhi(3)
  real(rt), intent(in   ) :: xlo(3), dx(3)
  real(rt), intent(inout) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
  integer,  intent(in   ) :: bc(3,2)

  integer :: q_lo(3), q_hi(3)

  q_lo = [q_l1, q_l2, q_l3]
  q_hi = [q_h1, q_h2, q_h3]

  call filccn(q_lo, q_hi, q, q_lo, q_hi, 1, domlo, domhi, dx, xlo, bc)

end subroutine filcc



subroutine hoextraptocc(q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,domlo,domhi,dx,xlo)

  use amrex_fort_module, only: rt => amrex_real
  use amrex_filcc_module, only : amrex_hoextraptocc_3d

  implicit none

  integer,  intent(in   ) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  integer,  intent(in   ) :: domlo(3), domhi(3)
  real(rt), intent(in   ) :: xlo(3), dx(3)
  real(rt), intent(inout) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)

  call amrex_hoextraptocc_3d(q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,domlo,domhi,dx,xlo)

end subroutine hoextraptocc

#endif
