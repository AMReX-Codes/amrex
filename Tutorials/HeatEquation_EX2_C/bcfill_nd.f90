subroutine phifill(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,time,bc)

  implicit none
  integer          :: phi_lo(3),phi_hi(3)
  integer          :: bc(*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
  ! no physical boundaries to fill because it is all periodic
  return

end subroutine phifill

