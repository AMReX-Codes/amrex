subroutine nullfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) bind(C, name="nullfill")
  implicit none
  integer          :: adv_lo(3),adv_hi(3)
  integer          :: bc(*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
  ! no physical boundaries to fill because it is all periodic
  return
end subroutine nullfill
