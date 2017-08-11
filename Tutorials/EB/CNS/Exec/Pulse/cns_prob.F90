
subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(*), probhi(*)
  ! nothing needs to be done here
end subroutine amrex_probinit


subroutine cns_initdata(level, time, lo, hi, u, ulo, uhi, dx, prob_lo) bind(C, name="cns_initdata")
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer, intent(in) :: level, lo(3), hi(3), ulo(3), uhi(3)
  real(rt), intent(in) :: time
  real(rt), intent(inout) :: u(ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3))
  real(rt), intent(in) :: dx(3), prob_lo(3)
  
end subroutine cns_initdata
