
subroutine PROBINIT (init,name,namlen,problo,probhi)

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  double precision, intent(in) :: problo(*), probhi(*)

  ! nothing to be done for this problem
end subroutine PROBINIT


subroutine initdata(level, time, lo, hi, &
     phi, phi_lo, phi_hi, &
     dx, problo, probhi)

end subroutine initdata
