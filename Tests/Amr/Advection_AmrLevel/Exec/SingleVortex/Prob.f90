
subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_fort_module

  implicit none

  integer(c_int), intent(in) :: init, namlen
  integer(c_int), intent(in) :: name(namlen)
  real(amrex_real), intent(in) :: problo(*), probhi(*)

  ! nothing needs to be done here,
  ! since there are no extra inputs to be read from probin file

end subroutine amrex_probinit


subroutine initdata(level, time, lo, hi, &
     phi, phi_lo, phi_hi, &
     dx, prob_lo) bind(C, name="initdata")

  use amrex_fort_module
  implicit none
  integer(c_int),   intent(in)    :: level, lo(3), hi(3), phi_lo(3), phi_hi(3)
  real(amrex_real), intent(in)    :: time
  real(amrex_real), intent(inout) :: phi(phi_lo(1):phi_hi(1), &
       &                                 phi_lo(2):phi_hi(2), &
       &                                 phi_lo(3):phi_hi(3))
  real(amrex_real), intent(in)    :: dx(3), prob_lo(3)

  integer          :: dm
  integer          :: i,j,k
  real(amrex_real) :: x,y,z,r2

  if ( phi_lo(3) == 0 .and. phi_hi(3) == 0 ) then
     dm = 2
  else
     dm = 3
  end if

  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)

        z = prob_lo(3) + (real(k,kind=amrex_real)+0.5d0) * dx(3)
        y = prob_lo(2) + (real(j,kind=amrex_real)+0.5d0) * dx(2)

        ! The compiler should automatically convert this innermost loop to SIMD
        do i = lo(1), hi(1)

           x = prob_lo(1) + (real(i,kind=amrex_real)+0.5d0) * dx(1)

           if ( dm == 2 ) then
              r2 = ((x-0.5d0)**2 + (y-0.75d0)**2) / 0.01d0
              phi(i,j,k) = 1.d0 + exp(-r2)
           else
              r2 = ((x-0.5d0)**2 + (y-0.75d0)**2 + (z-0.5d0)**2) / 0.01d0
              phi(i,j,k) = 1.d0 + exp(-r2)
           end if ! dm

        end do ! i
     end do ! j
  end do ! k
  !$omp end parallel do

end subroutine initdata
