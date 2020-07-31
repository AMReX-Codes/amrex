
subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  double precision, intent(in) :: problo(*), probhi(*)

  ! nothing needs to be done here

end subroutine amrex_probinit


subroutine initdata(level, time, lo, hi, &
     phi, phi_lo, phi_hi, &
     dx, prob_lo) bind(C, name="initdata")

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), phi_lo(3), phi_hi(3)
  double precision, intent(in) :: time
  double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1), &
       &                                 phi_lo(2):phi_hi(2), &
       &                                 phi_lo(3):phi_hi(3))
  double precision, intent(in) :: dx(3), prob_lo(3)

  integer          :: dm
  integer          :: i,j,k
  double precision :: x,y,z,r2
  
  if (phi_lo(3) .eq. 0 .and. phi_hi(3) .eq. 0) then
     dm = 2
  else
     dm = 3
  end if

  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
        do i=lo(1),hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
           r2 = sqrt((x-0.5d0)**2 + (y-0.75d0)**2)
           if(r2 > 0.15) then
              phi(i,j,k) = 0.d0
           elseif((abs(2*(x - 0.5d0)) < 0.05).and.((y - 0.75d0) + r2 < 0.25)) then
              phi(i,j,k) = 0.d0
           else
              phi(i,j,k) = 1.d0
           endif
        end do
     end do
  end do
  !$omp end parallel do

end subroutine initdata
