
subroutine PROBINIT (init,name,namlen,problo,probhi)

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  double precision, intent(in) :: problo(*), probhi(*)

  ! nothing to be done for this problem
end subroutine PROBINIT


subroutine initdata(level, time, lo, hi, &
     phi, phi_lo, phi_hi, &
     dx, prob_lo, prob_hi)
  integer, intent(in) :: level, lo(3), hi(3), phi_lo(3), phi_hi(3)
  double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1), &
       &                                 phi_lo(2):phi_hi(2), &
       &                                 phi_lo(3):phi_hi(3))
  double precision, intent(in) :: dx(3), prob_lo(3), prob_hi(3)

  integer          :: i,j,k
  double precision :: x,y,z,r1,r2,r3
  
  print *, lo, hi, phi_lo, phi_hi
  print *, dx
  print *, prob_lo, prob_hi

  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
        do i=lo(1),hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
           
           r1 = ((x-0.5d0)**2 + (y-0.5d0)**2 + (z-0.0d0)**2) / 0.01d0
           r2 = ((x-0.0d0)**2 + (y-0.0d0)**2 + (z-0.0d0)**2) / 0.01d0
           r3 = ((x+0.5d0)**2 + (y+0.5d0)**2 + (z-0.0d0)**2) / 0.01d0
           
           phi(i,j,k) = 1.d0 + exp(-r1) + exp(-r2) + exp(-r3)
           
        end do
     end do
  end do
  !$omp end parallel do

end subroutine initdata
