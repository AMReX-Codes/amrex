subroutine init_phi(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi) bind(C, name="init_phi")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3), philo(3), phihi(3)
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
  real(amrex_real), intent(in   ) :: dx(3) 
  real(amrex_real), intent(in   ) :: prob_lo(3) 
  real(amrex_real), intent(in   ) :: prob_hi(3) 

  integer          :: i,j,k
  double precision :: x,y,z,r2

  double precision :: L(3), shift(3)
  double precision :: x_p, y_p, z_p

  L = prob_hi - prob_lo
  shift = L/2.0

  do k = lo(3), hi(3)
     z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
     do j = lo(2), hi(2)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
        do i = lo(1), hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

           ! r2 = ((x-0.25d0)**2 + (y-0.25d0)**2 + (z-0.25d0)**2) / 0.01d0
           ! phi(i,j,k) = 1.d0 + exp(-r2)

           ! Periodic shift
           z_p = z + shift(3)
           y_p = y + shift(2)
           x_p = x + shift(1)
           if (x_p > prob_hi(1)) then
              x_p = x_p - L(1)
           endif
           if (y_p > prob_hi(2)) then
              y_p = y_p - L(2)
           endif
           if (z_p > prob_hi(3)) then
              z_p = z_p - L(3)
           endif
           r2 = ((x_p-0.0d0)**2 + (y_p-0.0d0)**2 + (z_p-0.0d0)**2) / 0.10d0
           phi(i,j,k) = 1.d0 + exp(-r2)
           
        end do
     end do
  end do

end subroutine init_phi
