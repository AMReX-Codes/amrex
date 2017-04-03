
module prob_module

  implicit none

  private

  public :: init_prob_data

contains

  subroutine init_prob_data(level, time, lo, hi, phi, phi_lo, phi_hi, &
       dx, prob_lo)

    use amrex_fort_module, only : amrex_spacedim, amrex_real
    
    implicit none
    integer, intent(in) :: level, lo(3), hi(3), phi_lo(3), phi_hi(3)
    real(amrex_real), intent(in) :: time
    real(amrex_real), intent(inout) :: phi(phi_lo(1):phi_hi(1), &
         &                                 phi_lo(2):phi_hi(2), &
         &                                 phi_lo(3):phi_hi(3))
    real(amrex_real), intent(in) :: dx(3), prob_lo(3)
    
    integer          :: i,j,k
    real(amrex_real) :: x,y,z,r2
    
    !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
             
             if ( amrex_spacedim .eq. 2) then
                r2 = ((x-0.5d0)**2 + (y-0.75d0)**2) / 0.01d0
                phi(i,j,k) = 1.d0 + exp(-r2)
             else
                r2 = ((x-0.5d0)**2 + (y-0.75d0)**2 + (z-0.5d0)**2) / 0.01d0
                phi(i,j,k) = 1.d0 + exp(-r2)
             end if
          end do
       end do
    end do
    !$omp end parallel do
    
  end subroutine init_prob_data

end module prob_module
