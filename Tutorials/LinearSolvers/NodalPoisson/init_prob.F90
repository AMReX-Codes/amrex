
module init_prob_module

  use amrex_fort_module, only : amrex_real
  use amrex_constants_module, only : two, three, four, M_PI
  implicit none

contains

  subroutine init_prob (lo, hi, rhs, rlo, rhi, phi, hlo, hhi, dx) &
       bind(c,name='init_prob')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, hlo, hhi
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(inout) :: phi(hlo(1):hhi(1),hlo(2):hhi(2),hlo(3):hhi(3))
    real(amrex_real), intent(in) :: dx(3)

    integer :: i,j,k
    real(amrex_real) :: x, y, z
    real(amrex_real), parameter :: tpi = two*M_PI
    real(amrex_real), parameter :: fpi = four*M_PI
    real(amrex_real), parameter :: fac = tpi*tpi*three


    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)

             x = i*dx(1)
             y = j*dx(2)
             z = k*dx(3)

             phi(i,j,k) = 1.d0 * (cos(tpi*x) * cos(tpi*y) * cos(tpi*z))  &
                &      + .25d0 * (cos(fpi*x) * cos(fpi*y) * cos(fpi*z))
                
             rhs(i,j,k) = -fac * (cos(tpi*x) * cos(tpi*y) * cos(tpi*z))  &
                  &       -fac * (cos(fpi*x) * cos(fpi*y) * cos(fpi*z))
          end do
       end do
    end do
  end subroutine init_prob

end module init_prob_module
