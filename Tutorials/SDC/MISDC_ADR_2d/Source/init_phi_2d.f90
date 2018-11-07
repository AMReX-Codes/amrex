subroutine init_phi(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi) bind(C, name="init_phi")
  !  Initialize the scalar field phi
  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2)
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in   ) :: dx(2) 
  real(amrex_real), intent(in   ) :: prob_lo(2) 
  real(amrex_real), intent(in   ) :: prob_hi(2) 

  integer          :: i,j
  double precision :: x,y,tupi
  tupi=3.14159265358979323846d0*2d0

  do j = lo(2), hi(2)
     y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
     do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)) * dx(1)

        phi(i,j) =sin(x*tupi)*sin(y*tupi)
     end do
  end do

end subroutine init_phi

subroutine err_phi(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi,a,d,r,time) bind(C, name="err_phi")
  !  Subtract the exact solution from phi.  This will only work for special initial conditions
  !  We use the exact discretized form for diffusion and reaction, and exact translation for advection
  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2)
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in   ) :: dx(2) 
  real(amrex_real), intent(in   ) :: prob_lo(2) 
  real(amrex_real), intent(in   ) :: prob_hi(2) 
  real(amrex_real), intent(in   ) :: a,d,r
  real(amrex_real), intent(in   ) :: time

  integer          :: i,j
  double precision :: x,y,sym,tupi, maxphi

  tupi=3.14159265358979323846d0*2d0

  !  Form the diffusion coefficient for the 2nd order Laplacian produces
  sym=d*(-2.0d0+2.0d0*cos(tupi*dx(1)))/(dx(1)*dx(1))
  sym=sym+d*(-2.0d0+2.0d0*cos(tupi*dx(2)))/(dx(2)*dx(2))
  sym=sym-r    !  Add reaction

  maxphi=-1.0
  do j = lo(2), hi(2)
     y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) +a*time
     do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)) * dx(1) + a*time
        phi(i,j) = phi(i,j)-sin(x*tupi)*sin(y*tupi)*exp(time*sym)
        if (abs(phi(i,j)) .gt. maxphi) maxphi=abs(phi(i,j))
     end do
  end do

  print *,'max error in phi=',maxphi
end subroutine err_phi
