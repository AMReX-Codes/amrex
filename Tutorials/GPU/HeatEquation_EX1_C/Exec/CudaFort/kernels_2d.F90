
module my_kernel_cuda_fort_module
    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    implicit none
contains

    !==================================
    !  CUDA Fortran offloading subroutines
    !==================================
AMREX_CUDA_FORT_DEVICE subroutine init_phi(lo, hi, phi, philo, phihi, dx, prob_lo) bind(C, name="init_phi")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2)
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in   ) :: dx(2) 
  real(amrex_real), intent(in   ) :: prob_lo(2) 

  integer          :: i,j
  double precision :: x,y,r2

  do j = lo(2), hi(2)
     y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
     do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

        r2 = ((x-0.25d0)**2 + (y-0.25d0)**2) / 0.01d0

        phi(i,j) = 1.d0 + exp(-r2)

     end do
  end do

end subroutine init_phi


AMREX_CUDA_FORT_DEVICE subroutine compute_flux (lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi, &
                         dx) bind(C, name="compute_flux")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), domlo(2), domhi(2)
  integer philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)

  ! local variables
  integer i,j

  ! x-fluxes
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)+1
     fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx(1)
  end do
  end do

  ! y-fluxes
  do j = lo(2), hi(2)+1
  do i = lo(1), hi(1)
     fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx(2)
  end do
  end do

end subroutine compute_flux


AMREX_CUDA_FORT_DEVICE subroutine update_phi (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
                       fluxx, fxlo, fxhi, fluxy, fylo, fyhi, &
                       dx, dt) bind(C, name="update_phi")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2))
  real(amrex_real), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2))
  real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
  real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  real(amrex_real), intent(in), value :: dt

  ! local variables
  integer i,j
  real(amrex_real) :: dtdx(2)

  dtdx = dt/dx

  do j = lo(2), hi(2)
  do i = lo(1), hi(1)

     phinew(i,j) = phiold(i,j) &
          + dtdx(1) * (fluxx(i+1,j  ) - fluxx(i,j)) &
          + dtdx(2) * (fluxy(i  ,j+1) - fluxy(i,j))

  end do
  end do

end subroutine update_phi

end module my_kernel_cuda_fort_module
