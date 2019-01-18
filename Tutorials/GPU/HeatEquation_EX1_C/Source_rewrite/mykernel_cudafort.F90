
module my_kernel_cuda_fort_module
    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    implicit none
contains

    !==================================
    !  CUDA Fortran offloading subroutines
    !==================================
AMREX_CUDA_FORT_DEVICE subroutine init_phi(lo, hi, phi, philo, phihi, dx, prob_lo) bind(C, name="init_phi")
      integer, intent(in) :: lo(3), hi(3), philo(3), phihi(3)
      real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
      real(amrex_real), intent(in   ) :: dx(3) 
      real(amrex_real), intent(in   ) :: prob_lo(3) 

      integer          :: i,j,k
      double precision :: x,y,z,r2

      do k = lo(3), hi(3)
         z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
         do j = lo(2), hi(2)
            y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
            do i = lo(1), hi(1)
               x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

               r2 = ((x-0.25d0)**2 + (y-0.25d0)**2 + (z-0.25d0)**2) / 0.01d0
               phi(i,j,k) = 1.d0 + exp(-r2)

            end do
         end do
      end do

    end subroutine init_phi


AMREX_CUDA_FORT_DEVICE subroutine compute_flux (lo, hi, domlo, domhi, phi, philo, phihi, &
                             fluxx, fxlo, fxhi, fluxy, fylo, fyhi, fluxz, fzlo, fzhi, &
                             dx) bind(C, name="compute_flux")
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer philo(3), phihi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)
      real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
      real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2), fxlo(3): fxhi(3))
      real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2), fylo(3): fyhi(3))
      real(amrex_real), intent(inout) :: fluxz( fzlo(1): fzhi(1), fzlo(2): fzhi(2), fzlo(3): fzhi(3))
      real(amrex_real), intent(in)    :: dx(3)
      
      ! local variables
      integer i,j,k

      ! x-fluxes
      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)+1
         fluxx(i,j,k) = ( phi(i,j,k) - phi(i-1,j,k) ) / dx(1)
      end do
      end do
      end do

      ! y-fluxes
      do k = lo(3), hi(3)
      do j = lo(2), hi(2)+1
      do i = lo(1), hi(1)
         fluxy(i,j,k) = ( phi(i,j,k) - phi(i,j-1,k) ) / dx(2)
      end do
      end do
      end do

      ! z-fluxes
      do k = lo(3), hi(3)+1
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         fluxz(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / dx(3)
      end do
      end do
      end do

    end subroutine compute_flux


AMREX_CUDA_FORT_DEVICE subroutine update_phi (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
                           fluxx, fxlo, fxhi, fluxy, fylo, fyhi, fluxz, fzlo, fzhi, &
                           dx, dt) bind(C, name="update_phi")
      integer lo(3), hi(3), polo(3), pohi(3), pnlo(3), pnhi(3), &
           fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)
      real(amrex_real), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2),polo(3):pohi(3))
      real(amrex_real), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2),pnlo(3):pnhi(3))
      real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
      real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
      real(amrex_real), intent(in   ) :: fluxz (fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
      real(amrex_real), intent(in)    :: dx(3)
      real(amrex_real), intent(in), value :: dt

      ! local variables
      integer i,j,k
      real(amrex_real) :: dtdx(3)

      dtdx = dt/dx

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)

         phinew(i,j,k) = phiold(i,j,k) &
              + dtdx(1) * (fluxx(i+1,j  ,k  ) - fluxx(i,j,k)) &
              + dtdx(2) * (fluxy(i  ,j+1,k  ) - fluxy(i,j,k)) &
              + dtdx(3) * (fluxz(i  ,j  ,k+1) - fluxz(i,j,k))

      end do
      end do
      end do

    end subroutine update_phi

end module my_kernel_cuda_fort_module

