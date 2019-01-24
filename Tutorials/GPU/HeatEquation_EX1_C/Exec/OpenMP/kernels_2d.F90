
module my_kernel_omp_module
    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    implicit none
contains

    !==================================
    !  OpenMP offloading subroutines
    !==================================

    subroutine init_phi_omp(lo, hi, phi, philo, phihi, dx, prob_lo)&
        bind(C, name="init_phi_omp")
      integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2)
      real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
      real(amrex_real), intent(in   ) :: dx(2), prob_lo(2) 

      integer          :: i,j
      double precision :: x,y,r2

      !$omp target teams distribute parallel do collapse(2) schedule(static,1) is_device_ptr(phi)
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
            r2 = ((x-0.25d0)**2 + (y-0.25d0)**2) / 0.01d0
            phi(i,j) = 1.d0 + exp(-r2)
         end do
      end do

    end subroutine init_phi_omp

    subroutine compute_flux_x_omp(lo, hi, fluxx, f_lo, f_hi, &
                                  phi, p_lo, p_hi, dxinv)    &
        bind(c,name="compute_flux_x_omp")
    integer(c_int), intent(in)   :: lo(2), hi(2), f_lo(2), f_hi(2), p_lo(2), p_hi(2)
    real(amrex_real), intent(in) :: phi(p_lo(1):p_hi(1), p_lo(2):p_hi(2))
    real(amrex_real), intent(in), value :: dxinv
    real(amrex_real), intent(inout) :: fluxx(f_lo(1):f_hi(1), f_lo(2):f_hi(2))
    
    integer(c_int) :: i,j
    !$omp target teams distribute parallel do collapse(2) schedule(static,1) is_device_ptr(phi, fluxx)
    ! x-fluxes
    do      j=lo(2), hi(2)
        do  i=lo(1), hi(1)
            fluxx(i,j) = dxinv* ( phi(i,j)-phi(i-1,j))
        end do
    end do

    end subroutine compute_flux_x_omp


    subroutine compute_flux_y_omp(lo, hi, fluxy, f_lo, f_hi, &
                                  phi, p_lo, p_hi, dyinv)    &
        bind(c,name="compute_flux_y_omp")
    integer(c_int), intent(in)   :: lo(2), hi(2), f_lo(2), f_hi(2), p_lo(2), p_hi(2)
    real(amrex_real), intent(in) :: phi(p_lo(1):p_hi(1), p_lo(2):p_hi(2))
    real(amrex_real), intent(in), value :: dyinv
    real(amrex_real), intent(inout) :: fluxy(f_lo(1):f_hi(1), f_lo(2):f_hi(2))
    
    integer(c_int) :: i,j
    !$omp target teams distribute parallel do collapse(2) schedule(static,1) is_device_ptr(phi, fluxy)
    do      j=lo(2), hi(2)
        do  i=lo(1), hi(1)
            fluxy(i,j) = dyinv* ( phi(i,j)-phi(i,j-1))
        end do
    end do

    end subroutine compute_flux_y_omp

    subroutine update_phi_omp(lo,hi,&
                              fluxx,fxlo,fxhi, &
                              fluxy,fylo,fyhi, &
                              phi_old,polo,pohi, &
                              phi_new,pnlo,pnhi, &
                              dt,dxinv,dyinv) &
          bind(C, name="update_phi_omp")

      integer(c_int), intent(in)      :: lo(2), hi(2), polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
      real(amrex_real), intent(in)    :: phi_old(polo(1):pohi(1),polo(2):pohi(2))
      real(amrex_real), intent(inout) :: phi_new(pnlo(1):pnhi(1),pnlo(2):pnhi(2))
      real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
      real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2))
      real(amrex_real), intent(in), value :: dt, dxinv, dyinv

      ! local variables
      integer(c_int) :: i,j

      !$omp target teams distribute parallel do collapse(2) schedule(static,1) is_device_ptr(phi_old, phi_new, fluxx, fluxy)
      do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             phi_new(i,j) = phi_old(i,j) &
                  + dxinv * dt * (fluxx(i+1,j  ) - fluxx(i,j)) &
                  + dyinv * dt * (fluxy(i  ,j+1) - fluxy(i,j))
          end do
      end do

    end subroutine update_phi_omp
end module my_kernel_omp_module
