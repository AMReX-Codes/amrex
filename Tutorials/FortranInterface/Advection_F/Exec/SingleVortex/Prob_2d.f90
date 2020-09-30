
module prob_module

  implicit none

  private

  public :: init_prob_data, init_part_data

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

  subroutine init_part_data(pc, lev, mfi, lo, hi, dx, prob_lo)

    use amrex_fort_module, only : amrex_spacedim, amrex_real
    use amrex_particlecontainer_module, only: amrex_particlecontainer, amrex_particle, &
         amrex_get_next_particle_id, amrex_get_cpu, amrex_set_particle_id, amrex_set_particle_cpu
    use amrex_multifab_module, only : amrex_mfiter
    
    implicit none
    type(amrex_particlecontainer), intent(inout) :: pc
    integer, intent(in) :: lev
    type(amrex_mfiter), intent(in) :: mfi
    integer, intent(in) :: lo(2), hi(2)
    real(amrex_real), intent(in) :: dx(2), prob_lo(2)

    integer          :: i,j
    real(amrex_real) :: x,y
    type(amrex_particle) :: p
    
    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

          p%pos(1) = x
          p%pos(2) = y

          p%vel = 0.d0

          call amrex_set_particle_id(amrex_get_next_particle_id(), p)
          call amrex_set_particle_cpu(amrex_get_cpu(), p)

          call pc%add_particle(lev, mfi, p)

       end do
    end do
    
  end subroutine init_part_data

  
end module prob_module
