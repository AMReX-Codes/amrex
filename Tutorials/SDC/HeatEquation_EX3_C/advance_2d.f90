
subroutine compute_flux (lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi,f, flo,fhi, &
                         dx,n) bind(C, name="compute_flux")

  !  Compute the flux for the nth rhs term

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), domlo(2), domhi(2)
  integer philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2),flo(2), fhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(inout) :: f   (flo(1):fhi(1),flo(2):fhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  integer, intent(in) :: n

  
  ! local variables
  integer i,j

  select case(n)
     case (0)  !  Explicit term (here it is advection)
        ! x-fluxes
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1
              fluxx(i,j) = ( phi(i,j) + phi(i-1,j) ) / 2.0d0 
              fluxx(i,j) = (-phi(i+1,j)+ 7.0d0*( phi(i,j) + phi(i-1,j) )-phi(i-2,j)) / 12.0d0
           end do
        end do
        
        ! y-fluxes
        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)
              fluxy(i,j) =(-phi(i,j+1)+7.0d0*( phi(i,j) + phi(i,j-1) )-phi(i,j-2)) / 12.0d0
           end do
        end do

        !  Function value is divergence of flux
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              f(i,j) =  (fluxx(i+1,j  ) - fluxx(i,j))/dx(1) &
                   + (fluxy(i  ,j+1) - fluxy(i,j))/dx(2)
              f(i,j)=0.0d0
           end do
           
        end do
     case (1)  !  First implicit piece (here it is diffusion)
        ! x-fluxes
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1
              fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx(1)
!              fluxx(i,j) = ( -phi(i+1,j)  +15.0d0*(phi(i,j) - phi(i-1,j)) + phi(i-2,j) ) /(12.0d0*dx(1))
           end do
        end do
        
        ! y-fluxes
        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)
              fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx(2)
!              fluxy(i,j) = ( -phi(i,j+1)  +15.0d0*(phi(i,j) - phi(i,j-1)) + phi(i,j-2) ) /(12.0d0*dx(2))
           end do
        end do

        !  Function value is divergence of flux
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              f(i,j) =  (fluxx(i+1,j  ) - fluxx(i,j))/dx(1) &
                   + (fluxy(i  ,j+1) - fluxy(i,j))/dx(2)
              f(i,j)=f(i,j)*0.1d0
           end do
        end do
     case default
        print *, 'bad case in advance_2d'
     end select
     
end subroutine compute_flux



