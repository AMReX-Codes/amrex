
subroutine SDC_feval_F(lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi,f, flo,fhi, &
                         dx,a,d,r,n) bind(C, name="SDC_feval_F")

  !  Compute the rhs terms
  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), domlo(2), domhi(2)
  integer philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2),flo(2), fhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(inout) :: f   (flo(1):fhi(1),flo(2):fhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  real(amrex_real), intent(in)    :: a,d,r
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
              f(i,j) =  a*((fluxx(i+1,j  ) - fluxx(i,j))/dx(1) &
                   + (fluxy(i  ,j+1) - fluxy(i,j))/dx(2))
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
              f(i,j) =  d*((fluxx(i+1,j  ) - fluxx(i,j))/dx(1) &
                   + (fluxy(i  ,j+1) - fluxy(i,j))/dx(2))
           end do
        end do

     case (2)  !  Second implicit piece (here it is reaction)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
!              f(i,j) =  r*phi(i,j)*(1.0d0-phi(i,j))*(0.5d0-phi(i,j))
              f(i,j) =  -r*phi(i,j)
             end do
          end do
     case default
        print *, 'bad case in advance_2d'
     end select
     
   end subroutine SDC_feval_F
   

   subroutine SDC_fcomp_reaction_F (lo, hi, domlo, domhi, phi, philo, phihi, &
        rhs, rhslo, rhshi, &
        f, flo,fhi, dtq,n) bind(C, name="SDC_fcomp_reaction_F")

  !  Solve for the reaction term
  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), domlo(2), domhi(2)
  integer philo(2), phihi(2)
  integer rhslo(2), rhshi(2)
  integer flo(2), fhi(2)
  real(amrex_real), intent(inout) :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in)    :: rhs  (rhslo(1):rhshi(1),rhslo(2):rhshi(2))  
  real(amrex_real), intent(in)     :: f   (flo(1):fhi(1),flo(2):fhi(2))
  real(amrex_real), intent(in)    :: dtq
  integer, intent(in)             :: n
  
  ! local variables
  integer i,j
  real(amrex_real) c

  !  Function 
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        c = 1.0d0-dtq*(1.0d0-phi(i,j))*(0.5d0-phi(i,j))
        c = 1.0d0+dtq
        phi(i,j) =  rhs(i,j)/c
     end do
  end do

end subroutine SDC_fcomp_reaction_F



