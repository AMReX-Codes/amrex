subroutine evaluate_f(lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi,f, flo,fhi, &
                         dx,v,nu,npiece) bind(C, name="evaluate_f")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), domlo(2), domhi(2),npiece
  integer philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2),flo(2), fhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(inout) :: f   (flo(1):fhi(1),flo(2):fhi(2))
  real(amrex_real), intent(in)    :: dx(2),v,nu


  !  Decide which piece it is we are doing
  !  npiece = 0 means explicit
  !  npiece = 1 means implicit
  !  npiece = 2 means 2nd implicit in misdc
  select case (npiece)
  case(0)
    call evaluate_fexp (lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi,f, flo,fhi, &
                         dx,v)
  case(1)
    call evaluate_fimp (lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi,f, flo,fhi, &
                         dx,nu)
  case(2)
    call evaluate_fimp (lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi,f, flo,fhi, &
                         dx,nu)
  case default
     print *, 'error statement here'
  end select   

end subroutine evaluate_f

subroutine evaluate_fexp (lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi,f, flo,fhi, &
                         dx,v) bind(C, name="evaluate_fexp")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), domlo(2), domhi(2)
  integer philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2),flo(2), fhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(inout) :: f   (flo(1):fhi(1),flo(2):fhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  real(amrex_real), intent(in)    :: v


  
  ! local variables
  integer i,j

  ! x-fluxes
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)+1
     fluxx(i,j) = ( phi(i,j) + phi(i-1,j) ) / 2.0d0
  end do
  end do

  ! y-fluxes
  do j = lo(2), hi(2)+1
  do i = lo(1), hi(1)
     fluxy(i,j) = ( phi(i,j) + phi(i,j-1) ) / 2.0d0
  end do
  end do

  !  Explicit function value
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        f(i,j) = v*((fluxx(i+1,j  ) - fluxx(i,j))/dx(1) &
          + (fluxy(i  ,j+1) - fluxy(i,j))/dx(2))
  end do
  end do

end subroutine evaluate_fexp

!  Subroutine to return the value of the implicit part of the function
subroutine evaluate_fimp (lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi,f, flo,fhi, &
                         dx,nu) bind(C, name="evaluate_fimp")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), domlo(2), domhi(2)
  integer philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2),flo(2), fhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(inout) :: f   (flo(1):fhi(1),flo(2):fhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  real(amrex_real), intent(in)    :: nu

  
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

  !  Function value
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        f(i,j) =  nu*((fluxx(i+1,j  ) - fluxx(i,j))/dx(1) &
          + (fluxy(i  ,j+1) - fluxy(i,j))/dx(2))

  end do
  end do

end subroutine evaluate_fimp



