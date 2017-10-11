
module amrex_mlabeclap_2d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlabeclap_adotx

contains

  subroutine amrex_mlabeclap_adotx (lo, hi, y, ylo, yhi, x, xlo, xhi, a, alo, ahi, &
       bx, bxlo, bxhi, dxinv, alpha, beta) &
       bind(c,name='amrex_mlabeclap_adotx')
    integer, dimension(1), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, alo, ahi, bxlo, bxhi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), value, intent(in) :: alpha, beta
    real(amrex_real), intent(inout) ::  y( ylo(1): yhi(1))
    real(amrex_real), intent(in   ) ::  x( xlo(1): xhi(1))
    real(amrex_real), intent(in   ) ::  a( alo(1): ahi(1))
    real(amrex_real), intent(in   ) :: bx(bxlo(1):bxhi(1))
    
    integer :: i
    real(amrex_real) :: dhx

    dhx = beta*dxinv(1)*dxinv(1)

    do i = lo(1), hi(1)
       y(i) = alpha*a(i)*x(i) &
            - dhx * (bX(i+1)*(x(i+1) - x(i  ))  &
            &      - bX(i  )*(x(i  ) - x(i-1)))
    end do
  end subroutine amrex_mlabeclap_adotx

end module amrex_mlabeclap_2d_module
