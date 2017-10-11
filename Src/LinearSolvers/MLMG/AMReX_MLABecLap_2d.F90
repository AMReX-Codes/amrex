
module amrex_mlabeclap_2d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlabeclap_adotx

contains

  subroutine amrex_mlabeclap_adotx (lo, hi, y, ylo, yhi, x, xlo, xhi, a, alo, ahi, &
       bx, bxlo, bxhi, by, bylo, byhi, dxinv, alpha, beta) &
       bind(c,name='amrex_mlabeclap_adotx')
    integer, dimension(2), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, alo, ahi, bxlo, bxhi, &
         bylo, byhi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), value, intent(in) :: alpha, beta
    real(amrex_real), intent(inout) ::  y( ylo(1): yhi(1), ylo(2): yhi(2))
    real(amrex_real), intent(in   ) ::  x( xlo(1): xhi(1), xlo(2): xhi(2))
    real(amrex_real), intent(in   ) ::  a( alo(1): ahi(1), alo(2): ahi(2))
    real(amrex_real), intent(in   ) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2))
    real(amrex_real), intent(in   ) :: by(bylo(1):byhi(1),bylo(2):byhi(2))
    
    integer :: i,j
    real(amrex_real) :: dhx, dhy

    dhx = beta*dxinv(1)*dxinv(1)
    dhy = beta*dxinv(2)*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          y(i,j) = alpha*a(i,j)*x(i,j) &
               - dhx * (bX(i+1,j)*(x(i+1,j) - x(i  ,j))  &
               &      - bX(i  ,j)*(x(i  ,j) - x(i-1,j))) &
               - dhy * (bY(i,j+1)*(x(i,j+1) - x(i,j  ))  &
               &      - bY(i,j  )*(x(i,j  ) - x(i,j-1)))
       end do
    end do
  end subroutine amrex_mlabeclap_adotx

end module amrex_mlabeclap_2d_module
