
module amrex_mlabeclap_1d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlabeclap_adotx, amrex_mlabeclap_normalize, amrex_mlabeclap_flux

contains

  subroutine amrex_mlabeclap_adotx (lo, hi, y, ylo, yhi, x, xlo, xhi, a, alo, ahi, &
       bx, bxlo, bxhi, dxinv, alpha, beta) bind(c,name='amrex_mlabeclap_adotx')
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


  subroutine amrex_mlabeclap_normalize (lo, hi, x, xlo, xhi, a, alo, ahi, &
       bx, bxlo, bxhi, dxinv, alpha, beta) bind(c,name='amrex_mlabeclap_normalize')
    integer, dimension(1), intent(in) :: lo, hi, xlo, xhi, alo, ahi, bxlo, bxhi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), value, intent(in) :: alpha, beta
    real(amrex_real), intent(inout) ::  x( xlo(1): xhi(1))
    real(amrex_real), intent(in   ) ::  a( alo(1): ahi(1))
    real(amrex_real), intent(in   ) :: bx(bxlo(1):bxhi(1))
    
    integer :: i
    real(amrex_real) :: dhx

    dhx = beta*dxinv(1)*dxinv(1)

    do i = lo(1), hi(1)
       x(i) = x(i) / (alpha*a(i) + dhx*(bX(i)+bX(i+1)))
    end do
  end subroutine amrex_mlabeclap_normalize


  subroutine amrex_mlabeclap_flux (lo, hi, fx, fxlo, fxhi, sol, slo, shi, bx, bxlo, bxhi, &
       dxinv, beta, face_only) bind(c, name='amrex_mlabeclap_flux')
    integer, dimension(1), intent(in) :: lo, hi, fxlo, fxhi, slo, shi, bxlo, bxhi
    real(amrex_real) :: dxinv(1)
    real(amrex_real), value, intent(in) :: beta
    integer, value, intent(in) :: face_only
    real(amrex_real), intent(inout) :: fx (fxlo(1):fxhi(1))
    real(amrex_real), intent(in   ) :: sol( slo(1): shi(1))
    real(amrex_real), intent(in   ) :: bx (bxlo(1):bxhi(1))
    
    integer :: i
    real(amrex_real) :: dhx

    dhx = beta*dxinv(1)

    if (face_only .eq. 1) then
       do i = lo(1), hi(1)+1, hi(1)+1-lo(1)
          fx(i) = -dhx * bx(i)*(sol(i) - sol(i-1))
       end do       
    else
       do i = lo(1), hi(1)+1
          fx(i) = -dhx * bx(i)*(sol(i) - sol(i-1))
       end do
    end if
  end subroutine amrex_mlabeclap_flux

end module amrex_mlabeclap_1d_module
