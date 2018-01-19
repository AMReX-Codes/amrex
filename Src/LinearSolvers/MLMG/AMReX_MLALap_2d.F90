
module amrex_mlalap_2d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlalap_adotx, amrex_mlalap_flux, amrex_mlalap_gsrb

contains

  subroutine amrex_mlalap_adotx (lo, hi, y, ylo, yhi, x, xlo, xhi, a, alo, ahi, &
       rc, re, rlo, rhi, dxinv, alpha, beta) bind(c,name='amrex_mlalap_adotx')
    integer, dimension(2), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, alo, ahi
    integer, intent(in) :: rlo, rhi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), value, intent(in) :: alpha, beta
    real(amrex_real), intent(inout) ::  y( ylo(1): yhi(1), ylo(2): yhi(2))
    real(amrex_real), intent(in   ) ::  x( xlo(1): xhi(1), xlo(2): xhi(2))
    real(amrex_real), intent(in   ) ::  a( alo(1): ahi(1), alo(2): ahi(2))    
    real(amrex_real), intent(in) :: rc(rlo:rhi)
    real(amrex_real), intent(in) :: re(rlo:rhi+1)
  end subroutine amrex_mlalap_adotx


  subroutine amrex_mlalap_flux (lo, hi, fx, fxlo, fxhi, fy, fylo, fyhi, &
       sol, slo, shi, dxinv, beta, face_only) bind(c, name='amrex_mlalap_flux')
    integer, dimension(2), intent(in) :: lo, hi, fxlo, fxhi, fylo, fyhi, slo, shi
    real(amrex_real) :: dxinv(2)
    real(amrex_real), value, intent(in) :: beta
    integer, value, intent(in) :: face_only
    real(amrex_real), intent(inout) :: fx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    real(amrex_real), intent(inout) :: fy (fylo(1):fyhi(1),fylo(2):fyhi(2))
    real(amrex_real), intent(in   ) :: sol( slo(1): shi(1), slo(2): shi(2))    
  end subroutine amrex_mlalap_flux


  subroutine amrex_mlalap_gsrb (lo, hi, phi, hlo, hhi, rhs, rlo, rhi, &
       f0,f0lo,f0hi,f1,f1lo,f1hi,f2,f2lo,f2hi,f3,f3lo,f3hi, &
       m0,m0lo,m0hi,m1,m1lo,m1hi,m2,m2lo,m2hi,m3,m3lo,m3hi, &
       a, alo, ahi, rc, re, blo, bhi, dxinv, alpha, beta, redblack) &
       bind(c,name='amrex_mlalap_gsrb')
    integer, dimension(2), intent(in) :: lo, hi, hlo, hhi, rlo, rhi, &
         f0lo,f0hi,f1lo,f1hi,f2lo,f2hi,f3lo,f3hi, &
         m0lo,m0hi,m1lo,m1hi,m2lo,m2hi,m3lo,m3hi, alo,ahi,blo, bhi
    integer, intent(in), value :: redblack
    real(amrex_real), intent(in) :: dxinv(2), alpha, beta
    real(amrex_real), intent(inout) :: phi(hlo(1):hhi(1),hlo(2):hhi(2))
    real(amrex_real), intent(in   ) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) :: f0(f0lo(1):f0hi(1),f0lo(2):f0hi(2))
    real(amrex_real), intent(in   ) :: f1(f1lo(1):f1hi(1),f1lo(2):f1hi(2))
    real(amrex_real), intent(in   ) :: f2(f2lo(1):f2hi(1),f2lo(2):f2hi(2))
    real(amrex_real), intent(in   ) :: f3(f3lo(1):f3hi(1),f3lo(2):f3hi(2))
    real(amrex_real), intent(in   ) ::  a( alo(1): ahi(1), alo(2): ahi(2))
    integer         , intent(in   ) :: m0(m0lo(1):m0hi(1),m0lo(2):m0hi(2))
    integer         , intent(in   ) :: m1(m1lo(1):m1hi(1),m1lo(2):m1hi(2))
    integer         , intent(in   ) :: m2(m2lo(1):m2hi(1),m2lo(2):m2hi(2))
    integer         , intent(in   ) :: m3(m3lo(1):m3hi(1),m3lo(2):m3hi(2))
    real(amrex_real), intent(in) :: rc(blo(1):bhi(1))
    real(amrex_real), intent(in) :: re(blo(1):bhi(1)+1)
  end subroutine amrex_mlalap_gsrb

end module amrex_mlalap_2d_module
