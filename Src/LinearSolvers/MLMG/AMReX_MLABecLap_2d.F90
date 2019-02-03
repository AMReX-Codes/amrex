
module amrex_mlabeclap_2d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlabeclap_flux

contains

  subroutine amrex_mlabeclap_flux (lo, hi, fx, fxlo, fxhi, fy, fylo, fyhi, &
       sol, slo, shi, bx, bxlo, bxhi, by, bylo, byhi, dxinv, beta, face_only) &
       bind(c, name='amrex_mlabeclap_flux')
    integer, dimension(2), intent(in) :: lo, hi, fxlo, fxhi, fylo, fyhi, &
         slo, shi, bxlo, bxhi, bylo, byhi
    real(amrex_real) :: dxinv(2)
    real(amrex_real), value, intent(in) :: beta
    integer, value, intent(in) :: face_only
    real(amrex_real), intent(inout) :: fx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    real(amrex_real), intent(inout) :: fy (fylo(1):fyhi(1),fylo(2):fyhi(2))
    real(amrex_real), intent(in   ) :: sol( slo(1): shi(1), slo(2): shi(2))    
    real(amrex_real), intent(in   ) :: bx (bxlo(1):bxhi(1),bxlo(2):bxhi(2))
    real(amrex_real), intent(in   ) :: by (bylo(1):byhi(1),bylo(2):byhi(2))
    
    integer :: i,j
    real(amrex_real) :: dhx, dhy

    dhx = beta*dxinv(1)
    dhy = beta*dxinv(2)

    if (face_only .eq. 1) then
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)+1, hi(1)+1-lo(1)
             fx(i,j) = -dhx * bx(i,j)*(sol(i,j) - sol(i-1,j))
          end do
       end do
       
       do    j = lo(2), hi(2)+1, hi(2)+1-lo(2)
          do i = lo(1), hi(1)
             fy(i,j) = -dhy * by(i,j)*(sol(i,j) - sol(i,j-1))
          end do
       end do
       
    else

       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             fx(i,j) = -dhx * bx(i,j)*(sol(i,j) - sol(i-1,j))
          end do
       end do
       
       do    j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             fy(i,j) = -dhy * by(i,j)*(sol(i,j) - sol(i,j-1))
          end do
       end do
       
    end if

  end subroutine amrex_mlabeclap_flux

end module amrex_mlabeclap_2d_module
