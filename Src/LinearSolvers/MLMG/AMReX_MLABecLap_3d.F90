
module amrex_mlabeclap_3d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlabeclap_adotx, amrex_mlabeclap_normalize, amrex_mlabeclap_flux

contains

  subroutine amrex_mlabeclap_adotx (lo, hi, y, ylo, yhi, x, xlo, xhi, a, alo, ahi, &
       bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, dxinv, alpha, beta) &
       bind(c,name='amrex_mlabeclap_adotx')
    integer, dimension(3), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, alo, ahi, bxlo, bxhi, &
         bylo, byhi, bzlo, bzhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), value, intent(in) :: alpha, beta
    real(amrex_real), intent(inout) ::  y( ylo(1): yhi(1), ylo(2): yhi(2), ylo(3): yhi(3))
    real(amrex_real), intent(in   ) ::  x( xlo(1): xhi(1), xlo(2): xhi(2), xlo(3): xhi(3))
    real(amrex_real), intent(in   ) ::  a( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3))
    real(amrex_real), intent(in   ) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
    real(amrex_real), intent(in   ) :: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
    real(amrex_real), intent(in   ) :: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
    
    integer :: i,j,k
    real(amrex_real) :: dhx, dhy, dhz

    dhx = beta*dxinv(1)*dxinv(1)
    dhy = beta*dxinv(2)*dxinv(2)
    dhz = beta*dxinv(3)*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             y(i,j,k) = alpha*a(i,j,k)*x(i,j,k) &
                  - dhx * (bX(i+1,j,k)*(x(i+1,j,k) - x(i  ,j,k))  &
                  &      - bX(i  ,j,k)*(x(i  ,j,k) - x(i-1,j,k))) &
                  - dhy * (bY(i,j+1,k)*(x(i,j+1,k) - x(i,j  ,k))  &
                  &      - bY(i,j  ,k)*(x(i,j  ,k) - x(i,j-1,k))) &
                  - dhz * (bZ(i,j,k+1)*(x(i,j,k+1) - x(i,j,k  ))  &
                  &      - bZ(i,j,k  )*(x(i,j,k  ) - x(i,j,k-1)))
          end do
       end do
    end do
  end subroutine amrex_mlabeclap_adotx


  subroutine amrex_mlabeclap_normalize (lo, hi, x, xlo, xhi, a, alo, ahi, &
       bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, dxinv, alpha, beta) &
       bind(c,name='amrex_mlabeclap_normalize')
    integer, dimension(3), intent(in) :: lo, hi, xlo, xhi, alo, ahi, bxlo, bxhi, &
         bylo, byhi, bzlo, bzhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), value, intent(in) :: alpha, beta
    real(amrex_real), intent(inout) ::  x( xlo(1): xhi(1), xlo(2): xhi(2), xlo(3): xhi(3))
    real(amrex_real), intent(in   ) ::  a( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3))
    real(amrex_real), intent(in   ) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
    real(amrex_real), intent(in   ) :: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
    real(amrex_real), intent(in   ) :: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
    
    integer :: i,j,k
    real(amrex_real) :: dhx, dhy, dhz

    dhx = beta*dxinv(1)*dxinv(1)
    dhy = beta*dxinv(2)*dxinv(2)
    dhz = beta*dxinv(3)*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             x(i,j,k) = x(i,j,k) / &
                  (alpha*a(i,j,k) + dhx*(bX(i,j,k)+bX(i+1,j,k)) &
                  &               + dhy*(bY(i,j,k)+bY(i,j+1,k)) &
                  &               + dhz*(bZ(i,j,k)+bZ(i,j,k+1)))
          end do
       end do
    end do
  end subroutine amrex_mlabeclap_normalize

  
  subroutine amrex_mlabeclap_flux (lo, hi, fx, fxlo, fxhi, fy, fylo, fyhi, &
       fz, fzlo, fzhi, sol, slo, shi, bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, &
       dxinv, beta, face_only) &
       bind(c, name='amrex_mlabeclap_flux')
    integer, dimension(3), intent(in) :: lo, hi, fxlo, fxhi, fylo, fyhi, fzlo, fzhi, &
         slo, shi, bxlo, bxhi, bylo, byhi, bzlo, bzhi
    real(amrex_real) :: dxinv(3)
    real(amrex_real), value, intent(in) :: beta
    integer, value, intent(in) :: face_only
    real(amrex_real), intent(inout) :: fx (fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
    real(amrex_real), intent(inout) :: fy (fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
    real(amrex_real), intent(inout) :: fz (fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
    real(amrex_real), intent(in   ) :: sol( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))    
    real(amrex_real), intent(in   ) :: bx (bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
    real(amrex_real), intent(in   ) :: by (bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
    real(amrex_real), intent(in   ) :: bz (bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
    
    integer :: i,j,k
    real(amrex_real) :: dhx, dhy, dhz

    dhx = beta*dxinv(1)
    dhy = beta*dxinv(2)
    dhz = beta*dxinv(3)

    if (face_only .eq. 1) then
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)+1, hi(1)+1-lo(1)
                fx(i,j,k) = -dhx * bx(i,j,k)*(sol(i,j,k) - sol(i-1,j,k))
             end do
          end do
       end do
       
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)+1, hi(2)+1-lo(2)
             do i = lo(1), hi(1)
                fy(i,j,k) = -dhy * by(i,j,k)*(sol(i,j,k) - sol(i,j-1,k))
             end do
          end do
       end do
       
       do       k = lo(3), hi(3)+1, hi(3)+1-lo(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                fz(i,j,k) = -dhz * bz(i,j,k)*(sol(i,j,k) - sol(i,j,k-1))
             end do
          end do
       end do

    else

       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)+1
                fx(i,j,k) = -dhx * bx(i,j,k)*(sol(i,j,k) - sol(i-1,j,k))
             end do
          end do
       end do
       
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)+1
             do i = lo(1), hi(1)
                fy(i,j,k) = -dhy * by(i,j,k)*(sol(i,j,k) - sol(i,j-1,k))
             end do
          end do
       end do
       
       do       k = lo(3), hi(3)+1
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                fz(i,j,k) = -dhz * bz(i,j,k)*(sol(i,j,k) - sol(i,j,k-1))
             end do
          end do
       end do
    end if

  end subroutine amrex_mlabeclap_flux

end module amrex_mlabeclap_3d_module
