module cns_nd_module

  use amrex_fort_module, only : rt=>amrex_real
  use cns_module, only : urho, umx, umy, umz, ueden, ueint, utemp, nvar
  implicit none
  private

  public :: cns_compute_temperature, cns_estdt

contains

  subroutine cns_compute_temperature (lo,hi,u,ulo,uhi) bind(c,name='cns_compute_temperature')
    use cns_physics_module, only : cv
    integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3)
    real(rt), intent(inout) :: u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),nvar)

    integer :: i,j,k
    real(rt) :: rhoinv

    do       k = lo(3),hi(3)
       do    j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhoinv = 1.d0/u(i,j,k,urho)
             u(i,j,k,ueint) = u(i,j,k,ueden) - 0.5d0*rhoinv* &
                  (u(i,j,k,umx)**2 + u(i,j,k,umy)**2 + u(i,j,k,umz)**2)
             u(i,j,k,utemp) = rhoinv*u(i,j,k,ueint)*(1.d0/cv)
          end do
       end do
    end do
  end subroutine cns_compute_temperature

  subroutine cns_estdt (lo,hi,u,ulo,uhi,dx,dt) bind(c,name='cns_estdt')
    use cns_physics_module, only : gamma
    use cns_module, only : smallp, smallr
    integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3)
    real(rt), intent(in) :: u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),nvar)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: dt

    integer :: i,j,k
    real(rt) :: rhoinv, vx, vy, vz, p, cs

    do       k = lo(3),hi(3)
       do    j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhoinv = 1.d0/max(smallr,u(i,j,k,urho))
             vx = u(i,j,k,umx)*rhoinv
             vy = u(i,j,k,umy)*rhoinv
             vz = u(i,j,k,umz)*rhoinv
             p = max((gamma-1.d0)*u(i,j,k,ueint),smallp)
             cs = sqrt(gamma*p*rhoinv)
             dt = min(dt,dx(1)/(abs(vx)+cs),dx(2)/(abs(vy)+cs),dx(3)/(abs(vz)+cs))
          end do
       end do
    end do
  end subroutine cns_estdt

  subroutine cns_compute_dsdt (lo,hi,ut,ulo,uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx) &
       bind(c,name='cns_compute_dsdt')
    integer, intent(in) :: lo(3),hi(3),ulo(3),uhi(3),fxlo(3),fxhi(3),fylo(3),fyhi(3),fzlo(3),fzhi(3)
    real(rt), intent(inout) :: ut( ulo(1): uhi(1), ulo(2): uhi(2), ulo(3): uhi(3),nvar)
    real(rt), intent(in   ) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),nvar)
    real(rt), intent(in   ) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),nvar)
    real(rt), intent(in   ) :: fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),nvar)
    real(rt), intent(in) :: dx(3)

    integer :: i,j,k,n
    real(rt) :: dxinv(3)

    dxinv = 1.d0/dx

    do n = 1, nvar
       do       k = lo(3),hi(3)
          do    j = lo(2),hi(2)
             do i = lo(1),hi(1)
                ut(i,j,k,n) = ut(i,j,k,n) &
                     + (fx(i,j,k,n)-fx(i+1,j,k,n)) * dxinv(1) &
                     + (fy(i,j,k,n)-fy(i,j+1,k,n)) * dxinv(2) &
                     + (fz(i,j,k,n)-fz(i,j,k+1,n)) * dxinv(3)
             end do
          end do
       end do
    end do
  end subroutine cns_compute_dsdt

end module cns_nd_module
