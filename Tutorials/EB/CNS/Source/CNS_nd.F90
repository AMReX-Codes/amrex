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
    integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3)
    real(rt), intent(in) :: u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),nvar)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: dt

    integer :: i,j,k
    real(rt) :: rhoinv, vx, vy, vz, p, cs

    do       k = lo(3),hi(3)
       do    j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhoinv = 1.d0/u(i,j,k,urho)
             vx = u(i,j,k,umx)*rhoinv
             vy = u(i,j,k,umy)*rhoinv
             vz = u(i,j,k,umz)*rhoinv
             p = (gamma-1.d0)*u(i,j,k,ueint)
             cs = sqrt(gamma*p*rhoinv)
             dt = min(dt,dx(1)/(abs(vx)+cs),dx(2)/(abs(vy)+cs),dx(3)/(abs(vz)+cs))
          end do
       end do
    end do
  end subroutine cns_estdt

end module cns_nd_module
