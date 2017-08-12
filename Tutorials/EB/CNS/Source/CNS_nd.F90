module cns_nd_module

  use amrex_fort_module, only : rt=>amrex_real
  use cns_module, only : urho, umx, umy, umz, ueden, ueint, utemp, nvar
  implicit none
  private

  public :: cns_compute_temperature

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

end module cns_nd_module
