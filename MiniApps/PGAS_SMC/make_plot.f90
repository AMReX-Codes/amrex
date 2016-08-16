module make_plot_module

  use iso_c_binding
  use variables_module
  use chemistry_module, only : nspecies

  implicit none

  integer, parameter :: pv_rho = 1
  integer, parameter :: pv_vx  = 2
  integer, parameter :: pv_vy  = 3
  integer, parameter :: pv_vz  = 4
  integer, parameter :: pv_T   = 5
  integer, parameter :: pv_Y   = 6

  private

  public :: make_plot_3d

contains

  subroutine make_plot_3d(lo, hi, pv, plo, phi, u, ulo, uhi) bind(c,name='make_plot_3d')
    integer, intent(in) :: lo(3), hi(3), plo(3), phi(3), ulo(3), uhi(3)
    double precision, intent(in)    ::  u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),ncons)
    double precision, intent(inout) :: pv(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),nspecies+5)

    integer :: i, j, k, n, iwrk, ierr
    double precision :: rho, rhoinv, rwrk, Y(nspecies), ei, Tt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rho = u(i,j,k,irho)
             rhoinv = 1.d0/rho
             pv(i,j,k,pv_rho) = rho
             pv(i,j,k,pv_vx ) = u(i,j,k,imx) * rhoinv
             pv(i,j,k,pv_vy ) = u(i,j,k,imy) * rhoinv
             pv(i,j,k,pv_vz ) = u(i,j,k,imz) * rhoinv
             
             do n = 1, nspecies
                Y(n) = u(i,j,k,iry1+n-1) * rhoinv
                pv(i,j,k,pv_Y+n-1) = Y(n)
             end do

             ei = rhoinv*u(i,j,k,iene) - 0.5d0*(pv(i,j,k,pv_vx)**2+pv(i,j,k,pv_vy)**2+pv(i,j,k,pv_vz)**2)
             Tt = 1.d3
             call get_T_given_ey(ei, Y, iwrk, rwrk, Tt, ierr)
             pv(i,j,k,pv_T) = Tt

          end do
       end do
    end do
  end subroutine make_plot_3d

end module make_plot_module
