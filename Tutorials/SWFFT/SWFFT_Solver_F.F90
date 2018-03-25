#include "AMReX_LO_BCTYPES.H"

module abl_module

  use amrex_fort_module
  use amrex_error_module
  implicit none

contains

  subroutine fort_init_rhs (lo, hi, rhs, rlo, rhi, dx) &
       bind(c,name="fort_init_rhs")
    integer, intent(in) :: lo(3), hi(3), rlo(3), rhi(3)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(in) ::  dx(3)

    integer :: i,j,k
    real(amrex_real) :: x,y,z
    real(amrex_real) :: pi, fpi, tpi, fac

    pi = 4.d0 * atan(1.d0)
    tpi = 2.0d0 * pi
    fpi = 4.0d0 * pi
    fac = 3.0d0 * tpi**2

    do k = lo(3), hi(3)
       z = (dble(k)+0.5d0)*dx(3)
       
       do j = lo(2), hi(2)
          y = (dble(j)+0.5d0)*dx(2)
          
          do i = lo(1), hi(1)
             x = (dble(i)+0.5d0)*dx(1)
             
             rhs(i,j,k) = -fac * (sin(tpi*x) * sin(tpi*y) * sin(tpi*z))  &
                  &       -fac * (sin(fpi*x) * sin(fpi*y) * sin(fpi*z))
          end do
       end do
    end do

  end subroutine fort_init_rhs

  subroutine fort_comp_asol (soln, slo, shi, dx) bind(c,name='fort_comp_asol')
    integer, intent(in) :: slo(3), shi(3)
    real(amrex_real), intent(inout) :: soln(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(amrex_real), intent(in) :: dx(3)

    integer :: i,j,k
    real(amrex_real) :: x,y,z
    real(amrex_real) :: pi, fpi, tpi

    pi = 4.d0 * atan(1.d0)
    tpi = 2.0d0 * pi
    fpi = 4.0d0 * pi

    do k = slo(3), shi(3)
       z = (dble(k)+0.5d0)*dx(3)
       do j = slo(2), shi(2)
          y = (dble(j)+0.5d0)*dx(2)
          do i = slo(1), shi(1)
             x = (dble(i)+0.5d0)*dx(1)
             
             soln(i,j,k) = 1.d0 * (sin(tpi*x) * sin(tpi*y) * sin(tpi*z))  &
                 &      + .25d0 * (sin(fpi*x) * sin(fpi*y) * sin(fpi*z))

          end do
       end do
    end do
    
  end subroutine fort_comp_asol

end module abl_module
