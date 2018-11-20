#include "AMReX_LO_BCTYPES.H"

module abl_module

  use amrex_fort_module
  use amrex_error_module
  implicit none

contains

#if (AMREX_SPACEDIM == 2)

  subroutine fort_init_phi_spatial (lo, hi, phi_spatial, rlo, rhi, domlo, domhi, problo, probhi, dx, prob_type) &
       bind(c,name="fort_init_phi_spatial")

    integer, intent(in) :: prob_type
    integer, intent(in) :: lo(2), hi(2), rlo(2), rhi(2), domlo(2), domhi(2)
    real(amrex_real), intent(inout) :: phi_spatial(rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in) :: problo(2), probhi(2)
    real(amrex_real), intent(in) ::  dx(2)

    integer :: i,j
    real(amrex_real) :: x,y
    real(amrex_real) :: Lx,Ly
    real(amrex_real) :: pi, fpi, tpi, fac

    Lx = probhi(1)-problo(1)
    Ly = probhi(2)-problo(2)

    pi = 4.d0 * atan(1.d0)
    tpi = 2.0d0 * pi
    fpi = 4.0d0 * pi
    fac = 1.d0

    do j = lo(2), hi(2)
       y = (dble(j)+0.5d0)*dx(2)/Ly

       do i = lo(1), hi(1)
          x = (dble(i)+0.5d0)*dx(1)/Lx

          SELECT CASE (prob_type)
          CASE (0)
             phi_spatial(i,j) = -fac * sin(tpi*x)
          CASE (1)
             phi_spatial(i,j) = -fac * (sin(tpi*x) * sin(tpi*y))  &
                  &       -fac * (sin(fpi*x) * sin(fpi*y))
          CASE (2)
             phi_spatial(i,j) = -fac * (sin(1.d0*tpi*x) * sin(1.d0*tpi*y))  &
                  &       -fac * (sin(2.d0*tpi*x) * sin(2.d0*tpi*y))  &
                  &       -fac * (sin(3.d0*tpi*x) * sin(3.d0*tpi*y))  &
                  &       -fac * (sin(4.d0*tpi*x) * sin(4.d0*tpi*y))
          CASE (3)
             phi_spatial = 0.d0
             if (SUM(lo-domlo) == 0) then
                phi_spatial(lo(1),lo(2)) = 1.d0
             endif
          CASE (4)
             phi_spatial(i,j) = exp(-5.d2*((x/Lx-0.5d0)**2 + (y/Ly-0.5d0)**2))
          CASE DEFAULT
             phi_spatial = 0.d0
          END SELECT

       end do
    end do

  end subroutine fort_init_phi_spatial

#elif (AMREX_SPACEDIM == 3)

  subroutine fort_init_phi_spatial (lo, hi, phi_spatial, rlo, rhi, domlo, domhi, problo, probhi, dx, prob_type) &
       bind(c,name="fort_init_phi_spatial")

    integer, intent(in) :: prob_type
    integer, intent(in) :: lo(3), hi(3), rlo(3), rhi(3), domlo(3), domhi(3)
    real(amrex_real), intent(inout) :: phi_spatial(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(in) :: problo(3), probhi(3)
    real(amrex_real), intent(in) ::  dx(3)

    integer :: i,j,k
    real(amrex_real) :: x,y,z
    real(amrex_real) :: Lx,Ly,Lz
    real(amrex_real) :: pi, fpi, tpi, fac

    Lx = probhi(1)-problo(1)
    Ly = probhi(2)-problo(2)
    Lz = probhi(3)-problo(3)

    pi = 4.d0 * atan(1.d0)
    tpi = 2.0d0 * pi
    fpi = 4.0d0 * pi
    ! fac = 3.0d0 * tpi**2 / (Lx**2 * Ly**2 * Lz**2)
    fac = 1.0d0

    do k = lo(3), hi(3)
       z = (dble(k)+0.5d0)*dx(3)/Lz

       do j = lo(2), hi(2)
          y = (dble(j)+0.5d0)*dx(2)/Ly

          do i = lo(1), hi(1)
             x = (dble(i)+0.5d0)*dx(1)/Lx

             SELECT CASE (prob_type)
             CASE (0)
                phi_spatial(i,j,k) = -fac * sin(tpi*x)
             CASE (1)
                phi_spatial(i,j,k) = -fac * (sin(tpi*x) * sin(tpi*y) * sin(tpi*z))  &
                     &       -fac * (sin(fpi*x) * sin(fpi*y) * sin(fpi*z))
             CASE (2)
                phi_spatial(i,j,k) = -fac * (sin(1.d0*tpi*x) * sin(1.d0*tpi*y) * sin(1.d0*tpi*z))  &
                     &       -fac * (sin(2.d0*tpi*x) * sin(2.d0*tpi*y) * sin(2.d0*tpi*z))  &
                     &       -fac * (sin(3.d0*tpi*x) * sin(3.d0*tpi*y) * sin(3.d0*tpi*z))  &
                     &       -fac * (sin(4.d0*tpi*x) * sin(4.d0*tpi*y) * sin(4.d0*tpi*z))
             CASE (3)
                phi_spatial = 0.d0
                if (SUM(lo-domlo) == 0) then
                   phi_spatial(lo(1),lo(2),lo(3)) = 1.d0
                endif
             CASE (4)
                phi_spatial(i,j,k) = exp(-5.d2*((x/Lx-0.5d0)**2 + (y/Ly-0.5d0)**2 + (z/Lz-0.5d0)**2))
             CASE DEFAULT
                phi_spatial = 0.d0
             END SELECT

          end do
       end do
    end do

  end subroutine fort_init_phi_spatial

#endif

end module abl_module
