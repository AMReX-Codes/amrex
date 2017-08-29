
module warpx_boosted_frame_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use constants

  implicit none

contains

!
! Given cell-centered data in the boosted reference fraem of the simulation, 
! this transforms E and B in place so that the multifab now contains values
! in the lab frame. This routine assumes that the simulation frame is moving
! in the positive z direction with respect to the lab frame.
! 
  subroutine warpx_lorentz_transform_3d(data, dlo, dhi, tlo, thi, gamma_boost, beta_boost) &
       bind(C, name="warpx_lorentz_transform_3d")

    integer(c_int),   intent(in)    :: dlo(3), dhi(3)
    integer(c_int),   intent(in)    :: tlo(3), thi(3)
    real(amrex_real), intent(inout) :: data(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3), 6)
    real(amrex_real), intent(in)    :: gamma_boost, beta_boost

    integer i, j, k
    real(amrex_real) ex_lab, ey_lab, bx_lab, by_lab
    
    do k = tlo(3), thi(3)
       do j = tlo(2), thi(2)
          do i = tlo(1), thi(1)
             
             ! note that ez and bz are not changed by the tranform
             ex_lab = gamma_boost * (data(i, j, k, 1) + beta_boost*clight*data(i, j, k, 5))
             by_lab = gamma_boost * (data(i, j, k, 5) + beta_boost*data(i, j, k, 1)/clight)

             data(i, j, k, 1) = ex_lab
             data(i, j, k, 5) = by_lab

             ey_lab = gamma_boost * (data(i, j, k, 2) - beta_boost*clight*data(i, j, k, 4))
             bx_lab = gamma_boost * (data(i, j, k, 4) - beta_boost*data(i, j, k, 2)/clight)

             data(i, j, k, 2) = ey_lab
             data(i, j, k, 4) = bx_lab

          end do
       end do
    end do

  end subroutine warpx_lorentz_transform_3d

end module warpx_boosted_frame_module
