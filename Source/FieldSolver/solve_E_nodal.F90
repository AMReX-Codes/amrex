module warpx_ES_solve_E_nodal

  use iso_c_binding
  use amrex_fort_module, only : amrex_real

  implicit none

contains

! This routine computes the node-centered electric field given a node-centered phi.
! The gradient is computed using 2nd-order centered differences. It assumes the
! Boundary conditions have already been set and that you have two rows of ghost cells
! for phi and one row of ghost cells for Ex, Ey, and Ez.
! Note that this routine includes the minus sign in E = - grad phi.
!
! Arguments:
!     lo, hi:     The corners of the valid box over which the gradient is taken
!     Ex, Ey, Ez: The electric field in the x, y, and z directions.
!     dx:         The cell spacing
!
  subroutine warpx_compute_E_nodal_3d (lo, hi, phi, Ex, Ey, Ez, dx) &
       bind(c,name='warpx_compute_E_nodal_3d')
    integer(c_int),   intent(in)    :: lo(3), hi(3)
    real(amrex_real), intent(in)    :: dx(3)
    real(amrex_real), intent(in   ) :: phi(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2)
    real(amrex_real), intent(inout) :: Ex (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
    real(amrex_real), intent(inout) :: Ey (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
    real(amrex_real), intent(inout) :: Ez (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

    integer :: i, j, k
    real(amrex_real) :: fac(3)

    fac = 0.5d0 / dx

    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             Ex(i,j,k) = fac(1) * (phi(i-1,j,k) - phi(i+1,j,k))
             Ey(i,j,k) = fac(2) * (phi(i,j-1,k) - phi(i,j+1,k))
             Ez(i,j,k) = fac(3) * (phi(i,j,k-1) - phi(i,j,k+1))

          end do
       end do
    end do

  end subroutine warpx_compute_E_nodal_3d

  subroutine warpx_compute_E_nodal_2d (lo, hi, phi, Ex, Ey, dx) &
       bind(c,name='warpx_compute_E_nodal_2d')
    integer(c_int),   intent(in)    :: lo(2), hi(2)
    real(amrex_real), intent(in)    :: dx(2)
    real(amrex_real), intent(in   ) :: phi(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2)
    real(amrex_real), intent(inout) :: Ex (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)
    real(amrex_real), intent(inout) :: Ey (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)

    integer :: i, j
    real(amrex_real) :: fac(2)

    fac = 0.5d0 / dx

    do j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)+1

          Ex(i,j) = fac(1) * (phi(i-1,j) - phi(i+1,j))
          Ey(i,j) = fac(2) * (phi(i,j-1) - phi(i,j+1))

       end do
    end do

  end subroutine warpx_compute_E_nodal_2d

end module warpx_ES_solve_E_nodal
