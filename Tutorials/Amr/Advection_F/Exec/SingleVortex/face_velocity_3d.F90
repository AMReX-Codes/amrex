module face_velocity_module

  implicit none
  private

  public :: get_face_velocity

contains

subroutine get_face_velocity(time, &
     vx, vxlo, vxhi, &
     vy, vylo, vyhi, &
     vz, vzlo, vzhi, &
     dx, prob_lo) bind(C, name="get_face_velocity")

use amrex_base_module

  implicit none

  !integer, intent(in) :: level
  double precision, intent(in) :: time
  integer, intent(in) :: vxlo(3),vxhi(3)
  integer, intent(in) :: vylo(3),vyhi(3)
  integer, intent(in) :: vzlo(3),vzhi(3)
  double precision, intent(out) :: vx(vxlo(1):vxhi(1),vxlo(2):vxhi(2),vxlo(3):vxhi(3))
  double precision, intent(out) :: vy(vylo(1):vyhi(1),vylo(2):vyhi(2),vylo(3):vyhi(3))
  double precision, intent(out) :: vz(vzlo(1):vzhi(1),vzlo(2):vzhi(2),vzlo(3):vzhi(3))
  double precision, intent(in) :: dx(3), prob_lo(3)

  integer :: i, j, k, plo(2), phi(2)
  double precision :: x, y, z
  double precision, pointer, contiguous :: psi(:,:)
  double precision, parameter :: M_PI = 3.141592653589793238462643383279502884197d0
  integer :: vx_l1,vx_l2,vx_l3,vy_l1,vy_l2,vy_l3,vz_l1,vz_l2,vz_l3
  integer :: vx_h1,vx_h2,vx_h3,vy_h1,vy_h2,vy_h3,vz_h1,vz_h2,vz_h3

  vx_l1=vxlo(1);vx_l2=vxlo(2);vx_l3=vxlo(3)
  vy_l1=vylo(1);vy_l2=vylo(2);vy_l3=vylo(3)
  vz_l1=vzlo(1);vz_l2=vzlo(2);vz_l3=vzlo(3) 
  vx_h1=vxhi(1);vx_h2=vxhi(2);vx_h3=vxhi(3)
  vy_h1=vyhi(1);vy_h2=vyhi(2);vy_h3=vyhi(3)
  vz_h1=vzhi(1);vz_h2=vzhi(2);vz_h3=vzhi(3)

  plo(1) = min(vx_l1-1, vy_l1-1)
  plo(2) = min(vx_l2-1, vy_l2-1)
  phi(1) = max(vx_h1  , vy_h1+1)
  phi(2) = max(vx_h2+1, vy_h2  )
  
  call bl_allocate(psi, plo(1), phi(1), plo(2), phi(2))

  ! streamfunction psi
  do j = plo(2), phi(2)
     y = (dble(j)+0.5d0)*dx(2) + prob_lo(2)
     do i = plo(1), phi(1)
        x = (dble(i)+0.5d0)*dx(1) + prob_lo(1)
        psi(i,j) =  sin(M_PI*x)**2 * sin(M_PI*y)**2 * cos (M_PI*time/2.d0) * (1.d0 / M_PI)
     end do
  end do
  
  ! x velocity
  do k = vx_l3, vx_h3
  do j = vx_l2, vx_h2
     y = (dble(j)+0.5d0) * dx(2) + prob_lo(2)
     do i = vx_l1, vx_h1
        x = dble(i) * dx(1) + prob_lo(1)
        vx(i,j,k) =  -( (psi(i,j+1)+psi(i-1,j+1)) - (psi(i,j-1)+psi(i-1,j-1)) ) * (0.25d0/dx(2))
     end do
  end do
  end do

  ! y velocity
  do k = vy_l3, vy_h3
  do j = vy_l2, vy_h2
     y = dble(j) * dx(2) + prob_lo(2)
     do i = vy_l1, vy_h1
        x = (dble(i)+0.5d0) * dx(1) + prob_lo(1)
        vy(i,j,k) = ( (psi(i+1,j)+psi(i+1,j-1)) - (psi(i-1,j)+psi(i-1,j-1)) ) * (0.25d0/dx(1))
     end do
  end do
  end do

  vz = 1.d0

  call bl_deallocate(psi)

end subroutine get_face_velocity

end module face_velocity_module
