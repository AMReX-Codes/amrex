
subroutine get_face_velocity(level, time, &
     vx, vx_l1, vx_l2, vx_h1, vx_h2, &
     vy, vy_l1, vy_l2, vy_h1, vy_h2, &
     dx, prob_lo)

  use amrex_mempool_module, only : bl_allocate, bl_deallocate

  implicit none

  integer, intent(in) :: level
  double precision, intent(in) :: time
  integer, intent(in) :: vx_l1, vx_l2, vx_h1, vx_h2
  integer, intent(in) :: vy_l1, vy_l2, vy_h1, vy_h2
  double precision, intent(out) :: vx(vx_l1:vx_h1,vx_l2:vx_h2)
  double precision, intent(out) :: vy(vy_l1:vy_h1,vy_l2:vy_h2)
  double precision, intent(in) :: dx(2), prob_lo(2)

  integer :: i, j, plo(2), phi(2)
  double precision :: x, y
  double precision, pointer, contiguous :: psi(:,:)
  double precision, parameter :: M_PI = 3.141592653589793238462643383279502884197d0

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
  do j = vx_l2, vx_h2
     y = (dble(j)+0.5d0) * dx(2) + prob_lo(2)
     do i = vx_l1, vx_h1
        x = dble(i) * dx(1) + prob_lo(1)
        vx(i,j) =  -( (psi(i,j+1)+psi(i-1,j+1)) - (psi(i,j-1)+psi(i-1,j-1)) ) * (0.25d0/dx(2))
     end do
  end do

  ! y velocity
  do j = vy_l2, vy_h2
     y = dble(j) * dx(2) + prob_lo(2)
     do i = vy_l1, vy_h1
        x = (dble(i)+0.5d0) * dx(1) + prob_lo(1)
        vy(i,j) = ( (psi(i+1,j)+psi(i+1,j-1)) - (psi(i-1,j)+psi(i-1,j-1)) ) * (0.25d0/dx(1))
     end do
  end do

  call bl_deallocate(psi)
  
end subroutine get_face_velocity

