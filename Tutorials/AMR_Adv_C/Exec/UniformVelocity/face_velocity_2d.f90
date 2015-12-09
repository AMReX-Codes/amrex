
subroutine get_face_velocity(level, time, &
     vx, vx_l1, vx_l2, vx_h1, vx_h2, &
     vy, vy_l1, vy_l2, vy_h1, vy_h2, &
     dx, prob_lo)

  use probdata_module, only : adv_vel

  implicit none

  integer, intent(in) :: level
  double precision, intent(in) :: time
  integer, intent(in) :: vx_l1, vx_l2, vx_h1, vx_h2
  integer, intent(in) :: vy_l1, vy_l2, vy_h1, vy_h2
  double precision, intent(out) :: vx(vx_l1:vx_h1,vx_l2:vx_h2)
  double precision, intent(out) :: vy(vy_l1:vy_h1,vy_l2:vy_h2)
  double precision, intent(in) :: dx(2), prob_lo(2)

  vx = adv_vel(1)
  vy = adv_vel(2)

end subroutine get_face_velocity
