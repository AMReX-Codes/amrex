
subroutine get_face_velocity(level, time, &
     vx, vx_l1, vx_l2, vx_l3, vx_h1, vx_h2, vx_h3, &
     vy, vy_l1, vy_l2, vy_l3, vy_h1, vy_h2, vy_h3, &
     vz, vz_l1, vz_l2, vz_l3, vz_h1, vz_h2, vz_h3, &
     dx, prob_lo) bind(C, name="get_face_velocity")

  use probdata_module, only : adv_vel

  implicit none

  integer, intent(in) :: level
  double precision, intent(in) :: time
  integer, intent(in) :: vx_l1, vx_l2, vx_l3, vx_h1, vx_h2, vx_h3
  integer, intent(in) :: vy_l1, vy_l2, vy_l3, vy_h1, vy_h2, vy_h3
  integer, intent(in) :: vz_l1, vz_l2, vz_l3, vz_h1, vz_h2, vz_h3
  double precision, intent(out) :: vx(vx_l1:vx_h1,vx_l2:vx_h2,vx_l3:vx_h3)
  double precision, intent(out) :: vy(vy_l1:vy_h1,vy_l2:vy_h2,vy_l3:vy_h3)
  double precision, intent(out) :: vz(vz_l1:vz_h1,vz_l2:vz_h2,vz_l3:vz_h3)
  double precision, intent(in) :: dx(3), prob_lo(3)

  vx = adv_vel(1)
  vy = adv_vel(2)
  vz = adv_vel(3)

end subroutine get_face_velocity
