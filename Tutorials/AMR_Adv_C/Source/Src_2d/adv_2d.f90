
subroutine advect(time, lo, hi, &
     &            phi_i,pi_l1,pi_l2,pi_h1,pi_h2, &
     &            phi_o,po_l1,po_l2,po_h1,po_h2, &
     &            vx,vx_l1,vx_l2,vx_h1,vx_h2, &
     &            vy,vy_l1,vy_l2,vy_h1,vy_h2, &
     &            fx,fx_l1,fx_l2,fx_h1,fx_h2, &
     &            fy,fy_l1,fy_l2,fy_h1,fy_h2, &
     &            dx,dt)
  
  use mempool_module, only : bl_allocate, bl_deallocate

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  double precision, intent(in) :: dx(2), dt, time
  integer, intent(in) :: pi_l1,pi_l2,pi_h1,pi_h2, &
       &                 po_l1,po_l2,po_h1,po_h2, &
       &                 vx_l1,vx_l2,vx_h1,vx_h2, &
       &                 vy_l1,vy_l2,vy_h1,vy_h2, &
       &                 fx_l1,fx_l2,fx_h1,fx_h2, &
       &                 fy_l1,fy_l2,fy_h1,fy_h2
  double precision, intent(in   ) :: phi_i(pi_l1:pi_h1,pi_l2:pi_h2)
  double precision, intent(inout) :: phi_o(po_l1:po_h1,po_l2:po_h2)
  double precision, intent(in   ) :: vx   (vx_l1:vx_h1,vx_l2:vx_h2)
  double precision, intent(in   ) :: vy   (vy_l1:vy_h1,vy_l2:vy_h2)
  double precision, intent(  out) :: fx   (fx_l1:fx_h1,fx_l2:fx_h2)
  double precision, intent(  out) :: fy   (fy_l1:fy_h1,fy_l2:fy_h2)


end subroutine advect
