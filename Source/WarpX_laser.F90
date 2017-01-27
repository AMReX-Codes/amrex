
module warpx_laser_module

  use iso_c_binding
  use bl_fort_module, only : c_real

  implicit none

contains

  subroutine warpx_laser_pusher(np,xp,yp,zp,uxp,uyp,uzp,gaminv,&
                                q,dt,laser_pusher_algo) &
       bind(C, name="warpx_laser_pusher")

    INTEGER(c_long), INTENT(IN)   :: np
    REAL(c_real),INTENT(INOUT)    :: gaminv(np)
    REAL(c_real),INTENT(INOUT)    :: xp(np),yp(np),zp(np)
    REAL(c_real),INTENT(INOUT)    :: uxp(np),uyp(np),uzp(np)
    REAL(c_real),INTENT(IN)       :: q,dt
    INTEGER(c_long), INTENT(IN)   :: laser_pusher_algo

    ! need to update velocity!

#if (BL_SPACEDIM == 3)    
    CALL pxr_pushxyz(np,xp,yp,zp,uxp,uyp,uzp,gaminv,dt)
#elif (BL_SPACEDIM == 2)
    CALL pxr_pushxz(np,xp,zp,uxp,uzp,gaminv,dt)
#endif

  end subroutine warpx_laser_pusher

end module warpx_laser_module
