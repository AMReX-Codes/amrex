module advect_module

  use amrex_base_module

  implicit none
  private

  public :: advect, advect_particles

contains

subroutine advect(time, lo, hi, &
     &            uin , ui_lo, ui_hi, &
     &            uout, uo_lo, uo_hi, &
     &            vx  , vx_lo, vx_hi, &
     &            vy  , vy_lo, vy_hi, &
     &            vz  , vz_lo, vz_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            flxz, fz_lo, fz_hi, &
     &            dx,dt) bind(C, name="advect")
  
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use compute_flux_module, only : compute_flux_3d

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  double precision, intent(in) :: dx(3), dt, time
  integer, intent(in) :: ui_lo(3), ui_hi(3)
  integer, intent(in) :: uo_lo(3), uo_hi(3)
  integer, intent(in) :: vx_lo(3), vx_hi(3)
  integer, intent(in) :: vy_lo(3), vy_hi(3)
  integer, intent(in) :: vz_lo(3), vz_hi(3)
  integer, intent(in) :: fx_lo(3), fx_hi(3)
  integer, intent(in) :: fy_lo(3), fy_hi(3)
  integer, intent(in) :: fz_lo(3), fz_hi(3)
  double precision, intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3))
  double precision, intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))
  double precision, intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
  double precision, intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2),vy_lo(3):vy_hi(3))
  double precision, intent(in   ) :: vz  (vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))
  double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
  double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
  double precision, intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))

  integer :: i, j, k
  integer :: glo(3), ghi(3)
  double precision :: dtdx(3), umax, vmax, wmax

  ! Some compiler may not support 'contiguous'.  Remove it in that case.
  double precision, dimension(:,:,:), pointer, contiguous :: &
       phix, phix_y, phix_z, phiy, phiy_x, phiy_z, phiz, phiz_x, phiz_y, slope

  dtdx = dt/dx

  glo = lo - 1
  ghi = hi + 1

  ! edge states
  call bl_allocate(phix  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phix_y,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phix_z,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiy  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiy_x,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiy_z,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiz  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiz_x,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiz_y,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  ! slope
  call bl_allocate(slope,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))  
  
  ! We like to allocate these **pointers** here and then pass them to a function
  ! to remove their pointerness for performance, because normally pointers could
  ! be aliasing.  We need to use pointers instead of allocatable arrays because
  ! we like to use BoxLib's bl_allocate to allocate memeory instead of the intrinsic
  ! allocate.  Bl_allocate is much faster than allocate inside OMP.  
  ! Note that one MUST CALL BL_DEALLOCATE.

  ! check if CFL condition is violated.
  umax = maxval(abs(vx))
  vmax = maxval(abs(vy))
  wmax = maxval(abs(vz))
  if ( umax*dt .ge. dx(1) .or. &
       vmax*dt .ge. dx(2) .or. &
       wmax*dt .ge. dx(3) ) then
     print *, "umax = ", umax, ", vmax = ", vmax, ", wmax = ", wmax, ", dt = ", dt, ", dx = ", dx
     call bl_error("CFL violation. Use smaller adv.cfl.")
  end if

  ! call a function to compute flux
  call compute_flux_3d(lo, hi, dt, dx, &
                       uin, ui_lo, ui_hi, &
                       vx, vx_lo, vx_hi, &
                       vy, vy_lo, vy_hi, &
                       vz, vz_lo, vz_hi, &
                       flxx, fx_lo, fx_hi, &
                       flxy, fy_lo, fy_hi, &
                       flxz, fz_lo, fz_hi, &
                       phix, phix_y, phix_z, &
                       phiy, phiy_x, phiy_z, &
                       phiz, phiz_x, phiz_y, &
                       slope, glo, ghi)

  ! Do a conservative update
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           uout(i,j,k) = uin(i,j,k) + &
                ( (flxx(i,j,k) - flxx(i+1,j,k)) * dtdx(1) &
                + (flxy(i,j,k) - flxy(i,j+1,k)) * dtdx(2) &
                + (flxz(i,j,k) - flxz(i,j,k+1)) * dtdx(3) )
        enddo
     enddo
  enddo

  ! Scale by face area in order to correctly reflx
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
           flxx(i,j,k) = flxx(i,j,k) * (dt * dx(2)*dx(3))
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)+1 
        do i = lo(1), hi(1)
           flxy(i,j,k) = flxy(i,j,k) * (dt * dx(1)*dx(3))
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)+1
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           flxz(i,j,k) = flxz(i,j,k) * (dt * dx(1)*dx(2))
        enddo
     enddo
  enddo

  call bl_deallocate(phix  )
  call bl_deallocate(phix_y)
  call bl_deallocate(phix_z)
  call bl_deallocate(phiy  )
  call bl_deallocate(phiy_x)
  call bl_deallocate(phiy_z)
  call bl_deallocate(phiz  )
  call bl_deallocate(phiz_x)
  call bl_deallocate(phiz_y)
  call bl_deallocate(slope)

end subroutine advect

subroutine advect_particles(particles, np, &
     ux, uxlo, uxhi, uy, uylo, uyhi, uz, uzlo, uzhi, dt, dx, plo) bind(c)
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use amrex_amrtracerparticlecontainer_module, only : amrex_tracerparticle
  
  integer,                    intent(in)            :: np
  type(amrex_tracerparticle), target, intent(inout) :: particles(np)
  integer,                    intent(in)            :: uxlo(3), uxhi(3)
  integer,                    intent(in)            :: uylo(3), uyhi(3)
  integer,                    intent(in)            :: uzlo(3), uzhi(3)
  real(amrex_real),           intent(in)            :: ux(uxlo(1):uxhi(1),uxlo(2):uxhi(2),uxlo(3):uxhi(3))
  real(amrex_real),           intent(in)            :: uy(uylo(1):uyhi(1),uylo(2):uyhi(2),uxlo(3):uxhi(3))
  real(amrex_real),           intent(in)            :: uz(uzlo(1):uzhi(1),uzlo(2):uzhi(2),uxlo(3):uxhi(3))
  real(amrex_real),           intent(in)            :: dt
  real(amrex_real),           intent(in)            :: plo(3)
  real(amrex_real),           intent(in)            :: dx(3)
  
  integer i, j, k, n, ipass
  real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
  real(amrex_real) lx, ly, lz
  real(amrex_real) inv_dx(AMREX_SPACEDIM)
  type(amrex_tracerparticle), pointer :: p
  
  inv_dx = 1.0d0/dx
  
  do ipass = 1, 2
     do n = 1, np
        p => particles(n)
        
        lx = (p%pos(1) - plo(1))*inv_dx(1) + 0.5d0
        ly = (p%pos(2) - plo(2))*inv_dx(2) + 0.5d0
        
        i = floor(lx)
        j = floor(ly)
        
        wx_hi = lx - i
        wy_hi = ly - j
        
        wx_lo = 1.0d0 - wx_hi
        wy_lo = 1.0d0 - wy_hi
        
        ! update the particle positions / velocities
        p%vel(1) = p%vel(1)
        p%vel(2) = p%vel(2)
        
        p%pos(1) = p%pos(1)
        p%pos(2) = p%pos(2)
     end do
  end do
end subroutine advect_particles

end module advect_module
