
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
       &            flxx, fx_lo, fx_hi, &
       &            flxy, fy_lo, fy_hi, &
       &            dx,dt)

    use compute_flux_module, only : compute_flux_2d

    integer, intent(in) :: lo(2), hi(2)
    real(amrex_real), intent(in) :: dx(2), dt, time
    integer, intent(in) :: ui_lo(2), ui_hi(2)
    integer, intent(in) :: uo_lo(2), uo_hi(2)
    integer, intent(in) :: vx_lo(2), vx_hi(2)
    integer, intent(in) :: vy_lo(2), vy_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    real(amrex_real), intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2))
    real(amrex_real), intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
    real(amrex_real), intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))
    real(amrex_real), intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2))
    real(amrex_real), intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    real(amrex_real), intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
    
    integer :: i, j
    integer :: glo(2), ghi(2)
    real(amrex_real) :: dtdx(2), umax, vmax
    
    real(amrex_real), dimension(:,:), pointer, contiguous :: phix_1d, phiy_1d, phix, phiy, slope

    dtdx = dt/dx
    
    glo = lo - 1
    ghi = hi + 1

    ! edge states
    call amrex_allocate(phix_1d, glo, ghi)
    call amrex_allocate(phiy_1d, glo, ghi)
    call amrex_allocate(phix   , glo, ghi)
    call amrex_allocate(phiy   , glo, ghi)
    ! slope
    call amrex_allocate(slope  , glo, ghi)

    ! We like to allocate these **pointers** here and then pass them to a function
    ! to remove their pointerness for performance, because normally pointers could
    ! be aliasing.  We need to use pointers instead of allocatable arrays because
    ! we like to use AMReX's amrex_allocate to allocate memeory instead of the intrinsic
    ! allocate.  Amrex_allocate is much faster than allocate inside OMP.  
    ! Note that one MUST CALL AMREX_DEALLOCATE.
    
    ! check if CFL condition is violated.
    umax = maxval(abs(vx))
    vmax = maxval(abs(vy))
    if ( umax*dt .ge. dx(1) .or. &
         vmax*dt .ge. dx(2) ) then
       print *, "umax = ", umax, ", vmax = ", vmax, ", dt = ", dt, ", dx = ", dx
       call amrex_error("CFL violation. Use smaller adv.cfl.")
    end if

    ! call a function to compute flux
    call compute_flux_2d(lo, hi, dt, dx, &
                         uin, ui_lo, ui_hi, &
                         vx, vx_lo, vx_hi, &
                         vy, vy_lo, vy_hi, &
                         flxx, fx_lo, fx_hi, &
                         flxy, fy_lo, fy_hi, &
                         phix_1d, phiy_1d, phix, phiy, slope, glo, ghi)

    ! Do a conservative update
    do    j = lo(2),hi(2)
       do i = lo(1),hi(1)
          uout(i,j) = uin(i,j) + &
               ( (flxx(i,j) - flxx(i+1,j)) * dtdx(1) &
               + (flxy(i,j) - flxy(i,j+1)) * dtdx(2) )
       enddo
    enddo
    
    ! Scale by face area in order to correctly reflx
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          flxx(i,j) = flxx(i,j) * ( dt * dx(2))
       enddo
    enddo
    
    ! Scale by face area in order to correctly reflx
    do    j = lo(2), hi(2)+1 
       do i = lo(1), hi(1)
          flxy(i,j) = flxy(i,j) * (dt * dx(1))
       enddo
    enddo
    
    call amrex_deallocate(phix_1d)
    call amrex_deallocate(phiy_1d)
    call amrex_deallocate(phix)
    call amrex_deallocate(phiy)
    call amrex_deallocate(slope)
    
  end subroutine advect

  subroutine advect_particles(particles, np, &
       ux, uxlo, uxhi, uy, uylo, uyhi, dt, dx, plo) bind(c)
    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    use amrex_particlecontainer_module, only : amrex_particle

    integer,                    intent(in)            :: np
    type(amrex_particle),       intent(inout)         :: particles(np)
    integer,                    intent(in)            :: uxlo(2), uxhi(2)
    integer,                    intent(in)            :: uylo(2), uyhi(2)
    real(amrex_real), target,   intent(in)            :: ux(uxlo(1):uxhi(1),uxlo(2):uxhi(2))
    real(amrex_real), target,   intent(in)            :: uy(uylo(1):uyhi(1),uylo(2):uyhi(2))
    real(amrex_real),           intent(in)            :: dt
    real(amrex_real),           intent(in)            :: plo(2)
    real(amrex_real),           intent(in)            :: dx(2)
    
    integer cell(2)
    integer cc_cell(2)
    integer e_cell(2)
    
    integer ipass, n, d
    real(amrex_real) w_lo(2), w_hi(2)
    real(amrex_real) e_lo(2), e_hi(2)
    real(amrex_real) length(2)
    real(amrex_real) inv_dx(2)
    real(amrex_real) vel

    type dataptr
       real(amrex_real), dimension(:,:), pointer, contiguous :: p
    end type dataptr

    type(dataptr), dimension(2) :: velocity

    if (np == 0) then
       return
    end if
    
    velocity(1)%p => ux
    velocity(2)%p => uy
    
    inv_dx = 1.0d0/dx

    do ipass = 1, 2
       do n = 1, np

          length = (particles(n)%pos - plo)*inv_dx          

          cc_cell = floor(length)
          cell    = floor(length + 0.5d0)
          
          w_hi = length + 0.5d0 - cell          
          w_lo = 1.d0 - w_hi

          do d = 1, 2
             e_cell = cell
             e_cell(d) = cc_cell(d) + 1
             
             e_hi = w_hi
             e_lo = w_lo
             e_hi(d) = length(d) - cc_cell(d)

             e_hi = max(0.d0,min(1.d0,e_hi))             
             e_lo(d) = 1.d0 - e_hi(d)

             vel = e_lo(1)*e_lo(2)*velocity(d)%p(e_cell(1)-1, e_cell(2)-1) + &
                   e_lo(1)*e_hi(2)*velocity(d)%p(e_cell(1)-1, e_cell(2)  ) + &
                   e_hi(1)*e_lo(2)*velocity(d)%p(e_cell(1)  , e_cell(2)-1) + &
                   e_hi(1)*e_hi(2)*velocity(d)%p(e_cell(1)  , e_cell(2)  )

             if (ipass == 1) then
                particles(n)%vel(d) = particles(n)%pos(d)
                particles(n)%pos(d) = particles(n)%pos(d) + 0.5d0*dt*vel
             else
                particles(n)%pos(d) = particles(n)%vel(d) + dt*vel
                particles(n)%vel(d) = vel
             end if
          end do                    
       end do
    end do
  end subroutine advect_particles
  
end module advect_module
