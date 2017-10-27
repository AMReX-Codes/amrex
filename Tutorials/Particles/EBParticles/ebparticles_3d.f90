  module eb_particle_module
    use amrex_fort_module, only: amrex_real, amrex_particle_real
    use iso_c_binding ,    only: c_int
    
    implicit none
    private
    
    public  particle_t
    
    type, bind(C)  :: particle_t
       real(amrex_particle_real) :: pos(3)     !< Position
       real(amrex_particle_real) :: vel(3)     !< Particle velocity
       real(amrex_particle_real) :: acc(3)     !< Particle acceleration
       integer(c_int)            :: id         !< Particle id
       integer(c_int)            :: cpu        !< Particle cpu
    end type particle_t
    
  end module eb_particle_module

  subroutine amrex_move_particles(particles, np, dt, prob_lo, prob_hi) &
       bind(c,name='amrex_move_particles')
    
    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    use eb_particle_module, only : particle_t

    implicit none

    integer,          intent(in   )         :: np
    type(particle_t), intent(inout), target :: particles(np)
    real(amrex_real), intent(in   )         :: dt
    real(amrex_real), intent(in   )         :: prob_lo(3), prob_hi(3)

    integer i
    type(particle_t), pointer :: p

    do i = 1, np
       
       p => particles(i)

!      update the particle positions / velocities
       p%vel(1) = p%vel(1) + p%acc(1) * dt 
       p%vel(2) = p%vel(2) + p%acc(2) * dt 
       p%vel(3) = p%vel(3) + p%acc(3) * dt 

       p%pos(1) = p%pos(1) + p%vel(1) * dt 
       p%pos(2) = p%pos(2) + p%vel(2) * dt 
       p%pos(3) = p%pos(3) + p%vel(3) * dt 

    end do

  end subroutine amrex_move_particles

  subroutine amrex_bounce_walls(particles, np, plo, dx, dt, &
       flag, fglo, fghi, bcent, blo, bhi, apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi) &
       bind(c,name='amrex_bounce_walls')
    
    use iso_c_binding
    use amrex_fort_module,       only : amrex_real
    use eb_particle_module,      only : particle_t
    use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell
    use amrex_error_module,      only : amrex_abort

    integer,          intent(in   )         :: np
    type(particle_t), intent(inout), target :: particles(np)
    integer, dimension(3), intent(in) :: axlo,axhi,aylo,ayhi,azlo,azhi,fglo,fghi,blo,bhi
    integer,          intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
    real(amrex_real), intent(in) :: bcent  (blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3)
    real(amrex_real), intent(in) :: apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(amrex_real), intent(in) :: apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(amrex_real), intent(in) :: apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(amrex_real), intent(in) :: plo(3)
    real(amrex_real), intent(in) :: dx(3), dt

    real(amrex_real) inv_dx(3)
    real(amrex_real) :: lx, ly, lz
    real(amrex_real) :: axm, axp, aym, ayp, azm, azp
    real(amrex_real) :: speed, vxnorm, vynorm, vznorm, dotp
    real(amrex_real) :: apnorm, apnorminv, anrmx, anrmy, anrmz
    real(amrex_real) :: bcentx, bcenty, bcentz, d
    integer i, j, k, n
    type(particle_t), pointer :: p

    inv_dx = 1.0d0 / dx

    do n = 1, np

       p => particles(n)

       lx = (p%pos(1) - plo(1))*inv_dx(1)
       ly = (p%pos(2) - plo(2))*inv_dx(2)
       lz = (p%pos(3) - plo(3))*inv_dx(3)
    
       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       if (is_regular_cell(flag(i,j,k))) then
          cycle
       end if

       if (is_covered_cell(flag(i,j,k))) then

          p%pos(1) = p%pos(1) - dt * p%vel(1)
          p%pos(2) = p%pos(2) - dt * p%vel(2)
          p%pos(3) = p%pos(3) - dt * p%vel(3)

          lx = (p%pos(1) - plo(1))*inv_dx(1)
          ly = (p%pos(2) - plo(2))*inv_dx(2)
          lz = (p%pos(3) - plo(3))*inv_dx(3)
          
          i = floor(lx)
          j = floor(ly)
          k = floor(lz)

          p%pos(1) = p%pos(1) + dt * p%vel(1)
          p%pos(2) = p%pos(2) + dt * p%vel(2)
          p%pos(3) = p%pos(3) + dt * p%vel(3)
       end if

! this is a cut cell. compute normal

       axm = apx(i,  j  ,k  )
       axp = apx(i+1,j  ,k  )
       aym = apy(i,  j  ,k  )
       ayp = apy(i,  j+1,k  )
       azm = apz(i,  j  ,k  )
       azp = apz(i,  j  ,k+1)

       apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2 + (azm-azp)**2)

       if (apnorm .eq. 0.d0) then
          print *, "bounce_walls: ", axm, axp, aym, ayp, azm, azp
          flush(6)
          call amrex_abort("bounce_walls: we are in trouble.")
       end if

       apnorminv = 1.d0 / apnorm
       anrmx = (axp-axm) * apnorminv   ! pointing to the wall
       anrmy = (ayp-aym) * apnorminv
       anrmz = (azp-azm) * apnorminv

! convert bcent to global coordinate system centered at plo
       bcentx = bcent(i, j, k, 1)*dx(1) + plo(1) + i*dx(1) + 0.5d0*dx(1)
       bcenty = bcent(i, j, k, 2)*dx(2) + plo(2) + j*dx(2) + 0.5d0*dx(2)
       bcentz = bcent(i, j, k, 3)*dx(3) + plo(3) + k*dx(3) + 0.5d0*dx(3)

! distance to boundary
       d = (p%pos(1) - bcentx) * anrmx + (p%pos(2) - bcenty) * anrmy  + &
           (p%pos(3) - bcentz) * anrmz

! if we are outside, bounce in
       if (d .lt. 0.d0) then
          p%pos(1) = p%pos(1) - 2.d0*d*anrmx
          p%pos(2) = p%pos(2) - 2.d0*d*anrmy
          p%pos(3) = p%pos(3) - 2.d0*d*anrmz

          speed = sqrt(p%vel(1)**2 + p%vel(2)**2 + p%vel(3)**2)
          vxnorm = p%vel(1) / speed
          vynorm = p%vel(2) / speed
          vznorm = p%vel(3) / speed

          dotp = anrmx*vxnorm + anrmy*vynorm + anrmz*vznorm

          vxnorm = 2.d0*dotp*anrmx - vxnorm
          vynorm = 2.d0*dotp*anrmy - vynorm
          vznorm = 2.d0*dotp*anrmz - vznorm

          p%vel(1) = -vxnorm*speed
          p%vel(2) = -vynorm*speed
          p%vel(3) = -vznorm*speed

       end if

    end do

  end subroutine amrex_bounce_walls
