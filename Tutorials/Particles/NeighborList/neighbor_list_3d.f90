  module short_range_particle_module
    use amrex_fort_module, only: amrex_real, amrex_particle_real
    use iso_c_binding ,    only: c_int, c_float, c_double
    
    implicit none
    private
    
    public  particle_t, neighbor_t
    
    type, bind(C)  :: particle_t
       real(amrex_particle_real) :: pos(3)     !< Position
       real(amrex_particle_real) :: vel(3)     !< Particle velocity
       real(amrex_particle_real) :: acc(3)     !< Particle acceleration
       integer(c_int)            :: id         !< Particle id
       integer(c_int)            :: cpu        !< Particle cpu
    end type particle_t
    
    type, bind(C)  :: neighbor_t
       real(amrex_particle_real) :: pos(3)     !< Position
       real(amrex_particle_real) :: vel(3)     !< Particle velocity
       real(amrex_particle_real) :: acc(3)     !< Particle acceleration
       integer(c_int)            :: id         !< Particle id
       integer(c_int)            :: cpu        !< Particle cpu
    end type neighbor_t
    
  end module short_range_particle_module
  
  subroutine amrex_move_particles(particles, np, dt, prob_lo, prob_hi) &
       bind(c,name='amrex_move_particles')
    
    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    use short_range_particle_module, only : particle_t

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

!      bounce off the walls in the x...
       do while (p%pos(1) .lt. prob_lo(1) .or. p%pos(1) .gt. prob_hi(1))
          if (p%pos(1) .lt. prob_lo(1)) then
             p%pos(1) = 2.d0*prob_lo(1) - p%pos(1)
          else
             p%pos(1) = 2.d0*prob_hi(1) - p%pos(1)
          end if
          p%vel(1) = -p%vel(1)
       end do

!      ... y .. 
       do while (p%pos(2) .lt. prob_lo(2) .or. p%pos(2) .gt. prob_hi(2))
          if (p%pos(2) .lt. prob_lo(2)) then
             p%pos(2) = 2.d0*prob_lo(2) - p%pos(2)
          else
             p%pos(2) = 2.d0*prob_hi(2) - p%pos(2)
          end if
          p%vel(2) = -p%vel(2)
       end do

!      ... and z directions
       do while (p%pos(3) .lt. prob_lo(3) .or. p%pos(3) .gt. prob_hi(3))
          if (p%pos(3) .lt. prob_lo(3)) then
             p%pos(3) = 2.d0*prob_lo(3) - p%pos(3)
          else
             p%pos(3) = 2.d0*prob_hi(3) - p%pos(3)
          end if
          p%vel(3) = -p%vel(3)
       end do

    end do

  end subroutine amrex_move_particles

  subroutine amrex_compute_forces(particles, np, neighbors, nn, & 
                                 cutoff, min_r) &
       bind(c,name='amrex_compute_forces')

    use iso_c_binding
    use amrex_fort_module,           only : amrex_real
    use short_range_particle_module, only : particle_t, neighbor_t
        
    integer,          intent(in   ) :: np, nn
    real(amrex_real), intent(in   ) :: cutoff, min_r
    type(particle_t), intent(inout) :: particles(np)
    type(neighbor_t), intent(inout) :: neighbors(nn)

    real(amrex_real) dx, dy, dz, r2, r, coef, mass
    integer i, j

    mass   = 1.d-2
    
    do i = 1, np

!      zero out the particle acceleration
       particles(i)%acc(1) = 0.d0
       particles(i)%acc(2) = 0.d0
       particles(i)%acc(3) = 0.d0

       do j = 1, np

          if (i .eq. j ) then
             cycle
          end if

          dx = particles(i)%pos(1) - particles(j)%pos(1)
          dy = particles(i)%pos(2) - particles(j)%pos(2)
          dz = particles(i)%pos(3) - particles(j)%pos(3)

          r2 = dx * dx + dy * dy + dz * dz

          if (r2 .gt. cutoff*cutoff) then
             cycle
          end if

          r2 = max(r2, min_r*min_r) 
          r = sqrt(r2)

          coef = (1.d0 - cutoff / r) / r2 / mass
          particles(i)%acc(1) = particles(i)%acc(1) + coef * dx
          particles(i)%acc(2) = particles(i)%acc(2) + coef * dy
          particles(i)%acc(3) = particles(i)%acc(3) + coef * dy

       end do

       do j = 1, nn

          dx = particles(i)%pos(1) - neighbors(j)%pos(1)
          dy = particles(i)%pos(2) - neighbors(j)%pos(2)
          dy = particles(i)%pos(3) - neighbors(j)%pos(3)

          r2 = dx * dx + dy * dy + dz * dz

          if (r2 .gt. cutoff*cutoff) then
             cycle
          end if

          r2 = max(r2, min_r*min_r) 
          r = sqrt(r2)

          coef = (1.d0 - cutoff / r) / r2 / mass
          particles(i)%acc(1) = particles(i)%acc(1) + coef * dx
          particles(i)%acc(2) = particles(i)%acc(2) + coef * dy
          particles(i)%acc(3) = particles(i)%acc(3) + coef * dz
          
      end do
    end do

  end subroutine amrex_compute_forces

  subroutine amrex_compute_forces_nl(rparticles, np, neighbors, & 
                                     nn, nl, size, cutoff, min_r) &
       bind(c,name='amrex_compute_forces_nl')

    use iso_c_binding
    use amrex_fort_module,           only : amrex_real
    use short_range_particle_module, only : particle_t
        
    integer,          intent(in   ) :: np, nn, size
    real(amrex_real), intent(in   ) :: cutoff, min_r
    type(particle_t), intent(inout) :: rparticles(np)
    type(particle_t), intent(inout) :: neighbors(nn)
    integer,          intent(in   ) :: nl(size)

    real(amrex_real) dx, dy, dz, r2, r, coef, mass
    integer i, j, index, nneighbors

    type(particle_t)                    :: particles(np+nn)
        
    particles(    1:np) = rparticles
    particles(np+1:   ) = neighbors

    mass   = 1.d-2
    
    index = 1
    do i = 1, np

!      zero out the particle acceleration
       particles(i)%acc(1) = 0.d0
       particles(i)%acc(2) = 0.d0
       particles(i)%acc(3) = 0.d0

       nneighbors = nl(index)
       index = index + 1

       do j = index, index + nneighbors - 1

          dx = particles(i)%pos(1) - particles(nl(j))%pos(1)
          dy = particles(i)%pos(2) - particles(nl(j))%pos(2)
          dz = particles(i)%pos(3) - particles(nl(j))%pos(3)

          r2 = dx * dx + dy * dy + dz * dz
          r2 = max(r2, min_r*min_r) 
          r = sqrt(r2)

          coef = (1.d0 - cutoff / r) / r2 / mass
          particles(i)%acc(1) = particles(i)%acc(1) + coef * dx
          particles(i)%acc(2) = particles(i)%acc(2) + coef * dy
          particles(i)%acc(3) = particles(i)%acc(3) + coef * dz

       end do

       index = index + nneighbors

    end do

    rparticles(:) = particles(1:np)
    neighbors(:)  = particles(np+1:)

  end subroutine amrex_compute_forces_nl
