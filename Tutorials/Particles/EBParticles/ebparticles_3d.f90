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

  end subroutine amrex_move_particles
