  module test_particle_module
    use amrex_fort_module, only: amrex_real, amrex_particle_real
    use iso_c_binding ,    only: c_int
    
    implicit none
    private
    
    public  particle_t
    
    type, bind(C)  :: particle_t
       real(amrex_particle_real) :: pos(3)     !< Position
       integer(c_int)            :: id         !< Particle id
       integer(c_int)            :: cpu        !< Particle cpu
    end type particle_t
    
  end module test_particle_module

subroutine move_particles(np, particles, vxp, vyp, vzp, dx)

  use amrex_fort_module, only : amrex_real
  use test_particle_module, only : particle_t
  implicit none

  integer,          intent(in), value :: np
  type(particle_t), intent(inout)     :: particles(np)
  real(amrex_real), intent(inout)     :: vxp(np), vyp(np), vzp(np)
  real(amrex_real), intent(in)        :: dx(3)

  integer                             :: ip
  real(amrex_real)                    :: velmag

!$acc parallel deviceptr(particles, vxp, vyp, vzp)
!$acc loop gang vector private(velmag)
  do ip = 1, np
     velmag = sqrt(vxp(ip)*vxp(ip) + vyp(ip)*vyp(ip) + vzp(ip)*vzp(ip))
     particles(ip)%pos(1) = particles(ip)%pos(1) + 0.5d0*vxp(ip)*dx(1) / velmag
     particles(ip)%pos(2) = particles(ip)%pos(2) + 0.5d0*vyp(ip)*dx(2) / velmag
     particles(ip)%pos(3) = particles(ip)%pos(3) + 0.5d0*vzp(ip)*dx(3) / velmag
  end do
!$acc end loop
!$acc end parallel

end subroutine move_particles
