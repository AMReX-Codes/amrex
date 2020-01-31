! Copyright 2019 Maxence Thevenet, Weiqun Zhang
!
! This file is part of WarpX.
!
! License: BSD-3-Clause-LBNL

module warpx_ES_push_particles

  use iso_c_binding
  use amrex_fort_module, only : amrex_real, amrex_particle_real

  implicit none

contains

!
! This routine updates the particle positions and velocities using the
! leapfrog time integration algorithm, given the electric fields at the
! particle positions. It also enforces specular reflection off the domain
! walls.
!
! Arguments:
!     particles : a pointer to the particle array-of-structs
!     ns        : the stride length of particle struct (the size of the struct in number of reals)
!     np        : the number of particles
!     vx_p      : the particle x-velocities
!     vy_p      : the particle y-velocities
!     vz_p      : the particle z-velocities
!     Ex_p      : the electric field in the x-direction at the particle positions
!     Ey_p      : the electric field in the y-direction at the particle positions
!     Ez_p      : the electric field in the z-direction at the particle positions
!     charge    : the charge of this particle species
!     mass      : the mass of this particle species
!     dt        : the time step
!     prob_lo   : the left-hand corner of the problem domain
!     prob_hi   : the right-hand corner of the problem domain
!
  subroutine warpx_push_leapfrog_3d(particles, ns, np,      &
                                    vx_p, vy_p, vz_p,       &
                                    Ex_p, Ey_p, Ez_p,       &
                                    charge, mass, dt,       &
                                    prob_lo, prob_hi)       &
       bind(c,name='warpx_push_leapfrog_3d')
    integer, value,   intent(in)     :: ns, np
    real(amrex_particle_real), intent(inout)  :: particles(ns,np)
    real(amrex_particle_real), intent(inout)  :: vx_p(np), vy_p(np), vz_p(np)
    real(amrex_real), intent(in)     :: Ex_p(np), Ey_p(np), Ez_p(np)
    real(amrex_real), intent(in)     :: charge
    real(amrex_real), intent(in)     :: mass
    real(amrex_real), intent(in)     :: dt
    real(amrex_real), intent(in)     :: prob_lo(3), prob_hi(3)

    integer n
    real(amrex_real) fac

    fac = charge * dt / mass

    do n = 1, np

       vx_p(n) = vx_p(n) + fac * Ex_p(n)
       vy_p(n) = vy_p(n) + fac * Ey_p(n)
       vz_p(n) = vz_p(n) + fac * Ez_p(n)

       particles(1, n) = particles(1, n) + dt * vx_p(n)
       particles(2, n) = particles(2, n) + dt * vy_p(n)
       particles(3, n) = particles(3, n) + dt * vz_p(n)

!      bounce off the walls in the x...
       do while (particles(1, n) .lt. prob_lo(1) .or. particles(1, n) .gt. prob_hi(1))
          if (particles(1, n) .lt. prob_lo(1)) then
             particles(1, n) = 2.d0*prob_lo(1) - particles(1, n)
          else
             particles(1, n) = 2.d0*prob_hi(1) - particles(1, n)
          end if
          vx_p(n) = -vx_p(n)
       end do

!      ... y...
       do while (particles(2, n) .lt. prob_lo(2) .or. particles(2, n) .gt. prob_hi(2))
          if (particles(2, n) .lt. prob_lo(2)) then
             particles(2, n) = 2.d0*prob_lo(2) - particles(2, n)
          else
             particles(2, n) = 2.d0*prob_hi(2) - particles(2, n)
          end if
          vy_p(n) = -vy_p(n)
       end do

!      ... and z directions
       do while (particles(3, n) .lt. prob_lo(3) .or. particles(3, n) .gt. prob_hi(3))
          if (particles(3, n) .lt. prob_lo(3)) then
             particles(3, n) = 2.d0*prob_lo(3) - particles(3, n)
          else
             particles(3, n) = 2.d0*prob_hi(3) - particles(3, n)
          end if
          vz_p(n) = -vz_p(n)
       end do

    end do

  end subroutine warpx_push_leapfrog_3d

  subroutine warpx_push_leapfrog_2d(particles, ns, np,      &
                                    vx_p, vy_p,             &
                                    Ex_p, Ey_p,             &
                                    charge, mass, dt,       &
                                    prob_lo, prob_hi)       &
       bind(c,name='warpx_push_leapfrog_2d')
    integer, value,   intent(in)     :: ns, np
    real(amrex_particle_real), intent(inout)  :: particles(ns,np)
    real(amrex_particle_real), intent(inout)  :: vx_p(np), vy_p(np)
    real(amrex_real), intent(in)     :: Ex_p(np), Ey_p(np)
    real(amrex_real), intent(in)     :: charge
    real(amrex_real), intent(in)     :: mass
    real(amrex_real), intent(in)     :: dt
    real(amrex_real), intent(in)     :: prob_lo(2), prob_hi(2)

    integer n
    real(amrex_real) fac

    fac = charge * dt / mass

    do n = 1, np

       vx_p(n) = vx_p(n) + fac * Ex_p(n)
       vy_p(n) = vy_p(n) + fac * Ey_p(n)

       particles(1, n) = particles(1, n) + dt * vx_p(n)
       particles(2, n) = particles(2, n) + dt * vy_p(n)

!      bounce off the walls in the x...
       do while (particles(1, n) .lt. prob_lo(1) .or. particles(1, n) .gt. prob_hi(1))
          if (particles(1, n) .lt. prob_lo(1)) then
             particles(1, n) = 2.d0*prob_lo(1) - particles(1, n)
          else
             particles(1, n) = 2.d0*prob_hi(1) - particles(1, n)
          end if
          vx_p(n) = -vx_p(n)
       end do

!      ... y...
       do while (particles(2, n) .lt. prob_lo(2) .or. particles(2, n) .gt. prob_hi(2))
          if (particles(2, n) .lt. prob_lo(2)) then
             particles(2, n) = 2.d0*prob_lo(2) - particles(2, n)
          else
             particles(2, n) = 2.d0*prob_hi(2) - particles(2, n)
          end if
          vy_p(n) = -vy_p(n)
       end do

    end do

  end subroutine warpx_push_leapfrog_2d

end module warpx_ES_push_particles
