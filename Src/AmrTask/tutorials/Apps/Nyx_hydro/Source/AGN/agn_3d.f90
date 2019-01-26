  subroutine nyx_compute_overlap(np, particles, ng, ghosts, delta_x) &
       bind(c,name='nyx_compute_overlap')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    use particle_mod      , only: agn_particle_t

    integer             , intent(in   ) :: np, ng
    type(agn_particle_t), intent(inout) :: particles(np)
    type(agn_particle_t), intent(in   ) :: ghosts(ng)
    real(amrex_real)    , intent(in   ) :: delta_x(3)

    real(amrex_real) r2, cutoff
    integer i, j

    cutoff = delta_x(1)
    
    do i = 1, np-1
       do j = i+1, np

          r2 = sum((particles(i)%pos - particles(j)%pos)**2)

          if (r2 <= cutoff*cutoff) then
	     print *, "found overlap particles", particles(i)%id, particles(j)%id
	     ! If one of the particles is already invalidated, don't do anything
	     if (particles(i)%id .eq. -1 .or. particles(j)%id .eq. -1) cycle
             ! We only remove the newer (aka of lower mass) particle
             if (particles(i)%mass .lt. particles(j)%mass)  then
                particles(i)%id = -1
             else
                particles(j)%id = -1
             end if
          end if

       end do
    end do

    do i = 1, np
       do j = 1, ng

          r2 = sum((particles(i)%pos - ghosts(j)%pos)**2)

          if (r2 <= cutoff*cutoff) then
             print *, "found overlap ghost particles", particles(i)%id, ghosts(j)%id
             ! We only remove a particle if it is both 1) valid 2) newer (aka of lower mass)
             if (particles(i)%mass .lt. ghosts(j)%mass) then
                particles(i)%id = -1
             end if
             
          end if

      end do
    end do

  end subroutine nyx_compute_overlap

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine agn_merge_particles(np, particles, ng, ghosts, delta_x) &
       bind(c,name='agn_merge_particles')

    use iso_c_binding

    use amrex_fort_module, only : amrex_real
    use fundamental_constants_module, only: Gconst
    use particle_mod      , only: agn_particle_t
    use agn_params_module , only : l_merge, cutoff_vel

    integer             , intent(in   ) :: np, ng
    type(agn_particle_t), intent(inout) :: particles(np)
    type(agn_particle_t), intent(in   ) :: ghosts(ng)
    real(amrex_real)    , intent(in   ) :: delta_x(3)

    real(amrex_real) :: r2, vrelsq, r, mergetime
    real(amrex_real) :: cutoff, larger_mass
    integer :: i, j, merger_count
    !logical, save :: cutoff_vel=.false.

    cutoff = l_merge * delta_x(1)
    print *, "cutoff", cutoff
    mergetime = 10. *3.154e13 /3.086e19 * cutoff/1e-3   
    merger_count = 0

    do i = 1, np
       do j = i+1, np

          ! Distance between particles
          r2 = sum((particles(i)%pos - particles(j)%pos)**2)

          if (r2 <= cutoff*cutoff) then

	     if (cutoff_vel .eqv. .true.) then

                 ! Relative velocity
                 vrelsq = sum((particles(i)%vel - particles(j)%vel)**2)
                 !larger_mass = max(particles(i)%mass, particles(j)%mass)
                 r = sqrt(r2)
                 !if ( (vrelsq * r) < Gconst * larger_mass) then
	         if ( sqrt(vrelsq) * mergetime < r) then

		   ! If one of the particles is already invalidated, don't do anything
                    if (particles(i)%id .eq. -1 .or. particles(j)%id .eq. -1) cycle
                    ! Merge lighter particle into heavier one.
                    ! Set particle ID of lighter particle to -1
                    if (particles(i)%mass >= particles(j)%mass) then
                        call agn_merge_pair(particles(i), particles(j))
                        particles(j)%id = -1
		        merger_count = merger_count +1
                    else
                       call agn_merge_pair(particles(j), particles(i))
                       particles(i)%id = -1
		       merger_count = merger_count +1 
                    end if
                 end if
	     else ! just merge
		print *, "found merging particles", particles(i)%id, particles(j)%id
		if (particles(i)%id .eq. -1 .or. particles(j)%id .eq. -1) cycle
                if (particles(i)%mass >= particles(j)%mass) then
                    call agn_merge_pair(particles(i), particles(j))
                    particles(j)%id = -1
                    merger_count = merger_count +1
                else
                    call agn_merge_pair(particles(j), particles(i))
                    particles(i)%id = -1
                    merger_count = merger_count +1
                end if

            end if

	  end if

       end do
    end do

    print *, "number of mergers at this timestep", merger_count

    !this is merging ghost particles
    do i = 1, np
       do j = 1, ng

          r2 = sum((particles(i)%pos - ghosts(j)%pos)**2)

          if (r2 <= cutoff*cutoff) then
             if (cutoff_vel .eqv. .true.) then
             	vrelsq = sum((particles(i)%vel - ghosts(j)%vel)**2)
             	!larger_mass = max(particles(i)%mass, ghosts(j)%mass)

             	!if ( (vrelsq * sqrt(r2)) < Gconst * larger_mass) then
	        if ( sqrt(vrelsq) * mergetime <  sqrt(r2) ) then
		    !print *, "found merging ghost particles", particles(i)%id, ghosts(j)%id
                    if (particles(i)%mass > ghosts(j)%mass) then
                       ! The bigger particle "i" is in the valid region,
                       ! so we put all the mass onto it.
                       call agn_merge_pair(particles(i), ghosts(j))
                    else
                      ! The bigger particle "j" is in the ghost region --
                      ! we will let its grid update that particle,
                      ! so here we just invalidate particle "i".
                      particles(i)%id = -1
                   end if
               end if

	   else

	   	if (particles(i)%mass > ghosts(j)%mass) then
                    call agn_merge_pair(particles(i), ghosts(j))
                else
                    particles(i)%id = -1
                end if
	   end if

        end if

      end do
    end do

  end subroutine agn_merge_particles

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine agn_merge_pair(particle_stay, particle_remove)

    use amrex_fort_module, only : amrex_real
    use particle_mod      , only: agn_particle_t

    type(agn_particle_t), intent(inout) :: particle_stay
    type(agn_particle_t), intent(in   ) :: particle_remove

    real(amrex_real) :: mom(3)

    ! Calculate total momentum of both particles before merging.
    mom = particle_stay%vel * particle_stay%mass + &
         particle_remove%vel * particle_remove%mass

    ! Combine masses.
    particle_stay%mass = particle_stay%mass + particle_remove%mass

    ! Put total momentum onto the one particle, and save new velocity.
    particle_stay%vel = mom / particle_stay%mass
   
  end subroutine agn_merge_pair

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine agn_particle_velocity(np, particles, &
       state_old, sold_lo, sold_hi, &
       state_new, snew_lo, snew_hi, &
       dx, add_energy) &
       bind(c,name='agn_particle_velocity')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    use meth_params_module, only : NVAR, UMX, UMY, UMZ, UEDEN
    use particle_mod      , only: agn_particle_t

    integer,              intent(in   )        :: np
    integer,              intent(in   )        :: sold_lo(3), sold_hi(3)
    integer,              intent(in   )        :: snew_lo(3), snew_hi(3)
    integer,              intent(in   )        :: add_energy
    type(agn_particle_t), intent(inout)        :: particles(np)
    real(amrex_real),     intent(in   )        :: state_old &
         (sold_lo(1):sold_hi(1),sold_lo(2):sold_hi(2),sold_lo(3):sold_hi(3),NVAR)
    real(amrex_real),     intent(in   )        :: state_new &
         (snew_lo(1):snew_hi(1),snew_lo(2):snew_hi(2),snew_lo(3):snew_hi(3),NVAR)
    real(amrex_real),     intent(in   )        :: dx(3)

    integer          :: i, j, k, n
    real(amrex_real) :: vol, weight(-1:1, -1:1, -1:1)
    real(amrex_real) :: mass, momx, momy, momz, E, deltaEnergy, frac
    
    vol = dx(1) * dx(2) * dx(3)

    do n = 1, np

       call get_weights(weight, particles(n)%pos, dx)

       i = particles(n)%pos(1) / dx(1)
       j = particles(n)%pos(2) / dx(2)
       k = particles(n)%pos(3) / dx(3)

       ! momx, momy, momz: momentum = volume x change in momentum density.
       momx = sum((state_new(i-1:i+1, j-1:j+1, k-1:k+1, UMX) - &
                   state_old(i-1:i+1, j-1:j+1, k-1:k+1, UMX)) * weight) * vol
       momy = sum((state_new(i-1:i+1, j-1:j+1, k-1:k+1, UMY) - &
                   state_old(i-1:i+1, j-1:j+1, k-1:k+1, UMY)) * weight) * vol
       momz = sum((state_new(i-1:i+1, j-1:j+1, k-1:k+1, UMZ) - &
                   state_old(i-1:i+1, j-1:j+1, k-1:k+1, UMZ)) * weight) * vol

       !print *, "momentums", n, momx, momy, momz
       mass = particles(n)%mass

       ! Update velocity of particle so as to reduce momentum in the amount
       ! of the difference between old and new state.
       particles(n)%vel(1) = particles(n)%vel(1) - momx / mass
       particles(n)%vel(2) = particles(n)%vel(2) - momy / mass
       particles(n)%vel(3) = particles(n)%vel(3) - momz / mass

       ! Update particle energy if particle isn't brand new
       if (add_energy .gt. 0) then
          ! E: total energy = volume x change in total energy density.
          E = sum((state_new(i-1:i+1, j-1:j+1, k-1:k+1, UEDEN) - &
                   state_old(i-1:i+1, j-1:j+1, k-1:k+1, UEDEN)) * weight) * vol
          deltaEnergy = - E / mass
          frac = deltaEnergy / (particles(n)%energy)
          print *, 'deltaEnergy =', deltaEnergy, 'fraction', frac
          particles(n)%energy = particles(n)%energy + deltaEnergy
       end if

    end do

  end subroutine agn_particle_velocity

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine agn_accrete_mass(np, particles, &
       state, density_lost, slo, shi, &
       dt, dx) &
       bind(c,name='agn_accrete_mass')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    use fundamental_constants_module, only: Gconst, pi, eddington_const, c_light
    use eos_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UEDEN
    use particle_mod      , only : agn_particle_t
    use agn_params_module , only : eps_rad, eps_coupling, bondi_boost, max_frac_removed, frac_kinetic, eps_kinetic

    integer,              intent(in   )        :: np, slo(3), shi(3)
    type(agn_particle_t), intent(inout)        :: particles(np)
    ! FIXME: Should state be in only?
    real(amrex_real),     intent(inout)        :: state &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),NVAR)
    real(amrex_real),     intent(inout)        :: density_lost &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(amrex_real),     intent(in   )        :: dt, dx(3)

    integer          :: i, j, k, n, d, d_mom
    real(amrex_real) :: avg_rho, avg_csq, avg_speedsq
    real(amrex_real) :: c, denom, mass, mdot, m_edd

    !    real(amrex_real), parameter :: bondi_const = bondi_boost * 4.0d0*pi * Gconst*Gconst
    real(amrex_real) :: bondi_const
    real(amrex_real) :: speedsq &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(amrex_real) :: csq &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(amrex_real) :: vol, weight(-1:1, -1:1, -1:1)
    real(amrex_real) :: mass_change, speed_jet, E_kinetic, rvec(3)

    bondi_const = bondi_boost * 4.0d0*pi * Gconst*Gconst

    ! Remember state holds primitive variables:
    ! velocity, not momentum.
    speedsq = state(:,:,:,UMX)**2 + state(:,:,:,UMY)**2 + state(:,:,:,UMZ)**2

    do i = slo(1), shi(1)
       do j = slo(2), shi(2)
          do k = slo(3), shi(3)
             ! state(i,j,k,UEINT) is internal energy divided by density
             if (state(i,j,k,URHO) .lt. 0.) then
                print *, 'density = ', state(i,j,k,URHO)
                continue
             end if
             if (state(i,j,k,UEINT) .lt. 0.) then
                print *, 'energy = ', state(i,j,k,UEINT)
                continue
             end if
             call nyx_eos_soundspeed(c, state(i,j,k,URHO), state(i,j,k,UEINT))
             csq(i, j, k) = c*c
          end do
       end do
    end do

    vol = dx(1) * dx(2) * dx(3)

    do n = 1, np

       call get_weights(weight, particles(n)%pos, dx)

       i = particles(n)%pos(1) / dx(1)
       j = particles(n)%pos(2) / dx(2)
       k = particles(n)%pos(3) / dx(3)

       avg_rho =     sum(  state(i-1:i+1, j-1:j+1, k-1:k+1, URHO) * weight)
       avg_speedsq = sum(speedsq(i-1:i+1, j-1:j+1, k-1:k+1) * weight)
       avg_csq =     sum(    csq(i-1:i+1, j-1:j+1, k-1:k+1) * weight)
       
       !print *, 'rho, speed, cs', avg_rho, sqrt(avg_speedsq), sqrt(avg_csq)

       denom = (avg_speedsq + avg_csq)**1.5

       mass = particles(n)%mass
 
       ! Bondi accretion 
       mdot = bondi_const * mass * mass * avg_rho / denom

       ! Eddington limit
       m_edd = eddington_const * mass / eps_rad
       print *,'MDOT, MEDD: ', mdot, m_edd
       mdot = min(mdot, m_edd)

       mass_change = mdot * dt

       ! From each stencil cell, we'll reduce density by
       ! mass_change * weight[cell] / vol.
       ! But the most density that will be removed from cell is
       ! max_frac_removed * density[cell].
       ! This corresponds to maximum mass removal of
       ! Hence the most mass that will be removed from cell is
       ! max_frac_removed * density[cell] * vol.
       ! That is the maximum allowable mass_change * weight[cell].
       mass_change = MIN(mass_change, &
            max_frac_removed * &
            MAXVAL(state(i-1:i+1, j-1:j+1, k-1:k+1, URHO) * vol / weight))

       particles(n)%mdot = mass_change / dt

       ! Increase the mass of the particle by mass_change
       particles(n)%mass = particles(n)%mass + mass_change * ( 1.d0 - eps_rad)

       density_lost(i-1:i+1, j-1:j+1, k-1:k+1) = &
            density_lost(i-1:i+1, j-1:j+1, k-1:k+1) + &
            mass_change * weight / vol

       ! Thermal feedback energy.
       particles(n)%energy = particles(n)%energy + (1. - frac_kinetic) * &
            (mass_change * c_light**2) * eps_coupling * eps_rad

       ! Kinetic feedback: energy as well as momentum (hence velocity).
       if (frac_kinetic > 0.) then
          ! Kinetic energy added to the particle.
          E_kinetic = frac_kinetic * (mass_change * c_light**2)
          particles(n)%energy = particles(n)%energy + E_kinetic * eps_kinetic

          ! Change momentum of gas.
          speed_jet = sqrt(2. * E_kinetic / (avg_rho * vol))
          call get_random_direction(rvec)
          do d = 1, 3
             d_mom = UMX + d-1
             ! Update energy density by adding momentum density * velocity.
             ! Momentum density is M*L^-2*T^-1.
             ! Velocity is L*T^-1.
             ! Energy density is M*L^-1*T^-2, product of the two above.
             state(i-1:i+1, j-1:j+1, k-1:k+1, UEDEN) = &
               state(i-1:i+1, j-1:j+1, k-1:k+1, UEDEN) + &
               speed_jet * &
               sqrt(weight) * state(i-1:i+1, j-1:j+1, k-1:k+1, d_mom) * rvec(d)
          end do
          do d = 1, 3
             d_mom = UMX + d-1
             ! Update momentum density by adding density * velocity.
             ! Density is M*L^-3.
             ! Velocity is L*T^-1.
             ! Momentum density is M*L^-2*T^-1, product of the two above.
             state(i-1:i+1, j-1:j+1, k-1:k+1, d_mom) = &
               state(i-1:i+1, j-1:j+1, k-1:k+1, d_mom) + &
               speed_jet * &
               sqrt(weight) * state(i-1:i+1, j-1:j+1, k-1:k+1, URHO) * rvec(d)
          end do
       end if

    end do

  end subroutine agn_accrete_mass

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine agn_release_energy(np, particles, &
       state, slo, shi, &
       diag_eos, dlo, dhi, &
       a, dx) &
       bind(c,name='agn_release_energy')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    use fundamental_constants_module, only: k_B, m_proton
    use eos_module
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, NDIAG, NE_COMP
    use particle_mod      , only: agn_particle_t
    !use eos_module, only : nyx_eos_given_RT
    use agn_params_module, only : T_min

    integer,              intent(in   )        :: np, slo(3), shi(3)
    integer,              intent(in   )        :: dlo(3), dhi(3)
    type(agn_particle_t), intent(inout)        :: particles(np)
    real(amrex_real),     intent(inout)        :: state &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),NVAR)
    real(amrex_real),     intent(inout)        :: diag_eos &
         (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),NDIAG)
    real(amrex_real),     intent(in   )        :: a
    real(amrex_real),     intent(in   )        :: dx(3)

    integer          :: i, j, k, n

    real(amrex_real) :: vol, weight(-1:1, -1:1, -1:1)
    real(amrex_real) :: avg_rho, avg_Ne, pressure, e, m_g

    vol = dx(1) * dx(2) * dx(3)

    do n = 1, np

       call get_weights(weight, particles(n)%pos, dx)

       i = particles(n)%pos(1) / dx(1)
       j = particles(n)%pos(2) / dx(2)
       k = particles(n)%pos(3) / dx(3)

       avg_rho = sum(state(i-1:i+1, j-1:j+1, k-1:k+1, URHO) * weight)
       m_g = avg_rho * vol

       avg_Ne = sum(diag_eos(i-1:i+1, j-1:j+1, k-1:k+1, NE_COMP) * weight)

       call nyx_eos_given_RT(e, pressure, avg_rho, T_min, avg_Ne, a)

       print *, 'AGN particle at ', particles(n)%pos, ':', i, j, k
       print 50, particles(n)%pos, i, j, k, particles(n)%mass, &
            particles(n)%energy
50     format (1x, 'AGN particle at ', 3F8.3, 3I4, ' m=', E12.5, ' e=', E12.5)

       if (particles(n)%energy > m_g * e) then

          print *, 'RELEASING ENERGY of particle at ', particles(n)%pos
          print *, 'neighborhood mass: ', m_g
          print *, 'e = ', e
          print *, 'particle energy: ', particles(n)%energy
          print *, 'm_g * e = ', (m_g * e)

          state(i-1:i+1, j-1:j+1, k-1:k+1, UEDEN) = &
          state(i-1:i+1, j-1:j+1, k-1:k+1, UEDEN) + &
          particles(n)%energy * weight

          state(i-1:i+1, j-1:j+1, k-1:k+1, UEINT) = &
          state(i-1:i+1, j-1:j+1, k-1:k+1, UEINT) + &
          particles(n)%energy * weight

          particles(n)%energy = 0.
               
       end if

    end do

  end subroutine agn_release_energy

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine get_length_frac(frac, x, dx)

    use amrex_fort_module, only : amrex_real, amrex_particle_real
    use bl_constants_module, only : ZERO, ONE, HALF

    real(amrex_real), intent(out  )  :: frac(-1:1)
    real(amrex_particle_real), intent(in   )  :: x
    real(amrex_real), intent(in   )  :: dx

    integer :: i
    real(amrex_real) :: offset

    i = x / dx
    offset = x / dx - i
    if (offset < HALF) then ! offset in range 0 : 0.5
       frac(-1) = HALF - offset
       frac(0) = HALF + offset
       frac(1) = ZERO
    else ! offset in range 0.5 : 1
       frac(-1) = ZERO
       frac(0) = ONE + HALF - offset
       frac(1) = offset - HALF
    end if
    
  end subroutine get_length_frac

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine get_weights(weight, pos, dx)

    use amrex_fort_module, only : amrex_real, amrex_particle_real 
    use bl_constants_module, only : ZERO, ONE, HALF

    real(amrex_real), intent(out  )  :: weight(-1:1, -1:1, -1:1)
    real(amrex_particle_real), intent(in   )  :: pos(3)
    real(amrex_real), intent(in   )  :: dx(3)

    integer :: d, i, j, k
    real(amrex_real) :: frac(-1:1, 3)

    do d = 1, 3
       call get_length_frac(frac(:, d), pos(d), dx(d))
    end do

    do k = -1, 1
       do j = -1, 1
          do i = -1, 1
             weight(i, j, k) = frac(i, 1) * frac(j, 2) * frac(k, 3)
          end do
       end do
    end do

  end subroutine get_weights

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine get_random_direction(vec)

    use amrex_fort_module, only : amrex_real
    use fundamental_constants_module, only : pi
    use uniform01_module, only : get_uniform01

    real(amrex_real), intent(out)  :: vec(3)
    real(amrex_real) :: u, v, theta, horiz

    u = get_uniform01()
    v = get_uniform01()

    theta = 2. * pi * u
    horiz = 2. * sqrt(v * (1. - v))

    vec(1) = cos(theta) * horiz
    vec(2) = sin(theta) * horiz
    vec(3) = 2. * v - 1.

  end subroutine get_random_direction
