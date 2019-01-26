module comoving_nd_module

  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real

  contains

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_integrate_comoving_a(old_a,new_a,dt) &
         bind(C, name="fort_integrate_comoving_a")

        use fundamental_constants_module, only: Hubble_const
        use comoving_module             , only: comoving_h, comoving_OmM, comoving_type

        implicit none

        real(rt), intent(in   ) :: old_a, dt
        real(rt), intent(  out) :: new_a

        real(rt), parameter :: xacc = 1.0d-8
        real(rt) :: H_0, OmL
        real(rt) :: Delta_t, prev_soln
        real(rt) :: start_a, end_a, start_slope, end_slope
        integer  :: iter, j, nsteps

        if (comoving_h .eq. 0.0d0) then
          new_a = old_a
          return
        endif

        H_0 = comoving_h * Hubble_const
        OmL = 1.d0 - comoving_OmM 

        prev_soln = 2.0d0 ! 0<a<1 so a=2 will do as "wrong" solution
        do iter = 1, 11  ! max allowed iterations
          nsteps = 2**(iter-1)
          Delta_t = dt/nsteps
          end_a = old_a

          do j = 1, nsteps
            ! This uses RK2 to integrate the ODE:
            !   da / dt = H_0 * sqrt(OmM/a + OmL*a^2)
            start_a = end_a

            ! Compute the slope at the old time
            if (comoving_type > 0) then
                start_slope = H_0*dsqrt(comoving_OmM / start_a + OmL*start_a**2)
            else
                start_slope = comoving_h
            end if

            ! Compute a provisional value of ln(a) at the new time 
            end_a = start_a + start_slope * Delta_t

            ! Compute the slope at the new time
            if (comoving_type > 0) then
                end_slope = H_0*dsqrt(comoving_OmM / end_a + OmL*end_a**2)
            else
                end_slope = comoving_h 
            end if
       
            ! Now recompute a at the new time using the average of the two slopes
            end_a = start_a + 0.5d0 * (start_slope + end_slope) * Delta_t
          enddo

          new_a  = end_a
          if (abs(1.0d0-new_a/prev_soln) .le. xacc) return
          prev_soln = new_a

        enddo

      end subroutine fort_integrate_comoving_a

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_integrate_comoving_a_to_z(old_a,z_value,dt) &
         bind(C, name="fort_integrate_comoving_a_to_z")

        use fundamental_constants_module, only: Hubble_const
        use comoving_module             , only: comoving_h, comoving_OmM, comoving_type

        implicit none

        real(rt), intent(in   ) :: old_a, z_value
        real(rt), intent(inout) :: dt

        real(rt), parameter :: xacc = 1.0d-8
        real(rt) :: H_0, OmL
        real(rt) :: Delta_t
        real(rt) :: start_a, end_a, start_slope, end_slope
        real(rt) :: a_value
        integer  :: j, nsteps

        if (comoving_h .eq. 0.0d0) &
          call amrex_error("fort_integrate_comoving_a_to_z: Shouldn't be setting plot_z_values if not evolving a")

        H_0 = comoving_h * Hubble_const
        OmL = 1.d0 - comoving_OmM 
      
        ! Translate the target "z" into a target "a"
        a_value = 1.d0 / (1.d0 + z_value)

        ! Use lots of steps if we want to nail the z_value
        nsteps = 1024

        ! We integrate a, but stop when a = a_value (or close enough)
        Delta_t = dt/nsteps
        end_a = old_a
        do j = 1, nsteps
            ! This uses RK2 to integrate the ODE:
            !   da / dt = H_0 * sqrt(OmM/a + OmL*a^2)
            start_a = end_a

            ! Compute the slope at the old time
            if (comoving_type > 0) then
                start_slope = H_0*dsqrt(comoving_OmM / start_a + OmL*start_a**2)
            else
                start_slope = comoving_h
            end if

            ! Compute a provisional value of ln(a) at the new time 
            end_a = start_a + start_slope * Delta_t

            ! Compute the slope at the new time
            if (comoving_type > 0) then
                end_slope = H_0*dsqrt(comoving_OmM / end_a + OmL*end_a**2)
            else
                end_slope = comoving_h 
            end if
       
            ! Now recompute a at the new time using the average of the two slopes
            end_a = start_a + 0.5d0 * (start_slope + end_slope) * Delta_t

            ! We have crossed from a too small to a too big in this step
            if ( (end_a - a_value) * (start_a - a_value) < 0) then
                dt = ( (  end_a - a_value) * dble(j  ) + &
                       (a_value - start_a) * dble(j+1) ) / (end_a - start_a) * Delta_t
                exit
            end if
        end do

      end subroutine fort_integrate_comoving_a_to_z
! :::
! ::: ----------------------------------------------------------------
! :::

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_integrate_comoving_a_to_a(old_a,a_value,dt) &
         bind(C, name="fort_integrate_comoving_a_to_a")

        use fundamental_constants_module, only: Hubble_const
        use comoving_module             , only: comoving_h, comoving_OmM, comoving_type

        implicit none

        real(rt), intent(in   ) :: old_a, a_value
        real(rt), intent(inout) :: dt

        real(rt), parameter :: xacc = 1.0d-8
        real(rt) :: H_0, OmL
        real(rt) :: Delta_t
        real(rt) :: start_a, end_a, start_slope, end_slope
        integer  :: j, nsteps

        real(rt) :: max_dt

        if (comoving_h .eq. 0.0d0) &
          call amrex_error("fort_integrate_comoving_a_to_z: Shouldn't be setting plot_z_values if not evolving a")

        H_0 = comoving_h * Hubble_const
        OmL = 1.d0 - comoving_OmM 
      

        ! Use lots of steps if we want to nail the z_value
!        nsteps = 1024

        ! Use enough steps if we want to be close to the a_value
        nsteps = 2048

        ! We integrate a, but stop when a = a_value (or close enough)
        Delta_t = dt/nsteps
        end_a = old_a
        do j = 1, nsteps
            ! This uses RK2 to integrate the ODE:
            !   da / dt = H_0 * sqrt(OmM/a + OmL*a^2)
            start_a = end_a

            ! Compute the slope at the old time
            if (comoving_type > 0) then
                start_slope = H_0*dsqrt(comoving_OmM / start_a + OmL*start_a**2)
            else
                start_slope = comoving_h
            end if

            ! Compute a provisional value of ln(a) at the new time 
            end_a = start_a + start_slope * Delta_t

            ! Compute the slope at the new time
            if (comoving_type > 0) then
                end_slope = H_0*dsqrt(comoving_OmM / end_a + OmL*end_a**2)
            else
                end_slope = comoving_h 
            end if
       
            ! Now recompute a at the new time using the average of the two slopes
            end_a = start_a + 0.5d0 * (start_slope + end_slope) * Delta_t

            ! We have crossed from a too small to a too big in this step
            if ( (end_a - a_value) * (start_a - a_value) < 0) then
                dt = ( (  end_a - a_value) * dble(j  ) + &
                       (a_value - start_a) * dble(j+1) ) / (end_a - start_a) * Delta_t
                exit
            end if
        end do

      end subroutine fort_integrate_comoving_a_to_a
! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_est_maxdt_comoving_a(old_a,dt)

        use fundamental_constants_module, only: Hubble_const
        use comoving_module             , only: comoving_h, comoving_OmM, comoving_type

        implicit none

        real(rt), intent(in   ) :: old_a
        real(rt), intent(inout) :: dt

        real(rt) :: H_0, OmL
        real(rt) :: max_dt

        OmL = 1.d0 - comoving_OmM 

        ! This subroutine computes dt based on not changing a by more than 5% 
        ! if we use forward Euler integration
        !   d(ln(a)) / dt = H_0 * sqrt(OmM/a^3 + OmL)

        H_0 = comoving_h * Hubble_const

        if (H_0 .ne. 0.0d0) then
           if (comoving_type > 0) then
              max_dt = (0.05d0) / H_0 / dsqrt(comoving_OmM / old_a**3 + OmL)
           else
              max_dt = (0.05d0) / abs(comoving_h)
           end if
           dt = min(dt,max_dt) 

        else 

           ! dt is unchanged

        end if

      end subroutine fort_est_maxdt_comoving_a

! :::
! ::: ----------------------------------------------------------------
! :::

      ! This might only work for t=0=> a=.00625, although constant canceled
      subroutine fort_est_lindt_comoving_a(old_a,new_a,dt)

        use fundamental_constants_module, only: Hubble_const
        use comoving_module             , only: comoving_h, comoving_OmM, comoving_type

        implicit none

        real(rt), intent(in   ) :: old_a,new_a
        real(rt), intent(inout) :: dt

        real(rt) :: H_0, OmL
        real(rt) :: lin_dt

        OmL = 1.d0 - comoving_OmM 

        ! This subroutine computes dt based on not changing a by more than 5% 
        ! if we use forward Euler integration
        !   d(ln(a)) / dt = H_0 * sqrt(OmM/a^3 + OmL)

        H_0 = comoving_h * Hubble_const

        ! Could definately be optimized better
        if (H_0 .ne. 0.0d0) then

           lin_dt= ((new_a/(.75**(2/3)*(OmL+ comoving_OmM)**(1/3)))**(.75)  - &
                ((old_a/(.75**(2/3)*(OmL+ comoving_OmM)**(1/3)))**(.75) ) ) /H_0
           dt=lin_dt
           
        else 

           ! dt is unchanged

        end if

      end subroutine fort_est_lindt_comoving_a

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_estdt_comoving_a(old_a,new_a,dt,change_allowed,fixed_da,final_a,dt_modified) &
         bind(C, name="fort_estdt_comoving_a")

        use comoving_module             , only: comoving_h

        implicit none

        real(rt), intent(in   ) :: old_a, change_allowed, fixed_da, final_a
        real(rt), intent(inout) :: dt
        real(rt), intent(  out) :: new_a
        integer , intent(  out) :: dt_modified
        real(rt) a_value
        real(rt) max_dt
        max_dt = dt

        if (comoving_h .ne. 0.0d0) then

           if( fixed_da .le. 0.0d0) then
              ! First call this to make sure dt that we send to integration routine isnt outrageous
              call fort_est_maxdt_comoving_a(old_a,dt)
              
              ! Initial call to see if existing dt will work
              call fort_integrate_comoving_a(old_a,new_a,dt)
              
              ! Make sure a isn't growing too fast
              call enforce_percent_change(old_a,new_a,dt,change_allowed)
           else
              ! First call this to make sure dt that we send to integration routine isnt outrageous
              new_a = (old_a +  fixed_da);
              call fort_est_lindt_comoving_a(old_a,new_a,dt)             
              call fort_est_maxdt_comoving_a(old_a,dt)

              ! Then integrate old_a to a_value using dt as a guess for the maximum dt
              ! Output dt is based on a fraction of the input dt
              call fort_integrate_comoving_a_to_a(old_a,new_a,dt)
           endif           

           ! Make sure we don't go past final_a (if final_a is set)
           if (final_a > 0.0d0) &
              call enforce_final_a(old_a,new_a,dt,final_a)

           dt_modified = 1

        else

           ! dt is unchanged by this call

           dt_modified = 0

        endif 

      end subroutine fort_estdt_comoving_a

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine enforce_percent_change(old_a,new_a,dt,change_allowed)

        implicit none

        real(rt), intent(in   ) :: old_a, change_allowed
        real(rt), intent(inout) :: dt
        real(rt), intent(inout) :: new_a

        integer          :: i
        real(rt) :: factor

        factor = ( (new_a - old_a) / old_a ) / change_allowed

        ! Only go into this process if percent change exceeds change_allowed

        if (factor > 1.d0) then

           do i = 1, 100
              factor = ( (new_a - old_a) / old_a ) / change_allowed

              ! Note: 0.99 is just a fudge factor so we don't get bogged down.
              if (factor > 1.d0) then
                 dt = (1.d0 / factor) * dt * 0.99d0
                 call fort_integrate_comoving_a(old_a,new_a,dt)
              else if (i.lt.100) then
                 call fort_integrate_comoving_a(old_a,new_a,dt)
                 ! We're done
                 return 
              else
                 call amrex_error("Too many iterations in enforce_percent_change")
              end if
           end do

        else
           ! We don't need to do anything
           return 
        end if

      end subroutine enforce_percent_change

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine enforce_final_a(old_a,new_a,dt,final_a)

        implicit none

        real(rt), intent(in   ) :: old_a, final_a
        real(rt), intent(inout) :: dt
        real(rt), intent(inout) :: new_a

        integer  :: i
        real(rt) :: factor
        real(rt), parameter :: eps = 1.d-10

        if (old_a > final_a) then
           call amrex_error("Oops -- old_a > final_a")
        end if

        ! Only go into this process if new_a is past final_a
        if (new_a > final_a) then

           do i = 1, 100
              if ( (new_a > (final_a+eps)) .or. (new_a < final_a) ) then
                 factor = (final_a - old_a) / (new_a - old_a)
                 dt = dt * factor 
                 call fort_integrate_comoving_a(old_a,new_a,dt)
              else if (i.lt.100) then
                 ! We're done
                 return 
              else
                 call amrex_error("Too many iterations in enforce_final_a")
              end if
           end do

        else
           ! We don't need to do anything
           return 
        end if

      end subroutine enforce_final_a

! ! :::
! ! ::: ----------------------------------------------------------------
! ! :::

!       subroutine fort_get_omb(frac) &
!          bind(C, name="fort_get_omb")

!         use comoving_module, only: comoving_OmB, comoving_OmM

!         real(rt) :: frac

!         frac = comoving_OmB / comoving_OmM

!       end subroutine fort_get_omb


! ! :::
! ! ::: ----------------------------------------------------------------
! ! :::

!       subroutine fort_get_omm(omm) &
!          bind(C, name="fort_get_omm")

!         use comoving_module, only: comoving_OmM

!         real(rt) :: omm

!         omm = comoving_OmM

!       end subroutine fort_get_omm

! ! :::
! ! ::: ----------------------------------------------------------------
! ! :::

!       subroutine fort_get_hubble(hubble) &
!          bind(C, name="fort_get_hubble")

!         use comoving_module, only: comoving_h

!         real(rt) :: hubble

!         hubble = comoving_h

!       end subroutine fort_get_hubble


! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_set_omb(omb) &
         bind(C, name="fort_set_omb")

        use comoving_module, only: comoving_OmB

        real(rt), intent(in) :: omb

        comoving_OmB = omb

      end subroutine fort_set_omb


! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_set_omm(omm) &
         bind(C, name="fort_set_omm")

        use comoving_module, only: comoving_OmM

        real(rt), intent(in) :: omm

        comoving_OmM = omm

      end subroutine fort_set_omm

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_set_hubble(hubble) &
         bind(C, name="fort_set_hubble")

        use comoving_module, only: comoving_h

        real(rt), intent(in) :: hubble

        comoving_h = hubble

      end subroutine fort_set_hubble

end module comoving_nd_module
