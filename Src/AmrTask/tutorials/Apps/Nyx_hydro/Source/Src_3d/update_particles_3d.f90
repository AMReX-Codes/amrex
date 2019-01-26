
  subroutine update_dm_particles(np, particles, accel, accel_lo, accel_hi, &
                                 plo, dx, dt, a_prev, a_cur, do_move) &
       bind(c,name='update_dm_particles')

    use iso_c_binding
    use amrex_error_module
    use amrex_fort_module, only : amrex_real
    use particle_mod      , only: dm_particle_t

    integer,             intent(in   )        :: np
    type(dm_particle_t), intent(inout)        :: particles(np)
    integer,             intent(in   )        :: accel_lo(3), accel_hi(3)
    real(amrex_real),    intent(in   )        :: accel &
         (accel_lo(1):accel_hi(1),accel_lo(2):accel_hi(2),accel_lo(3):accel_hi(3),3)
    real(amrex_real),    intent(in   )        :: plo(3),dx(3),dt,a_prev,a_cur
    integer,             intent(in   )        :: do_move

    integer          :: i, j, k, n, nc
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) half_dt, a_cur_inv, dt_a_cur_inv
    real(amrex_real) inv_dx(3)
    real(amrex_real) part_accel

    inv_dx = 1.0d0/dx

    half_dt       = 0.5d0 * dt
    a_cur_inv    = 1.d0 / a_cur;
    dt_a_cur_inv = dt * a_cur_inv;

    do n = 1, np

       if (particles(n)%id .le. 0) then
          cycle
       end if

       lx = (particles(n)%pos(1) - plo(1))*inv_dx(1) + 0.5d0
       ly = (particles(n)%pos(2) - plo(2))*inv_dx(2) + 0.5d0
       lz = (particles(n)%pos(3) - plo(3))*inv_dx(3) + 0.5d0

       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       if (i-1 .lt. accel_lo(1) .or. i .gt. accel_hi(1) .or. &
           j-1 .lt. accel_lo(2) .or. j .gt. accel_hi(2) .or. &
           k-1 .lt. accel_lo(3) .or. k .gt. accel_hi(3)) then
          print *,'PARTICLE ID ', particles(n)%id,' REACHING OUT OF BOUNDS AT (I,J,K) = ',i,j,k
          call amrex_error('Aborting in move_kick_drift')
       end if

       wx_hi = lx - i
       wy_hi = ly - j
       wz_hi = lz - k

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       ! moveKickDrift: update velocity by half dt:  (a u)^half = (a u)^old  + dt/2 grav^old
       ! moveKick    t: update velocity by half dt:  (a u)^new  = (a u)^half + dt/2 grav^new
       do nc = 1, 3

          part_accel = &
                wx_lo*wy_lo*wz_lo*accel(i-1, j-1, k-1, nc) + &
                wx_lo*wy_lo*wz_hi*accel(i-1, j-1, k  , nc) + &
                wx_lo*wy_hi*wz_lo*accel(i-1, j,   k-1, nc) + &
                wx_lo*wy_hi*wz_hi*accel(i-1, j,   k  , nc) + &
                wx_hi*wy_lo*wz_lo*accel(i,   j-1, k-1, nc) + &
                wx_hi*wy_lo*wz_hi*accel(i,   j-1, k  , nc) + &
                wx_hi*wy_hi*wz_lo*accel(i,   j,   k-1, nc) + &
                wx_hi*wy_hi*wz_hi*accel(i,   j,   k  , nc)

          particles(n)%vel(nc) = a_prev*particles(n)%vel(nc) + half_dt * part_accel
          particles(n)%vel(nc) = particles(n)%vel(nc) * a_cur_inv

       end do

       ! moveKickDrift: Update position by full dt: x^new = x^old + dt u^half / a^half
       if (do_move .eq. 1) then
          do nc = 1, 3
             particles(n)%pos(nc) = particles(n)%pos(nc) + dt_a_cur_inv * particles(n)%vel(nc)
          end do
       end if

    end do

  end subroutine update_dm_particles

  subroutine update_agn_particles(np, particles, accel, accel_lo, accel_hi, &
                                  plo, dx, dt, a_prev, a_cur, do_move) &
       bind(c,name='update_agn_particles')

    use iso_c_binding
    use amrex_error_module
    use amrex_fort_module, only : amrex_real
    use particle_mod      , only: agn_particle_t

    integer,             intent(in   )        :: np
    type(agn_particle_t), intent(inout)        :: particles(np)
    integer,             intent(in   )        :: accel_lo(3), accel_hi(3)
    real(amrex_real),    intent(in   )        :: accel &
         (accel_lo(1):accel_hi(1),accel_lo(2):accel_hi(2),accel_lo(3):accel_hi(3),3)
    real(amrex_real),    intent(in   )        :: plo(3),dx(3),dt,a_prev,a_cur
    integer,             intent(in   )        :: do_move

    integer          :: i, j, k, n, nc
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) half_dt, a_cur_inv, dt_a_cur_inv
    real(amrex_real) inv_dx(3)
    real(amrex_real) part_accel

    inv_dx = 1.0d0/dx

    half_dt       = 0.5d0 * dt
    a_cur_inv    = 1.d0 / a_cur;
    dt_a_cur_inv = dt * a_cur_inv;

    do n = 1, np

       if (particles(n)%id .le. 0) then
          cycle
       end if

       lx = (particles(n)%pos(1) - plo(1))*inv_dx(1) + 0.5d0
       ly = (particles(n)%pos(2) - plo(2))*inv_dx(2) + 0.5d0
       lz = (particles(n)%pos(3) - plo(3))*inv_dx(3) + 0.5d0

       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       if (i-1 .lt. accel_lo(1) .or. i .gt. accel_hi(1) .or. &
           j-1 .lt. accel_lo(2) .or. j .gt. accel_hi(2) .or. &
           k-1 .lt. accel_lo(3) .or. k .gt. accel_hi(3)) then
          print *,'PARTICLE ID ', particles(n)%id,' REACHING OUT OF BOUNDS AT (I,J,K) = ',i,j,k
          call amrex_error('Aborting in move_kick_drift')
       end if

       wx_hi = lx - i
       wy_hi = ly - j
       wz_hi = lz - k

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       ! moveKickDrift: update velocity by half dt:  (a u)^half = (a u)^old  + dt/2 grav^old
       ! moveKick    t: update velocity by half dt:  (a u)^new  = (a u)^half + dt/2 grav^new
       do nc = 1, 3

          part_accel = &
                wx_lo*wy_lo*wz_lo*accel(i-1, j-1, k-1, nc) + &
                wx_lo*wy_lo*wz_hi*accel(i-1, j-1, k  , nc) + &
                wx_lo*wy_hi*wz_lo*accel(i-1, j,   k-1, nc) + &
                wx_lo*wy_hi*wz_hi*accel(i-1, j,   k  , nc) + &
                wx_hi*wy_lo*wz_lo*accel(i,   j-1, k-1, nc) + &
                wx_hi*wy_lo*wz_hi*accel(i,   j-1, k  , nc) + &
                wx_hi*wy_hi*wz_lo*accel(i,   j,   k-1, nc) + &
                wx_hi*wy_hi*wz_hi*accel(i,   j,   k  , nc)

          particles(n)%vel(nc) = a_prev*particles(n)%vel(nc) + half_dt * part_accel
          particles(n)%vel(nc) = particles(n)%vel(nc) * a_cur_inv

       end do

       ! moveKickDrift: Update position by full dt: x^new = x^old + dt u^half / a^half
       if (do_move .eq. 1) then
          do nc = 1, 3
             particles(n)%pos(nc) = particles(n)%pos(nc) + dt_a_cur_inv * particles(n)%vel(nc)
          end do
       end if

    end do

  end subroutine update_agn_particles

