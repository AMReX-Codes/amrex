  subroutine amrex_move_particles(particles, ns, np, dt, prob_lo, prob_hi) &
       bind(c,name='amrex_move_particles')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    integer,          intent(in   ), value :: ns, np
    real(amrex_real), intent(inout)        :: particles(ns, np)
    real(amrex_real), intent(in   )        :: dt
    real(amrex_real), intent(in   )        :: prob_lo(3), prob_hi(3)

    integer i

    do i = 1, np

!      update the particle positions / velocites
       particles(4, i) = particles(4, i) + particles(7, i) * dt 
       particles(5, i) = particles(5, i) + particles(8, i) * dt 
       particles(6, i) = particles(6, i) + particles(9, i) * dt 


       particles(1, i) = particles(1, i) + particles(4, i) * dt 
       particles(2, i) = particles(2, i) + particles(5, i) * dt 
       particles(3, i) = particles(3, i) + particles(6, i) * dt 

!      bounce off the walls in the x...
       do while (particles(1, i) .lt. prob_lo(1) .or. particles(1, i) .gt. prob_hi(1))
          if (particles(1, i) .lt. prob_lo(1)) then
             particles(1, i) = 2.d0*prob_lo(1) - particles(1, i)
          else
             particles(1, i) = 2.d0*prob_hi(1) - particles(1, i)
          end if
          particles(4, i) = -particles(4, i)
       end do

!      ... y... 
       do while (particles(2, i) .lt. prob_lo(2) .or. particles(2, i) .gt. prob_hi(2))
          if (particles(2, i) .lt. prob_lo(2)) then
             particles(2, i) = 2.d0*prob_lo(2) - particles(2, i)
          else
             particles(2, i) = 2.d0*prob_hi(2) - particles(2, i)
          end if
          particles(5, i) = -particles(5, i)
       end do

!      ... and z directions
       do while (particles(3, i) .lt. prob_lo(3) .or. particles(3, i) .gt. prob_hi(3))
          if (particles(3, i) .lt. prob_lo(3)) then
             particles(3, i) = 2.d0*prob_lo(3) - particles(3, i)
          else
             particles(3, i) = 2.d0*prob_hi(3) - particles(3, i)
          end if
          particles(6, i) = -particles(6, i)
       end do

    end do

  end subroutine amrex_move_particles


  subroutine amrex_compute_forces(particles, ns, np, ghosts, ng) &
       bind(c,name='amrex_compute_forces')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    integer,          intent(in   ), value :: ns, np, ng
    real(amrex_real), intent(inout)        :: particles(ns, np)
    real(amrex_real), intent(in   )        :: ghosts(3, ng)

    real(amrex_real) dx, dy, dz, r2, r, coef
    real(amrex_real) cutoff, min_r, mass
    integer i, j

    cutoff = 1.d-2
    min_r  = (cutoff/1.d2)
    mass   = 1.d-2
    
    do i = 1, np

!      zero out the particle acceleration
       particles(7, i) = 0.d0
       particles(8, i) = 0.d0
       particles(9, i) = 0.d0

       do j = 1, np

          if (i .eq. j ) then
             cycle
          end if

          dx = particles(1, i) - particles(1, j)
          dy = particles(2, i) - particles(2, j)
          dz = particles(3, i) - particles(3, j)

          r2 = dx * dx + dy * dy + dz * dz

          if (r2 .gt. cutoff*cutoff) then
             cycle
          end if

          r2 = max(r2, min_r*min_r) 
          r = sqrt(r2)

          coef = (1.d0 - cutoff / r) / r2 / mass
          particles(7, i) = particles(7, i) + coef * dx
          particles(8, i) = particles(8, i) + coef * dy
          particles(9, i) = particles(9, i) + coef * dz

       end do

       do j = 1, ng

          dx = particles(1, i) - ghosts(1, j)
          dy = particles(2, i) - ghosts(2, j)
          dz = particles(3, i) - ghosts(3, j)

          r2 = dx * dx + dy * dy + dz * dz

          if (r2 .gt. cutoff*cutoff) then
             cycle
          end if

          r2 = max(r2, min_r*min_r) 
          r = sqrt(r2)

          coef = (1.d0 - cutoff / r) / r2 / mass
          particles(7, i) = particles(7, i) + coef * dx
          particles(8, i) = particles(8, i) + coef * dy
          particles(9, i) = particles(9, i) + coef * dz
          
      end do
    end do

  end subroutine amrex_compute_forces
