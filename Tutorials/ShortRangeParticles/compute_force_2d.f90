  subroutine amrex_compute_forces(particles, ns, np, ghosts, ng, ax, ay) &
       bind(c,name='amrex_compute_forces')

    use iso_c_binding
    use amrex_fort_module, only : amrex_real
    integer,          intent(in   ), value :: ns, np, ng
    real(amrex_real), intent(in   )        :: particles(ns, np)
    real(amrex_real), intent(in   )        :: ghosts(2, ng)
    real(amrex_real), intent(  out)        :: ax(np)
    real(amrex_real), intent(  out)        :: ay(np)

    real(amrex_real) dx, dy, r2, r, coef
    real(amrex_real) cutoff, min_r, mass
    integer i, j

    cutoff = 1.d-2
    min_r  = (cutoff/1.d2)
    mass   = 1.d-2
    
    do i = 1, np

       ax(i) = 0.d0
       ay(i) = 0.d0

       do j = 1, np

          if (i .eq. j ) then
             cycle
          end if

          dx = particles(1, i) - particles(1, j)
          dy = particles(2, i) - particles(2, j)

          r2 = dx * dx + dy * dy

          if (r2 .gt. cutoff*cutoff) then
             cycle
          end if

          r2 = max(r2, min_r*min_r) 
          r = sqrt(r2)

          coef = (1.d0 - cutoff / r) / r2 / mass
          ax(i) = ax(i) + coef * dx
          ay(i) = ay(i) + coef * dx

       end do

       do j = 1, ng

          dx = particles(1, i) - ghosts(1, j)
          dy = particles(2, i) - ghosts(2, j)

          r2 = dx * dx + dy * dy

          if (r2 .gt. cutoff*cutoff) then
             cycle
          end if

          r2 = max(r2, min_r*min_r) 
          r = sqrt(r2)

          coef = (1.d0 - cutoff / r) / r2 / mass
          ax(i) = ax(i) + coef * dx
          ay(i) = ay(i) + coef * dx
          
      end do
    end do

  end subroutine amrex_compute_forces
