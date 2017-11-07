module amrex_particle_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real, amrex_particle_real

  implicit none

  private

  public :: amrex_particle_set_position, amrex_particle_get_position

contains

  subroutine amrex_particle_set_position (particles, ns, np, x) &
       bind(c,name='amrex_particle_set_position')
    integer(c_int),            intent(in   ), value :: ns, np
    real(amrex_particle_real), intent(inout)        :: particles(ns,np)
    real(amrex_real),          intent(in   )        :: x(np)

    integer :: i

    do i = 1, np
       particles(1,i) = x(i)
    end do
  end subroutine amrex_particle_set_position

  subroutine amrex_particle_get_position (particles, ns, np, x) &
       bind(c,name='amrex_particle_get_position')
    integer(c_int)  ,          intent(in   ), value :: ns, np
    real(amrex_particle_real), intent(in   )        :: particles(ns,np)
    real(amrex_real),          intent(  out)        :: x(np)

    integer :: i

    do i = 1, np
       x(i) = particles(1,i)
    end do
  end subroutine amrex_particle_get_position

  subroutine amrex_deposit_cic(particles, ns, np, nc, rho, lo, hi, plo, dx) &
       bind(c,name='amrex_deposit_cic')
    integer, value                :: ns, np, nc
    real(amrex_particle_real)     :: particles(ns,np)
    integer                       :: lo(1)
    integer                       :: hi(1)
    real(amrex_real)              :: rho(lo(1):hi(1), nc)
    real(amrex_real)              :: plo(1)
    real(amrex_real)              :: dx(1)

    integer i, n, comp
    real(amrex_real) wx_lo, wx_hi
    real(amrex_real) lx
    real(amrex_real) inv_dx(1)
    inv_dx = 1.0d0/dx

    do n = 1, np
       lx = (particles(1, n) - plo(1))*inv_dx(1) + 0.5d0
       i = floor(lx)
       wx_hi = lx - i
       wx_lo = 1.0d0 - wx_hi

       rho(i-1, 1) = rho(i-1, 1) + wx_lo*particles(2, n)
       rho(i  , 1) = rho(i  , 1) + wx_hi*particles(2, n)

       do comp = 2, nc
          rho(i-1, comp) = rho(i-1, comp) + wx_lo*particles(2,n)*particles(1 + comp, n)
          rho(i  , comp) = rho(i  , comp) + wx_hi*particles(2,n)*particles(1 + comp, n)
       end do
    end do

  end subroutine amrex_deposit_cic

  subroutine amrex_deposit_particle_dx_cic(particles, ns, np, nc, & 
                                           rho, lo, hi, plo, dx,  &
                                           dx_particle) &
       bind(c,name='amrex_deposit_particle_dx_cic')
    integer, value                :: ns, np, nc
    real(amrex_particle_real)     :: particles(ns,np)
    integer                       :: lo(1)
    integer                       :: hi(1)
    real(amrex_real)              :: rho(lo(1):hi(1), nc)
    real(amrex_real)              :: plo(1)
    real(amrex_real)              :: dx(1)
    real(amrex_real)              :: dx_particle(1)

    integer i, n, comp
    real(amrex_real) lx, hx
    integer lo_x, hi_x
    real(amrex_real) wx
    real(amrex_real) inv_dx(1)
    real (amrex_real) factor, weight

    factor = dx(1)/dx_particle(1)
    inv_dx = 1.0d0/dx

    do n = 1, np

       lx = (particles(1, n) - plo(1) - 0.5d0*dx_particle(1))*inv_dx(1)
       hx = (particles(1, n) - plo(1) + 0.5d0*dx_particle(1))*inv_dx(1)

       lo_x = floor(lx)
       hi_x = floor(hx)

       do i = lo_x, hi_x
          if (i < lo(1) .or. i > hi(1)) then
             cycle
          end if
          wx = min(hx - i, 1.d0) - max(lx - i, 0.d0)

          weight = wx*factor

          rho(i, 1) = rho(i, 1) + weight*particles(2, n)

          do comp = 2, nc
             rho(i, comp) = rho(i, comp) + weight*particles(2, n)*particles(1+comp, n) 
          end do

       end do
    end do

  end subroutine amrex_deposit_particle_dx_cic

  subroutine amrex_interpolate_cic(particles, ns, np, acc, lo, hi, ncomp, plo, dx) &
       bind(c,name='amrex_interpolate_cic')
    integer, value                :: ns, np, ncomp
    real(amrex_particle_real)     :: particles(ns,np)
    integer                       :: lo(1)
    integer                       :: hi(1)
    real(amrex_real)              :: acc(lo(1):hi(1), ncomp)
    real(amrex_real)              :: plo(1)
    real(amrex_real)              :: dx(1)
    real(amrex_real)              :: acceleration(ncomp)

    integer i, n, nc
    real(amrex_real) wx_lo, wx_hi
    real(amrex_real) lx
    real(amrex_real) inv_dx(1)
    inv_dx = 1.0d0/dx

    do n = 1, np
       lx = (particles(1, n) - plo(1))*inv_dx(1) + 0.5d0
       i = floor(lx)
       wx_hi = lx - i
       wx_lo = 1.0d0 - wx_hi

       do nc = 1, ncomp
          acceleration(nc) = wx_lo*acc(i-1, nc) + &
                             wx_hi*acc(i,   nc)
       
          if (abs(acceleration(nc) - 5.d0) .ge. 1.0d-9) then
             print *, particles(1, n)
          end if

       end do
    end do

  end subroutine amrex_interpolate_cic

end module amrex_particle_module
