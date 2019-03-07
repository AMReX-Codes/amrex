module amrex_particle_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real, amrex_particle_real

  implicit none

  private

  public :: amrex_particle_set_position, amrex_particle_get_position, &
       amrex_deposit_cic, amrex_interpolate_cic

contains

  subroutine amrex_particle_set_position (particles, ns, np, x, y, z) &
       bind(c,name='amrex_particle_set_position')
    integer(c_int)  ,          intent(in   ), value :: ns, np
    real(amrex_particle_real), intent(inout)        :: particles(ns,np)
    real(amrex_real),          intent(in   )        :: x(np), y(np), z(np)

    integer :: i

    do i = 1, np
       particles(1,i) = x(i)
       particles(2,i) = y(i)
       particles(3,i) = z(i)
    end do
  end subroutine amrex_particle_set_position

  subroutine amrex_particle_get_position (particles, ns, np, x, y, z) &
       bind(c,name='amrex_particle_get_position')
    integer(c_int)  ,          intent(in   ), value :: ns, np
    real(amrex_particle_real), intent(in   )        :: particles(ns,np)
    real(amrex_real),          intent(  out)        :: x(np), y(np), z(np)

    integer :: i

    do i = 1, np
       x(i) = particles(1,i)
       y(i) = particles(2,i)
       z(i) = particles(3,i)
    end do
  end subroutine amrex_particle_get_position

  subroutine amrex_deposit_cic(particles, ns, np, nc, rho, lo, hi, plo, dx) &
       bind(c,name='amrex_deposit_cic')
    integer, value                :: ns, np, nc
    real(amrex_particle_real)     :: particles(ns,np)
    integer                       :: lo(3)
    integer                       :: hi(3)
    real(amrex_real)              :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),nc)
    real(amrex_real)              :: plo(3)
    real(amrex_real)              :: dx(3)

    integer i, j, k, n, comp
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) inv_dx(3)
    inv_dx = 1.0d0/dx

    do n = 1, np
       lx = (particles(1, n) - plo(1))*inv_dx(1) + 0.5d0
       ly = (particles(2, n) - plo(2))*inv_dx(2) + 0.5d0
       lz = (particles(3, n) - plo(3))*inv_dx(3) + 0.5d0

       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       wx_hi = lx - i
       wy_hi = ly - j
       wz_hi = lz - k

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       rho(i-1, j-1, k-1, 1) = rho(i-1, j-1, k-1, 1) + wx_lo*wy_lo*wz_lo*particles(4, n)
       rho(i-1, j-1, k  , 1) = rho(i-1, j-1, k  , 1) + wx_lo*wy_lo*wz_hi*particles(4, n)
       rho(i-1, j,   k-1, 1) = rho(i-1, j,   k-1, 1) + wx_lo*wy_hi*wz_lo*particles(4, n)
       rho(i-1, j,   k  , 1) = rho(i-1, j,   k,   1) + wx_lo*wy_hi*wz_hi*particles(4, n)
       rho(i,   j-1, k-1, 1) = rho(i,   j-1, k-1, 1) + wx_hi*wy_lo*wz_lo*particles(4, n)
       rho(i,   j-1, k  , 1) = rho(i,   j-1, k  , 1) + wx_hi*wy_lo*wz_hi*particles(4, n)
       rho(i,   j,   k-1, 1) = rho(i,   j,   k-1, 1) + wx_hi*wy_hi*wz_lo*particles(4, n)
       rho(i,   j,   k  , 1) = rho(i,   j,   k  , 1) + wx_hi*wy_hi*wz_hi*particles(4, n)

       do comp = 2, nc
          rho(i-1, j-1, k-1, comp) = rho(i-1, j-1, k-1, comp) + wx_lo*wy_lo*wz_lo*particles(4, n)*particles(3+comp, n)
          rho(i-1, j-1, k  , comp) = rho(i-1, j-1, k  , comp) + wx_lo*wy_lo*wz_hi*particles(4, n)*particles(3+comp, n)
          rho(i-1, j,   k-1, comp) = rho(i-1, j,   k-1, comp) + wx_lo*wy_hi*wz_lo*particles(4, n)*particles(3+comp, n)
          rho(i-1, j,   k  , comp) = rho(i-1, j,   k,   comp) + wx_lo*wy_hi*wz_hi*particles(4, n)*particles(3+comp, n)
          rho(i,   j-1, k-1, comp) = rho(i,   j-1, k-1, comp) + wx_hi*wy_lo*wz_lo*particles(4, n)*particles(3+comp, n)
          rho(i,   j-1, k  , comp) = rho(i,   j-1, k  , comp) + wx_hi*wy_lo*wz_hi*particles(4, n)*particles(3+comp, n)
          rho(i,   j,   k-1, comp) = rho(i,   j,   k-1, comp) + wx_hi*wy_hi*wz_lo*particles(4, n)*particles(3+comp, n)
          rho(i,   j,   k  , comp) = rho(i,   j,   k  , comp) + wx_hi*wy_hi*wz_hi*particles(4, n)*particles(3+comp, n)
       end do

    end do

  end subroutine amrex_deposit_cic

  subroutine amrex_deposit_particle_dx_cic(particles, ns, np, nc, & 
                                           rho, lo, hi, plo, dx,  &
                                           dx_particle) &
       bind(c,name='amrex_deposit_particle_dx_cic')
    integer, value                :: ns, np, nc
    real(amrex_particle_real)     :: particles(ns,np)
    integer                       :: lo(3)
    integer                       :: hi(3)
    real(amrex_real)              :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),nc)
    real(amrex_real)              :: plo(3)
    real(amrex_real)              :: dx(3)
    real(amrex_real)              :: dx_particle(3)

    integer i, j, k, n, comp
    real(amrex_real) lx, ly, lz, hx, hy, hz
    integer lo_x, lo_y, lo_z, hi_x, hi_y, hi_z
    real(amrex_real) wx, wy, wz
    real(amrex_real) inv_dx(3)
    real (amrex_real) factor, weight

    factor = (dx(1)/dx_particle(1))*(dx(2)/dx_particle(2))*(dx(3)/dx_particle(3))
    inv_dx = 1.0d0/dx

    do n = 1, np

       lx = (particles(1, n) - plo(1) - 0.5d0*dx_particle(1))*inv_dx(1)
       ly = (particles(2, n) - plo(2) - 0.5d0*dx_particle(2))*inv_dx(2)
       lz = (particles(3, n) - plo(3) - 0.5d0*dx_particle(3))*inv_dx(3)

       hx = (particles(1, n) - plo(1) + 0.5d0*dx_particle(1))*inv_dx(1)
       hy = (particles(2, n) - plo(2) + 0.5d0*dx_particle(2))*inv_dx(2)
       hz = (particles(3, n) - plo(3) + 0.5d0*dx_particle(3))*inv_dx(3)

       lo_x = floor(lx)
       lo_y = floor(ly)
       lo_z = floor(lz)

       hi_x = floor(hx)
       hi_y = floor(hy)
       hi_z = floor(hz)

       do i = lo_x, hi_x
          if (i < lo(1) .or. i > hi(1)) then
             cycle
          end if
          wx = min(hx - i, 1.d0) - max(lx - i, 0.d0)
          do j = lo_y, hi_y
             if (j < lo(2) .or. j > hi(2)) then
                cycle
             end if
             wy = min(hy - j, 1.d0) - max(ly - j, 0.d0)
             do k = lo_z, hi_z
                if (k < lo(3) .or. k > hi(3)) then
                   cycle
                end if
                wz = min(hz - k, 1.d0) - max(lz - k, 0.d0)

                weight = wx*wy*wz*factor

                rho(i, j, k, 1) = rho(i, j, k, 1) + weight*particles(4, n)

                do comp = 2, nc
                   rho(i, j, k, comp) = rho(i, j, k, comp) + weight*particles(4, n)*particles(3+comp, n) 
                end do
             end do
          end do
       end do
    end do

  end subroutine amrex_deposit_particle_dx_cic

  subroutine amrex_interpolate_cic(particles, ns, np, acc, lo, hi, ncomp, plo, dx) &
       bind(c,name='amrex_interpolate_cic')
    integer, value                :: ns, np, ncomp
    real(amrex_particle_real)     :: particles(ns,np)
    integer                       :: lo(3)
    integer                       :: hi(3)
    real(amrex_real)              :: acc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), ncomp)
    real(amrex_real)              :: plo(3)
    real(amrex_real)              :: dx(3)
    real(amrex_real)              :: acceleration(ncomp)

    integer i, j, k, n, nc
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) inv_dx(3)
    inv_dx = 1.0d0/dx

    do n = 1, np
       lx = (particles(1, n) - plo(1))*inv_dx(1) + 0.5d0
       ly = (particles(2, n) - plo(2))*inv_dx(2) + 0.5d0
       lz = (particles(3, n) - plo(3))*inv_dx(3) + 0.5d0

       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       wx_hi = lx - i
       wy_hi = ly - j
       wz_hi = lz - k

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       do nc = 1, ncomp
          acceleration(nc) = wx_lo*wy_lo*wz_lo*acc(i-1, j-1, k-1, nc) + &
                             wx_lo*wy_lo*wz_hi*acc(i-1, j-1, k  , nc) + &
                             wx_lo*wy_hi*wz_lo*acc(i-1, j,   k-1, nc) + &
                             wx_lo*wy_hi*wz_hi*acc(i-1, j,   k  , nc) + &
                             wx_hi*wy_lo*wz_lo*acc(i,   j-1, k-1, nc) + &
                             wx_hi*wy_lo*wz_hi*acc(i,   j-1, k  , nc) + &
                             wx_hi*wy_hi*wz_lo*acc(i,   j,   k-1, nc) + &
                             wx_hi*wy_hi*wz_hi*acc(i,   j,   k  , nc)
       
#ifdef AMREX_DEBUG
!          if (abs(acceleration(nc) - 5.d0) .ge. 1.0d-9) then
!             print *, particles(1, n), particles(2, n), particles(3, n)
!          end if
#endif

       end do
    end do

  end subroutine amrex_interpolate_cic

end module amrex_particle_module
