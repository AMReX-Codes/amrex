subroutine push_particles(particles, ns, np) &
     bind(c,name='push_particles')
  
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  
  implicit none

  integer              :: ns, np
  real(amrex_real)     :: particles(ns, np)
  integer i

  real(amrex_real) mass, charge, dt, fac
  mass = 1.d0
  charge = 1.d0
  dt = 1.d-6
  fac = charge*dt / mass

  do n = 1, np

     particles(5, n) = particles(5, n) + fac*particles(8,  n)
     particles(6, n) = particles(6, n) + fac*particles(9,  n)
     particles(7, n) = particles(7, n) + fac*particles(10, n)

     particles(1, n) = particles(1, n) + dt * particles(5, n)
     particles(2, n) = particles(2, n) + dt * particles(6, n)
     particles(3, n) = particles(3, n) + dt * particles(7, n)

  end if

end subroutine push_particles

subroutine deposit_cic(particles, ns, np, &
     counts, offsets, ngrids, gid, &
     rho, lo, hi, plo, dx) &
     bind(c,name='deposit_cic')
  
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  
  implicit none
  
  integer, value       :: ns, np
  real(amrex_real)     :: particles(ns,np)
  integer              :: lo(3)
  integer              :: hi(3)
  real(amrex_real)     :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
  real(amrex_real)     :: plo(3)
  real(amrex_real)     :: dx(3)

  integer i, j, k, n
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
     
     rho(i-1, j-1, k-1) = rho(i-1, j-1, k-1) + wx_lo*wy_lo*wz_lo*particles(4, n)
     rho(i-1, j-1, k)   = rho(i-1, j-1, k)   + wx_lo*wy_lo*wz_hi*particles(4, n)
     rho(i-1, j,   k-1) = rho(i-1, j,   k-1) + wx_lo*wy_hi*wz_lo*particles(4, n)
     rho(i-1, j,   k)   = rho(i-1, j,   k)   + wx_lo*wy_hi*wz_hi*particles(4, n)
     rho(i,   j-1, k-1) = rho(i,   j-1, k-1) + wx_hi*wy_lo*wz_lo*particles(4, n)
     rho(i,   j-1, k)   = rho(i,   j-1, k)   + wx_hi*wy_lo*wz_hi*particles(4, n)
     rho(i,   j,   k-1) = rho(i,   j,   k-1) + wx_hi*wy_hi*wz_lo*particles(4, n)
     rho(i,   j,   k)   = rho(i,   j,   k)   + wx_hi*wy_hi*wz_hi*particles(4, n)
     
  end do

end subroutine deposit_cic

subroutine interpolate_cic(particles, ns, np, &
     counts, offsets, ngrids, gid, & 
     acc, lo, hi, plo, dx) &
     bind(c,name='interpolate_cic')
  
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  
  implicit none

  integer, value       :: ns, np
  real(amrex_real)     :: particles(ns,np)
  integer              :: lo(3)
  integer              :: hi(3)
  real(amrex_real)     :: acc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 3)
  real(amrex_real)     :: plo(3)
  real(amrex_real)     :: dx(3)
  
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
   
     do nc = 1, 3
        particles(7+nc, n)  = wx_lo*wy_lo*wz_lo*acc(i-1, j-1, k-1, nc) + &
             wx_lo*wy_lo*wz_hi*acc(i-1, j-1, k  , nc) + &
             wx_lo*wy_hi*wz_lo*acc(i-1, j,   k-1, nc) + &
             wx_lo*wy_hi*wz_hi*acc(i-1, j,   k  , nc) + &
             wx_hi*wy_lo*wz_lo*acc(i,   j-1, k-1, nc) + &
             wx_hi*wy_lo*wz_hi*acc(i,   j-1, k  , nc) + &
             wx_hi*wy_hi*wz_lo*acc(i,   j,   k-1, nc) + &
             wx_hi*wy_hi*wz_hi*acc(i,   j,   k  , nc)       
     end do
  end do

end subroutine interpolate_cic
