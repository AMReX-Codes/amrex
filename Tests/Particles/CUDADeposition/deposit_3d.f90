subroutine deposit_cic(particles, ns, np, rho, lo, hi, plo, dx) &
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
