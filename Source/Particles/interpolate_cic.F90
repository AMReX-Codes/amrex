! Copyright 2019 Maxence Thevenet, Weiqun Zhang
!
! This file is part of WarpX.
!
! License: BSD-3-Clause-LBNL

module warpx_ES_interpolate_cic

  use iso_c_binding
  use amrex_fort_module, only : amrex_real, amrex_particle_real

  implicit none

contains

! This routine interpolates the electric field to the particle positions
! using cloud-in-cell interpolation. The electric fields are assumed to be
! node-centered.
!
! Arguments:
!     particles : a pointer to the particle array-of-structs
!     ns        : the stride length of particle struct (the size of the struct in number of reals)
!     np        : the number of particles
!     Ex_p      : the electric field in the x-direction at the particle positions (output)
!     Ey_p      : the electric field in the y-direction at the particle positions (output)
!     Ez_p      : the electric field in the z-direction at the particle positions (output)
!     Ex, Ey, Ez: Fabs conting the electric field on the mesh
!     lo        : a pointer to the lo corner of this valid box, in index space
!     hi        : a pointer to the hi corner of this valid box, in index space
!     plo       : the real position of the left-hand corner of the problem domain
!     dx        : the mesh spacing
!     ng        : the number of ghost cells for the E-field
!
  subroutine warpx_interpolate_cic_3d(particles, ns, np,      &
                                      Ex_p, Ey_p, Ez_p,       &
                                      Ex,   Ey,   Ez,         &
                                      lo, hi, plo, dx, ng)    &
       bind(c,name='warpx_interpolate_cic_3d')
    integer, value,   intent(in)     :: ns, np
    real(amrex_particle_real), intent(in)     :: particles(ns,np)
    real(amrex_real), intent(inout)  :: Ex_p(np), Ey_p(np), Ez_p(np)
    integer,          intent(in)     :: ng
    integer,          intent(in)     :: lo(3)
    integer,          intent(in)     :: hi(3)
    real(amrex_real), intent(in)     :: Ex(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng)
    real(amrex_real), intent(in)     :: Ey(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng)
    real(amrex_real), intent(in)     :: Ez(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng)
    real(amrex_real), intent(in)     :: plo(3)
    real(amrex_real), intent(in)     :: dx(3)

    integer i, j, k, n
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) inv_dx(3)
    inv_dx = 1.0d0/dx

    do n = 1, np
       lx = (particles(1, n) - plo(1))*inv_dx(1)
       ly = (particles(2, n) - plo(2))*inv_dx(2)
       lz = (particles(3, n) - plo(3))*inv_dx(3)

       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       wx_hi = lx - i
       wy_hi = ly - j
       wz_hi = lz - k

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi
       wz_lo = 1.0d0 - wz_hi

       Ex_p(n) = wx_lo*wy_lo*wz_lo*Ex(i,   j,   k  ) + &
                 wx_lo*wy_lo*wz_hi*Ex(i,   j,   k+1) + &
                 wx_lo*wy_hi*wz_lo*Ex(i,   j+1, k  ) + &
                 wx_lo*wy_hi*wz_hi*Ex(i,   j+1, k+1) + &
                 wx_hi*wy_lo*wz_lo*Ex(i+1, j,   k  ) + &
                 wx_hi*wy_lo*wz_hi*Ex(i+1, j,   k+1) + &
                 wx_hi*wy_hi*wz_lo*Ex(i+1, j+1, k  ) + &
                 wx_hi*wy_hi*wz_hi*Ex(i+1, j+1, k+1)

       Ey_p(n) = wx_lo*wy_lo*wz_lo*Ey(i,   j,   k  ) + &
                 wx_lo*wy_lo*wz_hi*Ey(i,   j,   k+1) + &
                 wx_lo*wy_hi*wz_lo*Ey(i,   j+1, k  ) + &
                 wx_lo*wy_hi*wz_hi*Ey(i,   j+1, k+1) + &
                 wx_hi*wy_lo*wz_lo*Ey(i+1, j,   k  ) + &
                 wx_hi*wy_lo*wz_hi*Ey(i+1, j,   k+1) + &
                 wx_hi*wy_hi*wz_lo*Ey(i+1, j+1, k  ) + &
                 wx_hi*wy_hi*wz_hi*Ey(i+1, j+1, k+1)

       Ez_p(n) = wx_lo*wy_lo*wz_lo*Ez(i,   j,   k  ) + &
                 wx_lo*wy_lo*wz_hi*Ez(i,   j,   k+1) + &
                 wx_lo*wy_hi*wz_lo*Ez(i,   j+1, k  ) + &
                 wx_lo*wy_hi*wz_hi*Ez(i,   j+1, k+1) + &
                 wx_hi*wy_lo*wz_lo*Ez(i+1, j,   k  ) + &
                 wx_hi*wy_lo*wz_hi*Ez(i+1, j,   k+1) + &
                 wx_hi*wy_hi*wz_lo*Ez(i+1, j+1, k  ) + &
                 wx_hi*wy_hi*wz_hi*Ez(i+1, j+1, k+1)

    end do

  end subroutine warpx_interpolate_cic_3d


  subroutine warpx_interpolate_cic_2d(particles, ns, np,      &
                             Ex_p, Ey_p,             &
                             Ex,   Ey,               &
                             lo, hi, plo, dx, ng)    &
       bind(c,name='warpx_interpolate_cic_2d')
    integer, value,   intent(in)     :: ns, np
    real(amrex_particle_real), intent(in)     :: particles(ns,np)
    real(amrex_real), intent(inout)  :: Ex_p(np), Ey_p(np)
    integer,          intent(in)     :: ng
    integer,          intent(in)     :: lo(2)
    integer,          intent(in)     :: hi(2)
    real(amrex_real), intent(in)     :: Ex(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng)
    real(amrex_real), intent(in)     :: Ey(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng)
    real(amrex_real), intent(in)     :: plo(2)
    real(amrex_real), intent(in)     :: dx(2)

    integer i, j, n
    real(amrex_real) wx_lo, wy_lo, wx_hi, wy_hi
    real(amrex_real) lx, ly
    real(amrex_real) inv_dx(2)
    inv_dx = 1.0d0/dx

    do n = 1, np
       lx = (particles(1, n) - plo(1))*inv_dx(1)
       ly = (particles(2, n) - plo(2))*inv_dx(2)

       i = floor(lx)
       j = floor(ly)

       wx_hi = lx - i
       wy_hi = ly - j

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi

       Ex_p(n) = wx_lo*wy_lo*Ex(i,   j  ) + &
                 wx_lo*wy_hi*Ex(i,   j+1) + &
                 wx_hi*wy_lo*Ex(i+1, j  ) + &
                 wx_hi*wy_hi*Ex(i+1, j+1)

       Ey_p(n) = wx_lo*wy_lo*Ey(i,   j  ) + &
                 wx_lo*wy_hi*Ey(i,   j+1) + &
                 wx_hi*wy_lo*Ey(i+1, j  ) + &
                 wx_hi*wy_hi*Ey(i+1, j+1)

    end do

  end subroutine warpx_interpolate_cic_2d


  subroutine warpx_interpolate_cic_two_levels_3d(particles, ns, np,      &
                                                 Ex_p, Ey_p, Ez_p,       &
                                                 Ex,   Ey,   Ez,         &
                                                 lo,   hi,   dx,         &
                                                 cEx,  cEy,  cEz,        &
                                                 mask,                   &
                                                 clo,  chi,  cdx,        &
                                                 plo,  ng,   lev)        &
       bind(c,name='warpx_interpolate_cic_two_levels_3d')
    integer, value,   intent(in)     :: ns, np
    real(amrex_particle_real), intent(in)     :: particles(ns,np)
    real(amrex_real), intent(inout)  :: Ex_p(np), Ey_p(np), Ez_p(np)
    integer,          intent(in)     :: ng, lev
    integer,          intent(in)     :: lo(3), hi(3)
    integer,          intent(in)     :: clo(3), chi(3)
    real(amrex_real), intent(in)     :: Ex(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng)
    real(amrex_real), intent(in)     :: Ey(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng)
    real(amrex_real), intent(in)     :: Ez(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng)
    real(amrex_real), intent(in)     :: cEx(clo(1)-ng:chi(1)+ng, clo(2)-ng:chi(2)+ng, clo(3)-ng:chi(3)+ng)
    real(amrex_real), intent(in)     :: cEy(clo(1)-ng:chi(1)+ng, clo(2)-ng:chi(2)+ng, clo(3)-ng:chi(3)+ng)
    real(amrex_real), intent(in)     :: cEz(clo(1)-ng:chi(1)+ng, clo(2)-ng:chi(2)+ng, clo(3)-ng:chi(3)+ng)
    integer(c_int),   intent(in)     :: mask (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(amrex_real), intent(in)     :: plo(3)
    real(amrex_real), intent(in)     :: dx(3), cdx(3)

    integer i, j, k, n
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) inv_dx(3), inv_cdx(3)
    inv_dx  = 1.0d0/dx
    inv_cdx = 1.0d0/cdx

    do n = 1, np

       lx = (particles(1, n) - plo(1))*inv_dx(1)
       ly = (particles(2, n) - plo(2))*inv_dx(2)
       lz = (particles(3, n) - plo(3))*inv_dx(3)

       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

! use the coarse E if near the level boundary
       if (lev .eq. 1 .and. mask(i,j,k) .eq. 1) then

          lx = (particles(1, n) - plo(1))*inv_cdx(1)
          ly = (particles(2, n) - plo(2))*inv_cdx(2)
          lz = (particles(3, n) - plo(3))*inv_cdx(3)

          i = floor(lx)
          j = floor(ly)
          k = floor(lz)

          wx_hi = lx - i
          wy_hi = ly - j
          wz_hi = lz - k

          wx_lo = 1.0d0 - wx_hi
          wy_lo = 1.0d0 - wy_hi
          wz_lo = 1.0d0 - wz_hi

          Ex_p(n) = wx_lo*wy_lo*wz_lo*cEx(i,   j,   k  ) + &
                    wx_lo*wy_lo*wz_hi*cEx(i,   j,   k+1) + &
                    wx_lo*wy_hi*wz_lo*cEx(i,   j+1, k  ) + &
                    wx_lo*wy_hi*wz_hi*cEx(i,   j+1, k+1) + &
                    wx_hi*wy_lo*wz_lo*cEx(i+1, j,   k  ) + &
                    wx_hi*wy_lo*wz_hi*cEx(i+1, j,   k+1) + &
                    wx_hi*wy_hi*wz_lo*cEx(i+1, j+1, k  ) + &
                    wx_hi*wy_hi*wz_hi*cEx(i+1, j+1, k+1)

          Ey_p(n) = wx_lo*wy_lo*wz_lo*cEy(i,   j,   k  ) + &
                    wx_lo*wy_lo*wz_hi*cEy(i,   j,   k+1) + &
                    wx_lo*wy_hi*wz_lo*cEy(i,   j+1, k  ) + &
                    wx_lo*wy_hi*wz_hi*cEy(i,   j+1, k+1) + &
                    wx_hi*wy_lo*wz_lo*cEy(i+1, j,   k  ) + &
                    wx_hi*wy_lo*wz_hi*cEy(i+1, j,   k+1) + &
                    wx_hi*wy_hi*wz_lo*cEy(i+1, j+1, k  ) + &
                    wx_hi*wy_hi*wz_hi*cEy(i+1, j+1, k+1)

          Ez_p(n) = wx_lo*wy_lo*wz_lo*cEz(i,   j,   k  ) + &
                    wx_lo*wy_lo*wz_hi*cEz(i,   j,   k+1) + &
                    wx_lo*wy_hi*wz_lo*cEz(i,   j+1, k  ) + &
                    wx_lo*wy_hi*wz_hi*cEz(i,   j+1, k+1) + &
                    wx_hi*wy_lo*wz_lo*cEz(i+1, j,   k  ) + &
                    wx_hi*wy_lo*wz_hi*cEz(i+1, j,   k+1) + &
                    wx_hi*wy_hi*wz_lo*cEz(i+1, j+1, k  ) + &
                    wx_hi*wy_hi*wz_hi*cEz(i+1, j+1, k+1)

! otherwise use the fine
       else

          wx_hi = lx - i
          wy_hi = ly - j
          wz_hi = lz - k

          wx_lo = 1.0d0 - wx_hi
          wy_lo = 1.0d0 - wy_hi
          wz_lo = 1.0d0 - wz_hi

          Ex_p(n) = wx_lo*wy_lo*wz_lo*Ex(i,   j,   k  ) + &
               wx_lo*wy_lo*wz_hi*Ex(i,   j,   k+1) + &
               wx_lo*wy_hi*wz_lo*Ex(i,   j+1, k  ) + &
               wx_lo*wy_hi*wz_hi*Ex(i,   j+1, k+1) + &
               wx_hi*wy_lo*wz_lo*Ex(i+1, j,   k  ) + &
               wx_hi*wy_lo*wz_hi*Ex(i+1, j,   k+1) + &
               wx_hi*wy_hi*wz_lo*Ex(i+1, j+1, k  ) + &
               wx_hi*wy_hi*wz_hi*Ex(i+1, j+1, k+1)

          Ey_p(n) = wx_lo*wy_lo*wz_lo*Ey(i,   j,   k  ) + &
                    wx_lo*wy_lo*wz_hi*Ey(i,   j,   k+1) + &
                    wx_lo*wy_hi*wz_lo*Ey(i,   j+1, k  ) + &
                    wx_lo*wy_hi*wz_hi*Ey(i,   j+1, k+1) + &
                    wx_hi*wy_lo*wz_lo*Ey(i+1, j,   k  ) + &
                    wx_hi*wy_lo*wz_hi*Ey(i+1, j,   k+1) + &
                    wx_hi*wy_hi*wz_lo*Ey(i+1, j+1, k  ) + &
                    wx_hi*wy_hi*wz_hi*Ey(i+1, j+1, k+1)

          Ez_p(n) = wx_lo*wy_lo*wz_lo*Ez(i,   j,   k  ) + &
                    wx_lo*wy_lo*wz_hi*Ez(i,   j,   k+1) + &
                    wx_lo*wy_hi*wz_lo*Ez(i,   j+1, k  ) + &
                    wx_lo*wy_hi*wz_hi*Ez(i,   j+1, k+1) + &
                    wx_hi*wy_lo*wz_lo*Ez(i+1, j,   k  ) + &
                    wx_hi*wy_lo*wz_hi*Ez(i+1, j,   k+1) + &
                    wx_hi*wy_hi*wz_lo*Ez(i+1, j+1, k  ) + &
                    wx_hi*wy_hi*wz_hi*Ez(i+1, j+1, k+1)

       end if

    end do

  end subroutine warpx_interpolate_cic_two_levels_3d


  subroutine warpx_interpolate_cic_two_levels_2d(particles, ns, np,      &
                                                 Ex_p, Ey_p,             &
                                                 Ex,   Ey,               &
                                                 lo,   hi,   dx,         &
                                                 cEx,  cEy,              &
                                                 mask,                   &
                                                 clo,  chi,  cdx,        &
                                                 plo,  ng,   lev)        &
       bind(c,name='warpx_interpolate_cic_two_levels_2d')
    integer, value,   intent(in)     :: ns, np
    real(amrex_particle_real), intent(in)     :: particles(ns,np)
    real(amrex_real), intent(inout)  :: Ex_p(np), Ey_p(np)
    integer,          intent(in)     :: ng, lev
    integer,          intent(in)     :: lo(2), hi(2)
    integer,          intent(in)     :: clo(2), chi(2)
    real(amrex_real), intent(in)     :: Ex(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng)
    real(amrex_real), intent(in)     :: Ey(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng)
    real(amrex_real), intent(in)     :: cEx(clo(1)-ng:chi(1)+ng, clo(2)-ng:chi(2)+ng)
    real(amrex_real), intent(in)     :: cEy(clo(1)-ng:chi(1)+ng, clo(2)-ng:chi(2)+ng)
    integer(c_int),   intent(in)     :: mask (lo(1):hi(1),lo(2):hi(2))
    real(amrex_real), intent(in)     :: plo(2)
    real(amrex_real), intent(in)     :: dx(2), cdx(2)

    integer i, j, n
    real(amrex_real) wx_lo, wy_lo, wx_hi, wy_hi
    real(amrex_real) lx, ly
    real(amrex_real) inv_dx(2), inv_cdx(2)
    inv_dx  = 1.0d0/dx
    inv_cdx = 1.0d0/cdx

    do n = 1, np

       lx = (particles(1, n) - plo(1))*inv_dx(1)
       ly = (particles(2, n) - plo(2))*inv_dx(2)

       i = floor(lx)
       j = floor(ly)

! use the coarse E if near the level boundary
       if (lev .eq. 1 .and. mask(i,j) .eq. 1) then

          lx = (particles(1, n) - plo(1))*inv_cdx(1)
          ly = (particles(2, n) - plo(2))*inv_cdx(2)

          i = floor(lx)
          j = floor(ly)

          wx_hi = lx - i
          wy_hi = ly - j

          wx_lo = 1.0d0 - wx_hi
          wy_lo = 1.0d0 - wy_hi

          Ex_p(n) = wx_lo*wy_lo*cEx(i,   j  ) + &
                    wx_lo*wy_hi*cEx(i,   j+1) + &
                    wx_hi*wy_lo*cEx(i+1, j  ) + &
                    wx_hi*wy_hi*cEx(i+1, j+1)

          Ey_p(n) = wx_lo*wy_lo*cEy(i,   j  ) + &
                    wx_lo*wy_hi*cEy(i,   j+1) + &
                    wx_hi*wy_lo*cEy(i+1, j  ) + &
                    wx_hi*wy_hi*cEy(i+1, j+1)

! otherwise use the fine
       else

          wx_hi = lx - i
          wy_hi = ly - j

          wx_lo = 1.0d0 - wx_hi
          wy_lo = 1.0d0 - wy_hi

          Ex_p(n) = wx_lo*wy_lo*Ex(i,   j  ) + &
                    wx_lo*wy_hi*Ex(i,   j+1) + &
                    wx_hi*wy_lo*Ex(i+1, j  ) + &
                    wx_hi*wy_hi*Ex(i+1, j+1)

          Ey_p(n) = wx_lo*wy_lo*Ey(i,   j  ) + &
                    wx_lo*wy_hi*Ey(i,   j+1) + &
                    wx_hi*wy_lo*Ey(i+1, j  ) + &
                    wx_hi*wy_hi*Ey(i+1, j+1)

       end if

    end do

  end subroutine warpx_interpolate_cic_two_levels_2d

end module warpx_ES_interpolate_cic
