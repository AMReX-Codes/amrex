module warpx_ES_deposit_cic

  use iso_c_binding
  use amrex_fort_module, only : amrex_real, amrex_particle_real

  implicit none

contains

! This routine computes the charge density due to the particles using cloud-in-cell
! deposition. The Fab rho is assumed to be node-centered.
!
! Arguments:
!     particles : a pointer to the particle array-of-structs
!     ns        : the stride length of particle struct (the size of the struct in number of reals)
!     np        : the number of particles
!     weights   : the particle weights
!     charge    : the charge of this particle species
!     rho       : a Fab that will contain the charge density on exit
!     lo        : a pointer to the lo corner of this valid box for rho, in index space
!     hi        : a pointer to the hi corner of this valid box for rho, in index space
!     plo       : the real position of the left-hand corner of the problem domain
!     dx        : the mesh spacing
!     ng        : the number of ghost cells in rho
!
  subroutine warpx_deposit_cic_3d(particles, ns, np,                     &
                                  weights, charge, rho, lo, hi, plo, dx, &
                                  ng)                                    &
       bind(c,name='warpx_deposit_cic_3d')
    integer, value,   intent(in)     :: ns, np
    real(amrex_particle_real), intent(in)     :: particles(ns,np)
    real(amrex_particle_real), intent(in)     :: weights(np)
    real(amrex_real), intent(in)     :: charge
    integer,          intent(in)     :: lo(3)
    integer,          intent(in)     :: hi(3)
    integer,          intent(in)     :: ng
    real(amrex_real), intent(inout)  :: rho(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng)
    real(amrex_real), intent(in)     :: plo(3)
    real(amrex_real), intent(in)     :: dx(3)

    integer i, j, k, n
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) inv_dx(3)
    real(amrex_real) qp, inv_vol

    inv_dx = 1.0d0/dx

    inv_vol = inv_dx(1) * inv_dx(2) * inv_dx(3)

    do n = 1, np

       qp = weights(n) * charge * inv_vol

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

       rho(i,   j,   k  ) = rho(i,   j,   k  ) + wx_lo*wy_lo*wz_lo*qp
       rho(i,   j,   k+1) = rho(i,   j,   k+1) + wx_lo*wy_lo*wz_hi*qp
       rho(i,   j+1, k  ) = rho(i,   j+1, k  ) + wx_lo*wy_hi*wz_lo*qp
       rho(i,   j+1, k+1) = rho(i,   j+1, k+1) + wx_lo*wy_hi*wz_hi*qp
       rho(i+1, j,   k  ) = rho(i+1, j,   k  ) + wx_hi*wy_lo*wz_lo*qp
       rho(i+1, j,   k+1) = rho(i+1, j,   k+1) + wx_hi*wy_lo*wz_hi*qp
       rho(i+1, j+1, k  ) = rho(i+1, j+1, k  ) + wx_hi*wy_hi*wz_lo*qp
       rho(i+1, j+1, k+1) = rho(i+1, j+1, k+1) + wx_hi*wy_hi*wz_hi*qp

    end do

  end subroutine warpx_deposit_cic_3d

  subroutine warpx_deposit_cic_2d(particles, ns, np,                     &
                                  weights, charge, rho, lo, hi, plo, dx, &
                                  ng)                                    &
       bind(c,name='warpx_deposit_cic_2d')
    integer, value,   intent(in)     :: ns, np
    real(amrex_particle_real), intent(in)     :: particles(ns,np)
    real(amrex_particle_real), intent(in)     :: weights(np)
    real(amrex_real), intent(in)     :: charge
    integer,          intent(in)     :: lo(2)
    integer,          intent(in)     :: hi(2)
    integer,          intent(in)     :: ng
    real(amrex_real), intent(inout)  :: rho(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng)
    real(amrex_real), intent(in)     :: plo(2)
    real(amrex_real), intent(in)     :: dx(2)

    integer i, j, n
    real(amrex_real) wx_lo, wy_lo, wx_hi, wy_hi
    real(amrex_real) lx, ly
    real(amrex_real) inv_dx(2)
    real(amrex_real) qp, inv_vol

    inv_dx = 1.0d0/dx

    inv_vol = inv_dx(1) * inv_dx(2)

    do n = 1, np

       qp = weights(n) * charge * inv_vol

       lx = (particles(1, n) - plo(1))*inv_dx(1)
       ly = (particles(2, n) - plo(2))*inv_dx(2)

       i = floor(lx)
       j = floor(ly)

       wx_hi = lx - i
       wy_hi = ly - j

       wx_lo = 1.0d0 - wx_hi
       wy_lo = 1.0d0 - wy_hi

       rho(i,   j  ) = rho(i,   j  ) + wx_lo*wy_lo*qp
       rho(i,   j+1) = rho(i,   j+1) + wx_lo*wy_hi*qp
       rho(i+1, j  ) = rho(i+1, j  ) + wx_hi*wy_lo*qp
       rho(i+1, j+1) = rho(i+1, j+1) + wx_hi*wy_hi*qp

    end do

  end subroutine warpx_deposit_cic_2d

end module warpx_ES_deposit_cic
