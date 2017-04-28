
module warpx_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real

  implicit none

contains

  subroutine warpx_compute_E (lo, hi, &
       phi, phlo, phhi, &
       Ex,  Exlo, Exhi, &
       Ey,  Eylo, Eyhi, &
       Ez,  Ezlo, Ezhi, &
       dx) bind(c,name='warpx_compute_E')
    integer(c_int), intent(in) :: lo(3), hi(3), phlo(3), phhi(3), Exlo(3), Exhi(3),  &
         Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in   ) :: phi(phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3))
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3))
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3))

    integer :: i, j, k
    real(amrex_real) :: dxinv(3)

    dxinv = 1.0 / dx

    do k    = lo(3), hi(3)
       do j = lo(2), hi(2)

          do i = lo(1), hi(1)-1
             Ex(i,j,k) = dxinv(1) * (phi(i,j,k) - phi(i+1,j,k))
          end do

          if (j < hi(2)) then
             do i = lo(1), hi(1)
                Ey(i,j,k) = dxinv(2) * (phi(i,j,k) - phi(i,j+1,k))
             end do
          end if

          if (k < hi(3)) then
             do i = lo(1), hi(1)
                Ez(i,j,k) = dxinv(3) * (phi(i,j,k) - phi(i,j,k+1))
             end do
          end if

       end do
    end do

  end subroutine warpx_compute_E

  subroutine warpx_deposit_cic(particles, ns, np,                     &
                               weights, charge, rho, lo, hi, plo, dx) &
       bind(c,name='warpx_deposit_cic')
    integer, value       :: ns, np
    real(amrex_real)     :: particles(ns,np)
    real(amrex_real)     :: weights(np)
    real(amrex_real)     :: charge
    integer              :: lo(3)
    integer              :: hi(3)
    real(amrex_real)     :: rho(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
    real(amrex_real)     :: plo(3)
    real(amrex_real)     :: dx(3)

    integer i, j, k, n
    real(amrex_real) wx_lo, wy_lo, wz_lo, wx_hi, wy_hi, wz_hi
    real(amrex_real) lx, ly, lz
    real(amrex_real) inv_dx(3)
    real(amrex_real) qp

    inv_dx = 1.0d0/dx

    do n = 1, np

       qp = weights(n) * charge

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

  end subroutine warpx_deposit_cic

  subroutine warpx_interpolate_cic(particles, ns, np,           &
                                   acc, lo, hi, ncomp, plo, dx) &
       bind(c,name='warpx_interpolate_cic')
    integer, value       :: ns, np, ncomp
    real(amrex_real)     :: particles(ns,np)
    integer              :: lo(3)
    integer              :: hi(3)
    real(amrex_real)     :: acc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), ncomp)
    real(amrex_real)     :: plo(3)
    real(amrex_real)     :: dx(3)
    real(amrex_real)     :: acceleration(ncomp)

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
       
          if (abs(acceleration(nc) - 5.d0) .ge. 1.0d-9) then
             print *, particles(1, n), particles(2, n), particles(3, n)
          end if

       end do
    end do

  end subroutine warpx_interpolate_cic

end module warpx_module
