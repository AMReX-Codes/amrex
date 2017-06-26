
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
    real(amrex_real), intent(in)  :: dx(3)
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

! This routine computes the node-centered electric field given a node-centered phi.
! The gradient is computed using 2nd-order centered differences. It assumes the 
! Boundary conditions have already been set and that you have one row of ghost cells.
! Note that this routine includes the minus sign in E = - grad phi.
! This routine is used when running WarpX in electrostatic mode.
!
! Arguments:
!     lo, hi:     The corners of the valid box over which the gradient is taken
!     Ex, Ey, Ez: The electric field in the x, y, and z directions.
!     dx:         The cell spacing
!
  subroutine warpx_compute_E_nodal (lo, hi, phi, Ex, Ey, Ez, dx) &
       bind(c,name='warpx_compute_E_nodal')
    integer(c_int),   intent(in)    :: lo(3), hi(3)
    real(amrex_real), intent(in)    :: dx(3)
    real(amrex_real), intent(in   ) :: phi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
    real(amrex_real), intent(inout) :: Ex (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
    real(amrex_real), intent(inout) :: Ey (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
    real(amrex_real), intent(inout) :: Ez (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

    integer :: i, j, k
    real(amrex_real) :: fac(3)

    fac = 0.5d0 / dx

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             Ex(i,j,k) = fac(1) * (phi(i-1,j,k) - phi(i+1,j,k))
             Ey(i,j,k) = fac(2) * (phi(i,j-1,k) - phi(i,j+1,k))
             Ez(i,j,k) = fac(3) * (phi(i,j,k-1) - phi(i,j,k+1))

          end do
       end do
    end do

  end subroutine warpx_compute_E_nodal

  subroutine warpx_deposit_cic(particles, ns, np,                     &
                               weights, charge, rho, lo, hi, plo, dx) &
       bind(c,name='warpx_deposit_cic')
    integer, value       :: ns, np
    real(amrex_real)     :: particles(ns,np)
    real(amrex_real)     :: weights(np)
    real(amrex_real)     :: charge
    integer              :: lo(3)
    integer              :: hi(3)
    real(amrex_real)     :: rho(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1)
    real(amrex_real)     :: plo(3)
    real(amrex_real)     :: dx(3)

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

  end subroutine warpx_deposit_cic

  subroutine warpx_interpolate_cic(particles, ns, np,      &
                                   Ex_p, Ey_p, Ez_p,       &
                                   Ex,   Ey,   Ez,         &
                                   lo, hi, plo, dx)        &
       bind(c,name='warpx_interpolate_cic')
    integer, value       :: ns, np
    real(amrex_real)     :: particles(ns,np)
    real(amrex_real)     :: Ex_p(np), Ey_p(np), Ez_p(np)
    integer              :: lo(3)
    integer              :: hi(3)
    real(amrex_real)     :: Ex(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1)
    real(amrex_real)     :: Ey(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1)
    real(amrex_real)     :: Ez(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1)
    real(amrex_real)     :: plo(3)
    real(amrex_real)     :: dx(3)

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

  end subroutine warpx_interpolate_cic

  subroutine warpx_push_leapfrog(particles, ns, np,      &
                                 vx_p, vy_p, vz_p,       &                                 
                                 Ex_p, Ey_p, Ez_p,       &
                                 charge, mass, dt,       &
                                 prob_lo, prob_hi)       &
       bind(c,name='warpx_push_leapfrog')
    integer, value       :: ns, np
    real(amrex_real)     :: particles(ns,np)
    real(amrex_real)     :: vx_p(np), vy_p(np), vz_p(np)
    real(amrex_real)     :: Ex_p(np), Ey_p(np), Ez_p(np)
    real(amrex_real)     :: charge
    real(amrex_real)     :: mass
    real(amrex_real)     :: dt
    real(amrex_real)     :: prob_lo(3), prob_hi(3)
   
    integer n
    real(amrex_real) fac

    fac = charge * dt / mass

    do n = 1, np

       vx_p(n) = vx_p(n) + fac * Ex_p(n)
       vy_p(n) = vy_p(n) + fac * Ey_p(n)
       vz_p(n) = vz_p(n) + fac * Ez_p(n)

       particles(1, n) = particles(1, n) + dt * vx_p(n)
       particles(2, n) = particles(2, n) + dt * vy_p(n)
       particles(3, n) = particles(3, n) + dt * vz_p(n)
              
!      bounce off the walls in the x...
       do while (particles(1, n) .lt. prob_lo(1) .or. particles(1, n) .gt. prob_hi(1))
          if (particles(1, n) .lt. prob_lo(1)) then
             particles(1, n) = 2.d0*prob_lo(1) - particles(1, n)
          else
             particles(1, n) = 2.d0*prob_hi(1) - particles(1, n)
          end if
          vx_p(n) = -vx_p(n)
       end do

!      ... y... 
       do while (particles(2, n) .lt. prob_lo(2) .or. particles(2, n) .gt. prob_hi(2))
          if (particles(2, n) .lt. prob_lo(2)) then
             particles(2, n) = 2.d0*prob_lo(2) - particles(2, n)
          else
             particles(2, n) = 2.d0*prob_hi(2) - particles(2, n)
          end if
          vy_p(n) = -vy_p(n)
       end do

!      ... and z directions
       do while (particles(3, n) .lt. prob_lo(3) .or. particles(3, n) .gt. prob_hi(3))
          if (particles(3, n) .lt. prob_lo(3)) then
             particles(3, n) = 2.d0*prob_lo(3) - particles(3, n)
          else
             particles(3, n) = 2.d0*prob_hi(3) - particles(3, n)
          end if
          vz_p(n) = -vz_p(n)
       end do

    end do

  end subroutine warpx_push_leapfrog

  subroutine warpx_push_leapfrog_positions(particles, ns, np,      &
                                           vx_p, vy_p, vz_p, dt)   &
       bind(c,name='warpx_push_leapfrog_positions')
    integer, value       :: ns, np
    real(amrex_real)     :: particles(ns,np)
    real(amrex_real)     :: vx_p(np), vy_p(np), vz_p(np)
    real(amrex_real)     :: dt
    
    integer n

    do n = 1, np

       particles(1, n) = particles(1, n) + dt * vx_p(n)
       particles(2, n) = particles(2, n) + dt * vy_p(n)
       particles(3, n) = particles(3, n) + dt * vz_p(n)
              
    end do

  end subroutine warpx_push_leapfrog_positions

  subroutine warpx_compute_divb_3d (lo, hi, divB, dlo, dhi, &
       Bx, xlo, xhi, By, ylo, yhi, Bz, zlo, zhi, dx) &
       bind(c, name='warpx_compute_divb_3d')
    integer, intent(in) :: lo(3),hi(3),dlo(3),dhi(3),xlo(3),xhi(3),ylo(3),yhi(3),zlo(3),zhi(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in   ) :: Bx  (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(amrex_real), intent(in   ) :: By  (ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    real(amrex_real), intent(in   ) :: Bz  (zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3))
    real(amrex_real), intent(inout) :: divB(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

    integer :: i,j,k
    real(amrex_real) :: dxinv(3)

    dxinv = 1.d0/dx

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             divB(i,j,k) = dxinv(1) * (Bx(i+1,j  ,k  ) - Bx(i,j,k)) &
                  +        dxinv(2) * (By(i  ,j+1,k  ) - By(i,j,k)) &
                  +        dxinv(3) * (Bz(i  ,j  ,k+1) - Bz(i,j,k))
          end do
       end do
    end do
  end subroutine warpx_compute_divb_3d


  subroutine warpx_compute_divb_2d (lo, hi, divB, dlo, dhi, &
       Bx, xlo, xhi, By, ylo, yhi, Bz, zlo, zhi, dx) &
       bind(c, name='warpx_compute_divb_2d')
    integer, intent(in) :: lo(2),hi(2),dlo(2),dhi(2),xlo(2),xhi(2),ylo(2),yhi(2),zlo(2),zhi(2)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in   ) :: Bx  (xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in   ) :: By  (ylo(1):yhi(1),ylo(2):yhi(2))
    real(amrex_real), intent(in   ) :: Bz  (zlo(1):zhi(1),zlo(2):zhi(2))
    real(amrex_real), intent(inout) :: divB(dlo(1):dhi(1),dlo(2):dhi(2))

    integer :: i,k
    real(amrex_real) :: dxinv(3)

    dxinv = 1.d0/dx

    do    k = lo(2), hi(2)
       do i = lo(1), hi(1)
          divB(i,k) = dxinv(1) * (Bx(i+1,k  ) - Bx(i,k)) &
               +      dxinv(3) * (Bz(i  ,k+1) - Bz(i,k))
       end do
    end do
  end subroutine warpx_compute_divb_2d


  subroutine warpx_compute_dive_3d (lo, hi, dive, dlo, dhi, &
       Ex, xlo, xhi, Ey, ylo, yhi, Ez, zlo, zhi, dx) &
       bind(c, name='warpx_compute_dive_3d')
    integer, intent(in) :: lo(3),hi(3),dlo(3),dhi(3),xlo(3),xhi(3),ylo(3),yhi(3),zlo(3),zhi(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in   ) :: Ex  (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(amrex_real), intent(in   ) :: Ey  (ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    real(amrex_real), intent(in   ) :: Ez  (zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3))
    real(amrex_real), intent(inout) :: dive(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

    integer :: i,j,k
    real(amrex_real) :: dxinv(3)

    dxinv = 1.d0/dx

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dive(i,j,k) = dxinv(1) * (Ex(i,j,k) - Ex(i-1,j,k)) &
                  +        dxinv(2) * (Ey(i,j,k) - Ey(i,j-1,k)) &
                  +        dxinv(3) * (Ez(i,j,k) - Ez(i,j,k-1))
          end do
       end do
    end do
  end subroutine warpx_compute_dive_3d


  subroutine warpx_compute_dive_2d (lo, hi, dive, dlo, dhi, &
       Ex, xlo, xhi, Ey, ylo, yhi, Ez, zlo, zhi, dx) &
       bind(c, name='warpx_compute_dive_2d')
    integer, intent(in) :: lo(2),hi(2),dlo(2),dhi(2),xlo(2),xhi(2),ylo(2),yhi(2),zlo(2),zhi(2)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in   ) :: Ex  (xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in   ) :: Ey  (ylo(1):yhi(1),ylo(2):yhi(2))
    real(amrex_real), intent(in   ) :: Ez  (zlo(1):zhi(1),zlo(2):zhi(2))
    real(amrex_real), intent(inout) :: dive(dlo(1):dhi(1),dlo(2):dhi(2))

    integer :: i,k
    real(amrex_real) :: dxinv(3)

    dxinv = 1.d0/dx

    do    k = lo(2), hi(2)
       do i = lo(1), hi(1)
          dive(i,k) = dxinv(1) * (Ex(i,k) - Ex(i-1,k)) &
               +      dxinv(3) * (Ez(i,k) - Ez(i,k-1))
       end do
    end do
  end subroutine warpx_compute_dive_2d


  subroutine warpx_sync_current_2d (lo, hi, crse, clo, chi, fine, flo, fhi, dir) &
       bind(c, name='warpx_sync_current_2d')
    integer, intent(in) :: lo(2), hi(2), flo(2), fhi(2), clo(2), chi(2), dir
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2))
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2))

    integer :: i,j,ii,jj

    if (dir == 0) then
       do j = lo(2), hi(2)
          jj = j*2
          do i = lo(1), hi(1)
             ii = i*2
             crse(i,j) = 0.25d0 * (fine(ii,jj) + fine(ii+1,jj) &
                  + 0.5d0*(fine(ii,jj-1) + fine(ii+1,jj-1) + fine(ii,jj+1) + fine(ii+1,jj+1)) )
          end do
       end do
    else if (dir == 2) then
       do j = lo(2), hi(2)
          jj = j*2
          do i = lo(1), hi(1)
             ii = i*2
             crse(i,j) = 0.25d0 * (fine(ii,jj) + fine(ii,jj+1) &
                  + 0.5d0*(fine(ii-1,jj) + fine(ii-1,jj+1) + fine(ii+1,jj) + fine(ii+1,jj+1)) )
          end do
       end do
    else
       do j = lo(2), hi(2)
          jj = j*2
          do i = lo(1), hi(1)
             ii = i*2
             crse(i,j) = 0.25d0 * &
                     ( fine(ii,jj) + 0.5d0 *(fine(ii-1,jj  )+fine(ii+1,jj  ) &
                     &                     + fine(ii  ,jj-1)+fine(ii  ,jj+1)) &
                     &             + 0.25d0*(fine(ii-1,jj-1)+fine(ii+1,jj-1) &
                     &                     + fine(ii-1,jj+1)+fine(ii+1,jj+1)) ) 
          end do
       end do
    end if
  end subroutine warpx_sync_current_2d

  subroutine warpx_sync_current_3d (lo, hi, crse, clo, chi, fine, flo, fhi, dir) &
       bind(c, name='warpx_sync_current_3d')
    integer, intent(in) :: lo(3), hi(3), flo(3), fhi(3), clo(3), chi(3), dir
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))

    integer :: i,j,k, ii,jj,kk

    if (dir == 0) then
       do k = lo(3), hi(3)
          kk = k*2
          do j = lo(2), hi(2)
             jj = j*2
             do i = lo(1), hi(1)
                ii = i*2
                crse(i,j,k) = 0.125d0* &
                     ( fine(ii  ,jj,kk) + 0.5d0 *(fine(ii  ,jj-1,kk  )+fine(ii  ,jj+1,kk  ) &
                     &                          + fine(ii  ,jj  ,kk-1)+fine(ii  ,jj  ,kk+1)) &
                     &                  + 0.25d0*(fine(ii  ,jj-1,kk-1)+fine(ii  ,jj+1,kk-1) &
                     &                          + fine(ii  ,jj-1,kk+1)+fine(ii  ,jj+1,kk+1)) &
                     + fine(ii+1,jj,kk) + 0.5d0 *(fine(ii+1,jj-1,kk  )+fine(ii+1,jj+1,kk  ) &
                     &                          + fine(ii+1,jj  ,kk-1)+fine(ii+1,jj  ,kk+1)) &
                     &                  + 0.25d0*(fine(ii+1,jj-1,kk-1)+fine(ii+1,jj+1,kk-1) &
                     &                          + fine(ii+1,jj-1,kk+1)+fine(ii+1,jj+1,kk+1)) )
             end do
          end do
       end do
    else if (dir == 1) then
       do k = lo(3), hi(3)
          kk = k*2
          do j = lo(2), hi(2)
             jj = j*2
             do i = lo(1), hi(1)
                ii = i*2
                crse(i,j,k) = 0.125d0* &
                     ( fine(ii,jj  ,kk) + 0.5d0 *(fine(ii-1,jj  ,kk  )+fine(ii+1,jj  ,kk  ) &
                     &                          + fine(ii  ,jj  ,kk-1)+fine(ii  ,jj  ,kk+1)) &
                     &                  + 0.25d0*(fine(ii-1,jj  ,kk-1)+fine(ii+1,jj  ,kk-1) &
                     &                          + fine(ii-1,jj  ,kk+1)+fine(ii+1,jj  ,kk+1)) &
                     + fine(ii,jj+1,kk) + 0.5d0 *(fine(ii-1,jj+1,kk  )+fine(ii+1,jj+1,kk  ) &
                     &                          + fine(ii  ,jj+1,kk-1)+fine(ii  ,jj+1,kk+1)) &
                     &                  + 0.25d0*(fine(ii-1,jj+1,kk-1)+fine(ii+1,jj+1,kk-1) &
                     &                          + fine(ii-1,jj+1,kk+1)+fine(ii+1,jj+1,kk+1)) )
             end do
          end do
       end do
    else
       do k = lo(3), hi(3)
          kk = k*2
          do j = lo(2), hi(2)
             jj = j*2
             do i = lo(1), hi(1)
                ii = i*2
                crse(i,j,k) = 0.125d0* &
                     ( fine(ii,jj,kk  ) + 0.5d0 *(fine(ii-1,jj  ,kk  )+fine(ii+1,jj  ,kk  ) &
                     &                          + fine(ii  ,jj-1,kk  )+fine(ii  ,jj+1,kk  )) &
                     &                  + 0.25d0*(fine(ii-1,jj-1,kk  )+fine(ii+1,jj-1,kk  ) &
                     &                          + fine(ii-1,jj+1,kk  )+fine(ii+1,jj+1,kk  )) &
                     + fine(ii,jj,kk+1) + 0.5d0 *(fine(ii-1,jj  ,kk+1)+fine(ii+1,jj  ,kk+1) &
                     &                          + fine(ii  ,jj-1,kk+1)+fine(ii  ,jj+1,kk+1)) &
                     &                  + 0.25d0*(fine(ii-1,jj-1,kk+1)+fine(ii+1,jj-1,kk+1) &
                     &                          + fine(ii-1,jj+1,kk+1)+fine(ii+1,jj+1,kk+1)) )
             end do
          end do
       end do
    end if
  end subroutine warpx_sync_current_3d

  subroutine warpx_clean_evec_2d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                          Ex, Exlo, Exhi, &
       &                          Ey, Eylo, Eyhi, &
       &                          Ez, Ezlo, Ezhi, &
       &                          F,  flo,  fhi,  &
       dtdx) bind(c, name='warpx_clean_evec_2d')
    integer, intent(in) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
         Exlo(2), Exhi(2), Eylo(2), Eyhi(2), Ezlo(2), Ezhi(2), flo(2), fhi(2)
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2))
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2))
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2))
    real(amrex_real), intent(in   ) :: F  ( flo(1): fhi(1), flo(2): fhi(2))
    real(amrex_real), intent(in   ) :: dtdx(3)

    integer :: i,j
    
    do    j = xlo(2), xhi(2)
       do i = xlo(1), xhi(1)
          Ex(i,j) = Ex(i,j) + dtdx(1) * (F(i+1,j)-F(i,j))
       end do
    end do

    ! Ey doesn't change

    do    j = zlo(2), zhi(2)
       do i = zlo(1), zhi(1)
          Ez(i,j) = Ez(i,j) + dtdx(3) * (F(i,j+1)-F(i,j))
       end do
    end do

  end subroutine warpx_clean_evec_2d

  subroutine warpx_clean_evec_3d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                          Ex, Exlo, Exhi, &
       &                          Ey, Eylo, Eyhi, &
       &                          Ez, Ezlo, Ezhi, &
       &                          F,  flo,  fhi,  &
       dtdx) bind(c, name='warpx_clean_evec_3d')
    integer, intent(in) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
         Exlo(3), Exhi(3), Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3), flo(3), fhi(3)
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3))
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3))
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3))
    real(amrex_real), intent(in   ) :: F  ( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3))
    real(amrex_real), intent(in   ) :: dtdx(3)

    integer :: i,j,k

    do       k = xlo(3), xhi(3)
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             Ex(i,j,k) = Ex(i,j,k) + dtdx(1)*(F(i+1,j,k)-F(i,j,k))
          end do
       end do
    end do

    do       k = xlo(3), xhi(3)
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             Ey(i,j,k) = Ey(i,j,k) + dtdx(2)*(F(i,j+1,k)-F(i,j,k))
          end do
       end do
    end do


    do       k = xlo(3), xhi(3)
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             Ez(i,j,k) = Ez(i,j,k) + dtdx(3)*(F(i,j,k+1)-F(i,j,k))
          end do
       end do
    end do

  end subroutine warpx_clean_evec_3d


  subroutine warpx_sync_rho_2d (lo, hi, crse, clo, chi, fine, flo, fhi) &
       bind(c, name='warpx_sync_rho_2d')
    integer, intent(in) :: lo(2), hi(2), flo(2), fhi(2), clo(2), chi(2)
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2))
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2))
    
    integer :: i,j,ii,jj

    do j = lo(2), hi(2)
       jj = j*2
       do i = lo(1), hi(1)
          ii = i*2
          crse(i,j) = 0.25d0 * &
               ( fine(ii,jj) + 0.5d0 *(fine(ii-1,jj  )+fine(ii+1,jj  ) &
               &                     + fine(ii  ,jj-1)+fine(ii  ,jj+1)) &
               &             + 0.25d0*(fine(ii-1,jj-1)+fine(ii+1,jj-1) &
               &                     + fine(ii-1,jj+1)+fine(ii+1,jj+1)) ) 
       end do
    end do
  end subroutine warpx_sync_rho_2d

  subroutine warpx_sync_rho_3d (lo, hi, crse, clo, chi, fine, flo, fhi) &
       bind(c, name='warpx_sync_rho_3d')
    integer, intent(in) :: lo(3), hi(3), flo(3), fhi(3), clo(3), chi(3)
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))

    integer :: i,j,k,ii,jj,kk

    do k = lo(3), hi(3)
       kk = k*2
       do j = lo(2), hi(2)
          jj = j*2
          do i = lo(1), hi(1)
             ii = i*2
             crse(i,j,k) = 0.125d0 * &
                  (fine(ii,jj,kk) + 0.5d0  *(fine(ii-1,jj  ,kk  )+fine(ii+1,jj  ,kk  ) &
                  &                        + fine(ii  ,jj-1,kk  )+fine(ii  ,jj+1,kk  ) &
                  &                        + fine(ii  ,jj  ,kk-1)+fine(ii  ,jj  ,kk+1)) &
                  &               + 0.25d0 *(fine(ii-1,jj-1,kk  )+fine(ii+1,jj-1,kk  ) &
                  &                        + fine(ii-1,jj+1,kk  )+fine(ii+1,jj+1,kk  ) &
                  &                        + fine(ii-1,jj  ,kk-1)+fine(ii+1,jj  ,kk-1) &
                  &                        + fine(ii-1,jj  ,kk+1)+fine(ii+1,jj  ,kk+1) &
                  &                        + fine(ii  ,jj-1,kk-1)+fine(ii  ,jj+1,kk-1) &
                  &                        + fine(ii  ,jj-1,kk+1)+fine(ii  ,jj+1,kk+1)) &
                  &               + 0.125d0*(fine(ii-1,jj-1,kk-1)+fine(ii-1,jj-1,kk+1) &
                  &                        + fine(ii-1,jj+1,kk-1)+fine(ii-1,jj+1,kk+1) &
                  &                        + fine(ii+1,jj-1,kk-1)+fine(ii+1,jj-1,kk+1) &
                  &                        + fine(ii+1,jj+1,kk-1)+fine(ii+1,jj+1,kk+1)))
          end do
       end do
    end do

  end subroutine warpx_sync_rho_3d

end module warpx_module
