
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

  subroutine warpx_push_pml_bvec_3d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                             Ex, Exlo, Exhi, &
       &                             Ey, Eylo, Eyhi, &
       &                             Ez, Ezlo, Ezhi, &
       &                             Bx, Bxlo, Bxhi, &
       &                             By, Bylo, Byhi, &
       &                             Bz, Bzlo, Bzhi, &
       &                             sigx1, sigx1_lo, sigx1_hi, &
       &                             sigx2, sigx2_lo, sigx2_hi, &
       &                             sigy1, sigy1_lo, sigy1_hi, &
       &                             sigy2, sigy2_lo, sigy2_hi, &
       &                             sigz1, sigz1_lo, sigz1_hi, &
       &                             sigz2, sigz2_lo, sigz2_hi) &
       bind(c,name='warpx_push_pml_bvec_3d')
    integer, intent(in) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
         Exlo(3), Exhi(3), Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3), &
         Bxlo(3), Bxhi(3), Bylo(3), Byhi(3), Bzlo(3), Bzhi(3)
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigy1_lo, sigy1_hi, sigy2_lo, sigy2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(in   ) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3),2)
    real(amrex_real), intent(in   ) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3),2)
    real(amrex_real), intent(in   ) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3),2)
    real(amrex_real), intent(inout) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),Bxlo(3):Bxhi(3),2)
    real(amrex_real), intent(inout) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),Bylo(3):Byhi(3),2)
    real(amrex_real), intent(inout) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),Bzlo(3):Bzhi(3),2)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigy1(sigy1_lo:sigy1_hi)
    real(amrex_real), intent(in) :: sigy2(sigy2_lo:sigy2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)

    integer :: i, j, k

    do       k = xlo(3), xhi(3)
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             Bx(i,j,k,1) = sigy1(j)*Bx(i,j,k,1) - sigy2(j)*(Ez(i,j+1,k  ,1)+Ez(i,j+1,k  ,2) &
                  &                                        -Ez(i,j  ,k  ,1)-Ez(i,j  ,k  ,2))
             Bx(i,j,k,2) = sigz1(k)*Bx(i,j,k,2) + sigz2(k)*(Ey(i,j  ,k+1,1)+Ey(i,j  ,k+1,2) &
                  &                                        -Ey(i,j  ,k  ,1)-Ey(i,j  ,k  ,2))
          end do
       end do
    end do

    do       k = ylo(3), yhi(3)
       do    j = ylo(2), yhi(2)
          do i = ylo(1), yhi(1)
             By(i,j,k,1) = sigz1(k)*By(i,j,k,1) - sigz2(k)*(Ex(i  ,j,k+1,1)+Ex(i  ,j,k+1,2) &
                  &                                        -Ex(i  ,j,k  ,1)-Ex(i  ,j,k  ,2))
             By(i,j,k,2) = sigx1(i)*By(i,j,k,2) + sigx2(i)*(Ez(i+1,j,k  ,1)+Ez(i+1,j,k  ,2) &
                  &                                        -Ez(i  ,j,k  ,1)-Ez(i  ,j,k  ,2))
          end do
       end do
    end do

    do       k = zlo(3), zhi(3)
       do    j = zlo(2), zhi(2)
          do i = zlo(1), zhi(1)
             Bz(i,j,k,1) = sigx1(i)*Bz(i,j,k,1) - sigx2(i)*(Ey(i+1,j  ,k,1)+Ey(i+1,j  ,k,2) &
                  &                                        -Ey(i  ,j  ,k,1)-Ey(i  ,j  ,k,2))
             Bz(i,j,k,2) = sigy1(j)*Bz(i,j,k,2) + sigy2(j)*(Ex(i  ,j+1,k,1)+Ex(i  ,j+1,k,2) &
                  &                                        -Ex(i  ,j  ,k,1)-Ex(i  ,j  ,k,2))
          end do
       end do
    end do

  end subroutine warpx_push_pml_bvec_3d

  subroutine warpx_push_pml_evec_3d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                             Ex, Exlo, Exhi, &
       &                             Ey, Eylo, Eyhi, &
       &                             Ez, Ezlo, Ezhi, &
       &                             Bx, Bxlo, Bxhi, &
       &                             By, Bylo, Byhi, &
       &                             Bz, Bzlo, Bzhi, &
       &                             sigx1, sigx1_lo, sigx1_hi, &
       &                             sigx2, sigx2_lo, sigx2_hi, &
       &                             sigy1, sigy1_lo, sigy1_hi, &
       &                             sigy2, sigy2_lo, sigy2_hi, &
       &                             sigz1, sigz1_lo, sigz1_hi, &
       &                             sigz2, sigz2_lo, sigz2_hi) &
       bind(c,name='warpx_push_pml_evec_3d')
    integer, intent(in) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
         Exlo(3), Exhi(3), Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3), &
         Bxlo(3), Bxhi(3), Bylo(3), Byhi(3), Bzlo(3), Bzhi(3)
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigy1_lo, sigy1_hi, sigy2_lo, sigy2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3),2)
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3),2)
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3),2)
    real(amrex_real), intent(in   ) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),Bxlo(3):Bxhi(3),2)
    real(amrex_real), intent(in   ) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),Bylo(3):Byhi(3),2)
    real(amrex_real), intent(in   ) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),Bzlo(3):Bzhi(3),2)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigy1(sigy1_lo:sigy1_hi)
    real(amrex_real), intent(in) :: sigy2(sigy2_lo:sigy2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)

    integer :: i, j, k

    do       k = xlo(3), xhi(3)
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             Ex(i,j,k,1) = sigy1(j)*Ex(i,j,k,1) + sigy2(j)*(Bz(i,j  ,k  ,1)+Bz(i,j  ,k  ,2) &
                  &                                        -Bz(i,j-1,k  ,1)-Bz(i,j-1,k  ,2))
             Ex(i,j,k,2) = sigz1(k)*Ex(i,j,k,2) - sigz2(k)*(By(i,j  ,k  ,1)+By(i,j  ,k  ,2) &
                  &                                        -By(i,j  ,k-1,1)-By(i,j  ,k-1,2))
          end do
       end do
    end do

    do       k = ylo(3), yhi(3)
       do    j = ylo(2), yhi(2)
          do i = ylo(1), yhi(1)
             Ey(i,j,k,1) = sigz1(k)*Ey(i,j,k,1) + sigz2(k)*(Bx(i  ,j,k  ,1)+Bx(i  ,j,k  ,2) &
                  &                                        -Bx(i  ,j,k-1,1)-Bx(i  ,j,k-1,2))
             Ey(i,j,k,2) = sigx1(i)*Ey(i,j,k,2) - sigx2(i)*(Bz(i  ,j,k  ,1)+Bz(i  ,j,k  ,2) &
                  &                                        -Bz(i-1,j,k  ,1)-Bz(i-1,j,k  ,2))
          end do
       end do
    end do

    do       k = zlo(3), zhi(3)
       do    j = zlo(2), zhi(2)
          do i = zlo(1), zhi(1)
             Ez(i,j,k,1) = sigx1(i)*Ez(i,j,k,1) + sigx2(i)*(By(i  ,j  ,k,1)+By(i  ,j  ,k,2) &
                  &                                        -By(i-1,j  ,k,1)-By(i-1,j  ,k,2))
             Ez(i,j,k,2) = sigy1(j)*Ez(i,j,k,2) - sigy2(j)*(Bx(i  ,j  ,k,1)+Bx(i  ,j  ,k,2) &
                  &                                        -Bx(i  ,j-1,k,1)-Bx(i  ,j-1,k,2))
          end do
       end do
    end do

  end subroutine warpx_push_pml_evec_3d

  subroutine warpx_push_pml_bvec_2d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                             Ex, Exlo, Exhi, &
       &                             Ey, Eylo, Eyhi, &
       &                             Ez, Ezlo, Ezhi, &
       &                             Bx, Bxlo, Bxhi, &
       &                             By, Bylo, Byhi, &
       &                             Bz, Bzlo, Bzhi, &
       &                             sigx1, sigx1_lo, sigx1_hi, &
       &                             sigx2, sigx2_lo, sigx2_hi, &
       &                             sigz1, sigz1_lo, sigz1_hi, &
       &                             sigz2, sigz2_lo, sigz2_hi) &
       bind(c,name='warpx_push_pml_bvec_2d')
    integer, intent(in) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
         Exlo(2), Exhi(2), Eylo(2), Eyhi(2), Ezlo(2), Ezhi(2), &
         Bxlo(2), Bxhi(2), Bylo(2), Byhi(2), Bzlo(2), Bzhi(2)
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(in   ) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),2)
    real(amrex_real), intent(in   ) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),2)
    real(amrex_real), intent(in   ) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),2)
    real(amrex_real), intent(inout) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),2)
    real(amrex_real), intent(inout) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),2)
    real(amrex_real), intent(inout) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),2)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)

    integer :: i, k

    do    k = xlo(2), xhi(2)
       do i = xlo(1), xhi(1)
          Bx(i,k,2) = sigz1(k)*Bx(i,k,2) + sigz2(k)*(Ey(i,k+1,1)+Ey(i,k+1,2) &
               &                                    -Ey(i,k  ,1)-Ey(i,k  ,2))
       end do
    end do

    do    k = ylo(2), yhi(2)
       do i = ylo(1), yhi(1)
          By(i,k,1) = sigz1(k)*By(i,k,1) - sigz2(k)*(Ex(i  ,k+1,1)+Ex(i  ,k+1,2) &
               &                                    -Ex(i  ,k  ,1)-Ex(i  ,k  ,2))
          By(i,k,2) = sigx1(i)*By(i,k,2) + sigx2(i)*(Ez(i+1,k  ,1)+Ez(i+1,k  ,2) &
               &                                    -Ez(i  ,k  ,1)-Ez(i  ,k  ,2))
       end do
    end do

    do    k = zlo(2), zhi(2)
       do i = zlo(1), zhi(1)
          Bz(i,k,1) = sigx1(i)*Bz(i,k,1) - sigx2(i)*(Ey(i+1,k,1)+Ey(i+1,k,2) &
               &                                    -Ey(i  ,k,1)-Ey(i  ,k,2))
       end do
    end do

  end subroutine warpx_push_pml_bvec_2d

  subroutine warpx_push_pml_evec_2d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                             Ex, Exlo, Exhi, &
       &                             Ey, Eylo, Eyhi, &
       &                             Ez, Ezlo, Ezhi, &
       &                             Bx, Bxlo, Bxhi, &
       &                             By, Bylo, Byhi, &
       &                             Bz, Bzlo, Bzhi, &
       &                             sigx1, sigx1_lo, sigx1_hi, &
       &                             sigx2, sigx2_lo, sigx2_hi, &
       &                             sigz1, sigz1_lo, sigz1_hi, &
       &                             sigz2, sigz2_lo, sigz2_hi) &
       bind(c,name='warpx_push_pml_evec_2d')
    integer, intent(in) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
         Exlo(2), Exhi(2), Eylo(2), Eyhi(2), Ezlo(2), Ezhi(2), &
         Bxlo(2), Bxhi(2), Bylo(2), Byhi(2), Bzlo(2), Bzhi(2)
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),2)
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),2)
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),2)
    real(amrex_real), intent(in   ) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),2)
    real(amrex_real), intent(in   ) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),2)
    real(amrex_real), intent(in   ) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),2)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)

    integer :: i, k

    do    k = xlo(2), xhi(2)
       do i = xlo(1), xhi(1)
          Ex(i,k,2) = sigz1(k)*Ex(i,k,2) - sigz2(k)*(By(i,k  ,1)+By(i,k  ,2) &
               &                                    -By(i,k-1,1)-By(i,k-1,2))
       end do
    end do

    do    k = ylo(2), yhi(2)
       do i = ylo(1), yhi(1)
          Ey(i,k,1) = sigz1(k)*Ey(i,k,1) + sigz2(k)*(Bx(i  ,k  ,1)+Bx(i  ,k  ,2) &
               &                                    -Bx(i  ,k-1,1)-Bx(i  ,k-1,2))
          Ey(i,k,2) = sigx1(i)*Ey(i,k,2) - sigx2(i)*(Bz(i  ,k  ,1)+Bz(i  ,k  ,2) &
               &                                    -Bz(i-1,k  ,1)-Bz(i-1,k  ,2))
       end do
    end do

    do    k = zlo(2), zhi(2)
       do i = zlo(1), zhi(1)
          Ez(i,k,1) = sigx1(i)*Ez(i,k,1) + sigx2(i)*(By(i  ,k,1)+By(i  ,k,2) &
               &                                    -By(i-1,k,1)-By(i-1,k,2))
       end do
    end do

  end subroutine warpx_push_pml_evec_2d
  
end module warpx_module
