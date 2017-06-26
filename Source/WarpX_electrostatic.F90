module warpx_electrostatic_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real

  implicit none

contains

  subroutine warpx_sum_fine_to_crse_nodal_3d (lo, hi, lrat, crse, clo, chi, fine, flo, fhi) &
       bind(c, name="warpx_sum_fine_to_crse_nodal_3d")

    integer, intent(in)             ::   lo(3),  hi(3)
    integer, intent(in)             ::  clo(3), chi(3)
    integer, intent(in)             ::  flo(3), fhi(3)
    integer, intent(in)             ::  lrat(3)
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(amrex_real), intent(in)    :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    
    integer :: i, j, k, ii, jj, kk
    
    do k        = lo(3), hi(3)
       kk       = k * lrat(3)
       do j     = lo(2), hi(2)
          jj    = j * lrat(2)
          do i  = lo(1), hi(1)
             ii = i * lrat(1)
             crse(i,j,k)  =  fine(ii,jj,kk)                              + &
! These six fine nodes are shared by two coarse nodes...
                  0.5d0   * (fine(ii-1,jj,kk)     + fine(ii+1,jj,kk)     + & 
                             fine(ii,jj-1,kk)     + fine(ii,jj+1,kk)     + &
                             fine(ii,jj,kk-1)     + fine(ii,jj,kk+1))    + &
! ... these twelve are shared by four...
                  0.25d0  * (fine(ii,jj-1,kk-1)   + fine(ii,jj+1,kk-1)   + &
                             fine(ii,jj-1,kk+1)   + fine(ii,jj+1,kk+1)   + &
                             fine(ii-1,jj,kk-1)   + fine(ii+1,jj,kk-1)   + &
                             fine(ii-1,jj,kk+1)   + fine(ii+1,jj,kk+1)   + &
                             fine(ii-1,jj-1,kk)   + fine(ii+1,jj-1,kk)   + &
                             fine(ii-1,jj+1,kk)   + fine(ii+1,jj+1,kk))  + &
! ... and these eight are shared by eight
                  0.125d0 * (fine(ii-1,jj-1,kk-1) + fine(ii-1,jj-1,kk+1) + &
                             fine(ii-1,jj+1,kk-1) + fine(ii-1,jj+1,kk+1) + &
                             fine(ii+1,jj-1,kk-1) + fine(ii+1,jj-1,kk+1) + &
                             fine(ii+1,jj+1,kk-1) + fine(ii+1,jj+1,kk+1))
! ... note that we have 27 nodes in total...
             crse(i,j,k) = crse(i,j,k) / 8.d0
          end do
       end do
    end do

  end subroutine warpx_sum_fine_to_crse_nodal_3d

  subroutine warpx_sum_fine_to_crse_nodal_2d (lo, hi, lrat, crse, clo, chi, fine, flo, fhi) &
       bind(c, name="warpx_sum_fine_to_crse_nodal_2d")

    integer, intent(in)             ::   lo(2),  hi(2)
    integer, intent(in)             ::  clo(2), chi(2)
    integer, intent(in)             ::  flo(2), fhi(2)
    integer, intent(in)             ::  lrat(2)
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2))
    real(amrex_real), intent(in)    :: fine(flo(1):fhi(1),flo(2):fhi(2))
    
    integer :: i, j, ii, jj
    
    do j     = lo(2), hi(2)
       jj    = j * lrat(2)
       do i  = lo(1), hi(1)
          ii = i * lrat(1)
          crse(i,j)  =  fine(ii,jj)                             + &
! These four fine nodes are shared by two coarse nodes...
               0.5d0   * (fine(ii-1,jj)     + fine(ii+1,jj)     + & 
               fine(ii,jj-1)     + fine(ii,jj+1))               + &
! ... and these four are shared by four...
               0.25d0  * (fine(ii-1,jj-1)   + fine(ii-1,jj+1)   + &
               fine(ii-1,jj+1)   + fine(ii+1,jj+1))
! ... note that we have 9 nodes in total...
             crse(i,j) = crse(i,j) / 4.d0
       end do
    end do

  end subroutine warpx_sum_fine_to_crse_nodal_2d

  subroutine warpx_zero_out_bndry_3d (lo, hi, input_data, bndry_data, mask) &
       bind(c,name='warpx_zero_out_bndry_3d')

    integer(c_int),   intent(in   ) :: lo(3), hi(3)
    double precision, intent(inout) :: input_data(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision, intent(inout) :: bndry_data(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
    integer(c_int),   intent(in   ) :: mask (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    
    integer :: i, j, k
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)          
             if (mask(i,j,k) .eq. 1) then
                bndry_data(i,j,k) = input_data(i,j,k)
                input_data(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

  end subroutine warpx_zero_out_bndry_3d

  subroutine warpx_zero_out_bndry_2d (lo, hi, input_data, bndry_data, mask) &
       bind(c,name='warpx_zero_out_bndry_2d')

    integer(c_int),   intent(in   ) :: lo(2), hi(2)
    double precision, intent(inout) :: input_data(lo(1):hi(1),lo(2):hi(2))
    double precision, intent(inout) :: bndry_data(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)
    integer(c_int),   intent(in   ) :: mask (lo(1):hi(1),lo(2):hi(2))
    
    integer :: i, j
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)          
          if (mask(i,j) .eq. 1) then
             bndry_data(i,j) = input_data(i,j)
             input_data(i,j) = 0.d0
          end if
       end do
    end do

  end subroutine warpx_zero_out_bndry_2d

  subroutine warpx_build_mask_3d (lo, hi, tmp_mask, mask, ncells) &
       bind(c,name='warpx_build_mask_3d')
    integer(c_int),   intent(in   ) :: lo(3), hi(3)
    integer(c_int),   intent(in   ) :: ncells
    integer(c_int),   intent(in   ) :: tmp_mask(lo(1)-ncells:hi(1)+ncells,lo(2)-ncells:hi(2)+ncells,lo(3)-ncells:hi(3)+ncells)
    integer(c_int),   intent(inout) :: mask (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    integer :: i, j, k, ii, jj, kk, total

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             total = 0
             do ii = i-ncells, i+ncells
                do jj = j-ncells, j+ncells
                   do kk = k-ncells, k+ncells
                      total = total + tmp_mask(ii, jj, kk)
                   end do
                end do
             end do
             
             if (total .gt. 0) then
                mask(i,j,k) = 1
             else
                mask(i,j,k) = 0
             end if

          end do
       end do
    end do

  end subroutine warpx_build_mask_3d

  subroutine warpx_build_mask_2d (lo, hi, tmp_mask, mask, ncells) &
       bind(c,name='warpx_build_mask_2d')
    integer(c_int),   intent(in   ) :: lo(2), hi(2)
    integer(c_int),   intent(in   ) :: ncells
    integer(c_int),   intent(in   ) :: tmp_mask(lo(1)-ncells:hi(1)+ncells,lo(2)-ncells:hi(2)+ncells)
    integer(c_int),   intent(inout) :: mask (lo(1):hi(1),lo(2):hi(2))

    integer :: i, j, ii, jj, total

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          
          total = 0
          do ii = i-ncells, i+ncells
             do jj = j-ncells, j+ncells
                total = total + tmp_mask(ii, jj)
             end do
          end do
             
          if (total .gt. 0) then
             mask(i,j) = 1
          else
             mask(i,j) = 0
          end if

       end do
    end do

  end subroutine warpx_build_mask_2d

! This routine computes the node-centered electric field given a node-centered phi.
! The gradient is computed using 2nd-order centered differences. It assumes the 
! Boundary conditions have already been set and that you have two rows of ghost cells
! for phi and one row of ghost cells for Ex, Ey, and Ez.
! Note that this routine includes the minus sign in E = - grad phi.
!
! Arguments:
!     lo, hi:     The corners of the valid box over which the gradient is taken
!     Ex, Ey, Ez: The electric field in the x, y, and z directions.
!     dx:         The cell spacing
!
  subroutine warpx_compute_E_nodal_3d (lo, hi, phi, Ex, Ey, Ez, dx) &
       bind(c,name='warpx_compute_E_nodal_3d')
    integer(c_int),   intent(in)    :: lo(3), hi(3)
    real(amrex_real), intent(in)    :: dx(3)
    real(amrex_real), intent(in   ) :: phi(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2)
    real(amrex_real), intent(inout) :: Ex (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
    real(amrex_real), intent(inout) :: Ey (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
    real(amrex_real), intent(inout) :: Ez (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

    integer :: i, j, k
    real(amrex_real) :: fac(3)

    fac = 0.5d0 / dx

    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             
             Ex(i,j,k) = fac(1) * (phi(i-1,j,k) - phi(i+1,j,k))
             Ey(i,j,k) = fac(2) * (phi(i,j-1,k) - phi(i,j+1,k))
             Ez(i,j,k) = fac(3) * (phi(i,j,k-1) - phi(i,j,k+1))

          end do
       end do
    end do

  end subroutine warpx_compute_E_nodal_3d

  subroutine warpx_compute_E_nodal_2d (lo, hi, phi, Ex, Ey, dx) &
       bind(c,name='warpx_compute_E_nodal_2d')
    integer(c_int),   intent(in)    :: lo(2), hi(2)
    real(amrex_real), intent(in)    :: dx(2)
    real(amrex_real), intent(in   ) :: phi(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2)
    real(amrex_real), intent(inout) :: Ex (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)
    real(amrex_real), intent(inout) :: Ey (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)

    integer :: i, j
    real(amrex_real) :: fac(2)

    fac = 0.5d0 / dx

    do j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)+1
             
          Ex(i,j) = fac(1) * (phi(i-1,j) - phi(i+1,j))
          Ey(i,j) = fac(2) * (phi(i,j-1) - phi(i,j+1))

       end do
    end do

  end subroutine warpx_compute_E_nodal_2d


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
    real(amrex_real), intent(in)     :: particles(ns,np)
    real(amrex_real), intent(in)     :: weights(np)
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
    real(amrex_real), intent(in)     :: particles(ns,np)
    real(amrex_real), intent(in)     :: weights(np)
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
    real(amrex_real), intent(in)     :: particles(ns,np)
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
    real(amrex_real), intent(in)     :: particles(ns,np)
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
    real(amrex_real), intent(in)     :: particles(ns,np)
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
    real(amrex_real), intent(in)     :: particles(ns,np)
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


!
! This routine updates the particle positions and velocities using the 
! leapfrog time integration algorithm, given the electric fields at the
! particle positions. It also enforces specular reflection off the domain
! walls.
!
! Arguments:
!     particles : a pointer to the particle array-of-structs 
!     ns        : the stride length of particle struct (the size of the struct in number of reals)
!     np        : the number of particles
!     vx_p      : the particle x-velocities
!     vy_p      : the particle y-velocities
!     vz_p      : the particle z-velocities
!     Ex_p      : the electric field in the x-direction at the particle positions
!     Ey_p      : the electric field in the y-direction at the particle positions
!     Ez_p      : the electric field in the z-direction at the particle positions
!     charge    : the charge of this particle species
!     mass      : the mass of this particle species
!     dt        : the time step
!     prob_lo   : the left-hand corner of the problem domain
!     prob_hi   : the right-hand corner of the problem domain
!
  subroutine warpx_push_leapfrog_3d(particles, ns, np,      &
                                    vx_p, vy_p, vz_p,       &                                 
                                    Ex_p, Ey_p, Ez_p,       &
                                    charge, mass, dt,       &
                                    prob_lo, prob_hi)       &
       bind(c,name='warpx_push_leapfrog_3d')
    integer, value,   intent(in)     :: ns, np
    real(amrex_real), intent(inout)  :: particles(ns,np)
    real(amrex_real), intent(inout)  :: vx_p(np), vy_p(np), vz_p(np)
    real(amrex_real), intent(in)     :: Ex_p(np), Ey_p(np), Ez_p(np)
    real(amrex_real), intent(in)     :: charge
    real(amrex_real), intent(in)     :: mass
    real(amrex_real), intent(in)     :: dt
    real(amrex_real), intent(in)     :: prob_lo(3), prob_hi(3)
   
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

  end subroutine warpx_push_leapfrog_3d

  subroutine warpx_push_leapfrog_2d(particles, ns, np,      &
                                    vx_p, vy_p,             &
                                    Ex_p, Ey_p,             &
                                    charge, mass, dt,       &
                                    prob_lo, prob_hi)       &
       bind(c,name='warpx_push_leapfrog_2d')
    integer, value,   intent(in)     :: ns, np
    real(amrex_real), intent(inout)  :: particles(ns,np)
    real(amrex_real), intent(inout)  :: vx_p(np), vy_p(np)
    real(amrex_real), intent(in)     :: Ex_p(np), Ey_p(np)
    real(amrex_real), intent(in)     :: charge
    real(amrex_real), intent(in)     :: mass
    real(amrex_real), intent(in)     :: dt
    real(amrex_real), intent(in)     :: prob_lo(2), prob_hi(2)
   
    integer n
    real(amrex_real) fac

    fac = charge * dt / mass

    do n = 1, np

       vx_p(n) = vx_p(n) + fac * Ex_p(n)
       vy_p(n) = vy_p(n) + fac * Ey_p(n)

       particles(1, n) = particles(1, n) + dt * vx_p(n)
       particles(2, n) = particles(2, n) + dt * vy_p(n)

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

    end do

  end subroutine warpx_push_leapfrog_2d


!
! This routine advances the particle positions using the current
! velocity. This is needed to desynchronize the particle positions
! from the velocities after particle initialization.
!
! Arguments:
!     particles : a pointer to the particle array-of-structs 
!     ns        : the stride length of particle struct (the size of the struct in number of reals)
!     np        : the number of particles
!     xx_p      : the electric field in the x-direction at the particle positions
!     vy_p      : the electric field in the y-direction at the particle positions
!     vz_p      : the electric field in the z-direction at the particle positions
!     dt        : the time step
!     prob_lo   : the left-hand corner of the problem domain
!     prob_hi   : the right-hand corner of the problem domain
!
  subroutine warpx_push_leapfrog_positions_3d(particles, ns, np,      &
                                             vx_p, vy_p, vz_p, dt,   &
                                             prob_lo, prob_hi)       &
       bind(c,name='warpx_push_leapfrog_positions_3d')
    integer, value,   intent(in)    :: ns, np
    real(amrex_real), intent(inout) :: particles(ns,np)
    real(amrex_real), intent(inout) :: vx_p(np), vy_p(np), vz_p(np)
    real(amrex_real), intent(in)    :: dt
    real(amrex_real), intent(in)    :: prob_lo(3), prob_hi(3)

    integer n

    do n = 1, np

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

  end subroutine warpx_push_leapfrog_positions_3d

  subroutine warpx_push_leapfrog_positions_2d(particles, ns, np,      &
                                              vx_p, vy_p, dt,         &
                                              prob_lo, prob_hi)       &
       bind(c,name='warpx_push_leapfrog_positions_2d')
    integer, value,   intent(in)    :: ns, np
    real(amrex_real), intent(inout) :: particles(ns,np)
    real(amrex_real), intent(inout) :: vx_p(np), vy_p(np)
    real(amrex_real), intent(in)    :: dt
    real(amrex_real), intent(in)    :: prob_lo(2), prob_hi(2)

    integer n

    do n = 1, np

       particles(1, n) = particles(1, n) + dt * vx_p(n)
       particles(2, n) = particles(2, n) + dt * vy_p(n)

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

    end do

  end subroutine warpx_push_leapfrog_positions_2d

end module electrostatic_pic_module
