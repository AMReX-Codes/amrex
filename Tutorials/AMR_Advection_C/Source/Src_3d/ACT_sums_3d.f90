! ::
! :: ----------------------------------------------------------
! :: summass
! ::             MASS = sum{ vol(i,j,k)*rho(i,j,k) }
! :: INPUTS / OUTPUTS:
! ::  rho        => density field
! ::  rlo,rhi    => index limits of rho array
! ::  lo,hi      => index limits of grid interior
! ::  dx         => cell size
! ::  mass      <=  total mass
! :: ----------------------------------------------------------
! ::

      subroutine summass(rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,lo,hi,dx,mass)

        use bl_constants_module

        implicit none

        integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
        integer          :: lo(3), hi(3)
        double precision :: mass, dx(3)
        double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)
        
        integer          :: i, j, k
        double precision :: vol

        vol  = dx(1)*dx(2)*dx(3)
        mass = ZERO

        do k = lo(3), hi(3)
           do i = lo(1), hi(1)
              do j = lo(2), hi(2)
                 mass = mass + rho(i,j,k)
              enddo
           enddo
        enddo

        mass = mass * vol

      end subroutine summass
