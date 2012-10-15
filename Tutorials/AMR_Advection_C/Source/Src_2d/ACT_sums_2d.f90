
! :: ----------------------------------------------------------
! :: summass
! ::             MASS = sum{ vol(i,j)*rho(i,j) }
! ::
! :: INPUTS / OUTPUTS:
! ::  rho        => density field
! ::  rlo,rhi    => index limits of rho aray
! ::  lo,hi      => index limits of grid interior
! ::  dx         => cell size
! ::  mass      <=  total mass
! ::  tmp        => temp column array
! :: ----------------------------------------------------------
! ::

       subroutine summass(rho,r_l1,r_l2,r_h1,r_h2,lo,hi,dx,mass)

       use prob_params_module, only : coord_type
  
       implicit none
       integer r_l1,r_l2,r_h1,r_h2
       integer lo(2), hi(2)
       double precision mass, dx(2)
       double precision rho(r_l1:r_h1,r_l2:r_h2)
       double precision  tmp(lo(2):hi(2))

       integer i, j
       double precision vol

       vol = dx(1) * dx(2)

       do j = lo(2),hi(2)
          tmp(j) = 0.d0
       enddo

       do i = lo(1), hi(1)
          do j = lo(2), hi(2)
             tmp(j) = tmp(j) + vol * rho(i,j)
          enddo
       enddo

       mass = 0.d0
       do j = lo(2), hi(2)
          mass = mass + tmp(j)
       enddo

       end subroutine summass

