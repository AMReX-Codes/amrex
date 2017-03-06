#include <AMReX_REAL.H>
#include "AMReX_LO_BCTYPES.H"


subroutine comp_asol(asol, a_l1, a_l2, a_h1, a_h2, &
                     lo, hi, dx, ibnd, offset) bind(C, name="comp_asol")

  implicit none

  integer          :: lo(2), hi(2), ibnd
  integer          :: a_l1, a_l2, a_h1, a_h2
  double precision :: asol(a_l1:a_h1, a_l2:a_h2)
  double precision :: dx(2), offset(2)

  integer i,j
  REAL_T  x,y
  REAL_T  pi, fpi, tpi

  pi = 4.d0 * datan(1.d0)

  tpi = 2.0d0 * pi
  fpi = 4.0d0 * pi

  do j = a_l2, a_h2
     y = (dble(j)+offset(2))*dx(2)
     do i = a_l1, a_h1
        x = (dble(i)+offset(1))*dx(1)

        if (ibnd .eq. 0 .or. ibnd.eq. LO_NEUMANN) then
           asol(i,j) = 1.d0 * cos(tpi*x) * cos(tpi*y) &
                    + .25d0 * cos(fpi*x) * cos(fpi*y)
        else if (ibnd .eq. LO_DIRICHLET) then
           asol(i,j) = 1.d0 * sin(tpi*x) * sin(tpi*y) &
                    + .25d0 * sin(fpi*x) * sin(fpi*y)
        else
           print *, 'FORT_COMP_ASOL: unknown boundary type'
           stop
        endif
     end do
  end do
 
end subroutine comp_asol


subroutine set_alpha(alpha, a_l1, a_l2, a_h1, a_h2, lo, hi, dx) bind(C, name="set_alpha")

  implicit none

  integer          :: lo(2), hi(2)
  integer          :: a_l1, a_l2, a_h1, a_h2
  double precision :: alpha(a_l1:a_h1, a_l2:a_h2)
  double precision :: dx(2)
  
  integer i,j

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        alpha(i,j) = 1.d0
     end do
  end do

end subroutine set_alpha


subroutine set_cc_coef(coef, c_l1, c_l2, c_h1, c_h2, &
                       lo, hi, dx, sigma, w) bind(C, name="set_cc_coef")
  
  implicit none

  integer          :: lo(2), hi(2)
  integer          :: c_l1, c_l2, c_h1, c_h2
  double precision :: coef(c_l1:c_h1, c_l2:c_h2)
  double precision :: dx(2), sigma, w
 
  integer i,j
  double precision theta, x, y, r
 
  theta = 0.5d0*log(3.0) / (w + 1.d-50)
      
  do j = c_l2, c_h2
     y = (dble(j)+0.5d0)*dx(2)
     do i = c_l1, c_h1
        x = (dble(i)+0.5d0)*dx(1)
            
        r = sqrt((x-0.5d0)**2 + (y-0.5d0)**2)
            
        coef(i,j) = (sigma-1.d0)/2.d0*tanh(theta*(r-0.25d0)) &
                  + (sigma+1.d0)/2.d0
     end do
  end do
  
end subroutine set_cc_coef

