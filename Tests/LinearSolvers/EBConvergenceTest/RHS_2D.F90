
module rhs
  
private
public :: build_rhs
contains 
  subroutine build_rhs(lo, hi, rhs, rlo, rhi, &
                               a  , alo, ahi, &
                             bx , bxlo, bxhi, &
                             by , bylo, byhi, & 
                          dx, problo, probhi) &
                          bind(c, name='build_rhs')
  use amrex_fort_module, only : amrex_real
  implicit none
    integer, dimension(2), intent(in) :: lo, hi, rlo, rhi, alo, ahi, bxlo, bxhi, bylo, byhi
    real(amrex_real), dimension(2), intent(in) :: dx, problo, probhi
    real(amrex_real), intent(inout) :: rhs( rlo(1): rhi(1),  rlo(2): rhi(2))
    real(amrex_real), intent(in   ) :: a  ( alo(1): ahi(1),  alo(2): ahi(2))
    real(amrex_real), intent(in   ) :: bx (bxlo(1):bxhi(1), bxlo(2):bxhi(2))
    real(amrex_real), intent(in   ) :: by (bylo(1):byhi(1), bylo(2):byhi(2))
   
    real(amrex_real) :: dxb, dyb, denom, x, y, dyinv, dxinv, term1, term2, b, xl, xh, yl, yh
    integer :: i, j 
   
    dxinv   = 1.d0/dx(1)
    dyinv   = 1.d0/dx(2)
    yl = problo(2)
    xl = problo(1)
    yh = probhi(2)
    xh = probhi(1) 
    do j    = lo(2)+1, hi(2)-1
      y     = yl + (dble(j) + 0.5d0)*dx(2) -0.5d0
      do i  = lo(1)+1, hi(1)-1
        !Inside Domain
        dyb   = (by(i,j+1) - by(i,j))*dyinv
        b  = 0.25d0*(bx(i+1,j) + bx(i,j) + by(i,j+1) + by(i,j))
        dxb = (bx(i+1,j) - bx(i,j))*dxinv
        x   = xl + (dble(i) + 0.5d0)*dx(1) - 0.5d0; 
        denom = sqrt(x*x + y*y) 
        term1 = a(i,j)*x/denom
        term2 = (y*(y*dxb -x*dyb) - x*b)/(denom**3)
        rhs(i,j) = term1 - term2
      enddo
      !Boundary Case xlo and xhi, y inside
      b   = 0.25d0*(bx(lo(1)+1,j) + bx(lo(1),j) + by(lo(1),j+1) + by(lo(1),j))
      dyb = (by(lo(1), j+1) - by(lo(1),j))*dyinv
      x = xl - 0.5d0
      dxb = (bx(lo(1)+1,j)-bx(lo(1),j))*dxinv
      denom = sqrt(x**2 + y*y)
      term1 = a(lo(1),j)*x/denom
      term2 = (y*(y*dxb - x*dyb) - x*b)/(denom**3)
      rhs(lo(1),j) = term1 + term2 

      dyb = (by(hi(1), j+1) - by(hi(1), j))*dyinv
      x = xh - 0.5d0 
      b   = 0.25d0*(bx(hi(1)+1,j) + bx(hi(1),j) + by(hi(1),j+1) + by(hi(1),j))
      dxb = (bx(hi(1)+1, j) - bx(hi(1),j))*dxinv
      denom = sqrt((x)**2 + y*y)
      term1 = a(hi(1),j)*x/denom
      term2 = (y*(y*dxb - x*dyb) - x*b)/(denom**3)
      rhs(hi(1),j) = term1 - term2 
    enddo   

    !Boundary case ylo, yhi, x inside
    do i = lo(1)+1, hi(1)-1
      x = xl + (dble(i) + 0.5d0)*dx(1) - 0.5d0; 
     
      y = yl - 0.5d0 
      b   = 0.25d0*(bx(i+1,lo(2)) + bx(i,lo(2)) + by(i,lo(2)+1) + by(i,lo(2)))
      dyb = (by(i, lo(2)+1) - by(i,lo(2)))*dyinv
      dxb = (bx(i+1, lo(2)) - bx(i,lo(2)))*dxinv
      denom = sqrt(x**2 + y**2)
      term1 = a(i,lo(2))*x/denom
      term2 = (y*(y*dxb - x*dyb) - x*b)/(denom**3)
      rhs(i,lo(2)) = term1 - term2 

      y = yh -0.5d0
      b   = 0.25d0*(bx(i+1,hi(2)) + bx(i,hi(2)) + by(i,hi(2)+1) + by(i,hi(2)))
      dyb = (by(i,hi(2)+1) - by(i,hi(2)))*dyinv
      dxb = (bx(i+1,hi(2)) - bx(i,hi(2)))*dxinv
      denom = sqrt(x**2 + y**2)
      term1 = a(i,hi(2))*x/denom
      term2 = (y*(y*dxb - x*dyb) - x*b)/(denom**3)
      rhs(i,hi(2)) = term1 - term2      
    enddo          

    !Corner cases
    !xlo ylo 
    b   = 0.25d0*(bx(lo(1)+1,lo(2)) + bx(lo(1),lo(2)) + by(lo(1),lo(2)+1) + by(lo(1),lo(2)))
    x = xl - 0.5d0
    y = yl - 0.5d0 
    dxb = (bx(lo(1)+1, lo(2)) - bx(lo(1), lo(2)))*dxinv
    dyb = (by(lo(1), lo(2)+1) - by(lo(1), lo(2)))*dyinv
    denom = sqrt(x*x + y*y)
    term1 = a(lo(1),lo(2))*x/denom
    term2 = (y*(y*dxb - x*dyb) - x*b)/(denom**3)
    rhs(lo(1), lo(2)) = term1 - term2 
  
    !xlo yhi
    b   = 0.25d0*(bx(lo(1)+1,hi(2)) + bx(lo(1),hi(2)) + by(lo(1),hi(2)+1) + by(lo(1),hi(2)))
    y = yh -0.5d0
    dxb = (bx(lo(1)+1, hi(2)) - bx(lo(1), hi(2)))*dxinv
    dyb = (by(lo(1), hi(2)+1) - by(lo(1), hi(2)))*dyinv
    denom = sqrt(x*x + y*y)
    term1 = a(lo(1), hi(2))*x/denom
    term2 = (y*(y*dxb - x*dyb) - x*b)/(denom**3)
    rhs(lo(1), hi(2)) = term1 - term2

    !xhi yhi 
    b   = 0.25d0*(bx(hi(1)+1,hi(2)) + bx(hi(1),hi(2)) + by(hi(1),hi(2)+1) + by(hi(1),hi(2)))
    x = xh -0.5d0
    dxb = (bx(hi(1)+1, hi(2)) - bx(hi(1), hi(2)))*dxinv
    dyb = (by(hi(1), hi(2)+1) - by(hi(1), hi(2)))*dyinv
    denom = sqrt(x*x + y*y)
    term1 = a(hi(1), hi(2))*x/denom
    term2 = (y*(y*dxb - x*dyb) - x*b)/(denom**3)
    rhs(hi(1), hi(2)) = term1 - term2

    !xhi ylo 
    b   = 0.25d0*(bx(hi(1)+1,lo(2)) + bx(hi(1),lo(2)) + by(hi(1),lo(2)+1) + by(hi(1),lo(2)))
    y = yl -0.5d0
    dxb = (bx(hi(1)+1, lo(2)) - bx(hi(1), lo(2)))*dxinv
    dyb = (by(hi(1), lo(2)+1) - by(hi(1), lo(2)))*dyinv
    denom = sqrt(x*x + y*y)
    term1 = a(hi(1), lo(2))*x/denom
    term2 = (y*(y*dxb - x*dyb) - x*b)/(denom**3)
    rhs(hi(1), lo(2)) = term1 - term2
  
  end subroutine build_rhs
end module rhs
