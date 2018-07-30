
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
   
    real(amrex_real) :: dxb, dyb, denom, x, y, dyinv, dxinv, term1, term2, b
    integer :: i, j 
   
    dxinv   = 1.d0/dx(1)
    dyinv   = 1.d0/dx(2)
    do j    = lo(2), hi(2) 
      y     = problo(2) + (dble(j) + 0.5d0)*dx(2) -0.5d0
      dyb   = (by(i,j+1) - by(i,j))*dyinv
      do i  = lo(1), hi(1)

        b  = 0.25d0*(bx(i+1,j) + bx(i,j) + by(i,j+1) + by(i,j))
        dxb = (bx(i+1,j) - bx(i,j))*dxinv
        x   = problo(1) + (dble(i) + 0.5d0)*dx(1) - 0.5d0; 
        denom = sqrt(x*x + y*y) 
        term1 = a(i,j)*x/denom
        term2 = (y*(y*dxb -x*dyb) - x*b)/(denom**3)
        rhs(i,j) = term1 - term2
      enddo
    enddo             
  end subroutine build_rhs
end module rhs
