
subroutine fort_set_rhs(rhs, lo, hi, nc, dx, a, b, sigma, w) bind(c)
  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in), value :: nc
  double precision, intent(in), value :: a, b, sigma, w
  double precision, intent(inout) :: rhs(lo(1):hi(1),lo(2):hi(2),nc)
  double precision, intent(in) :: dx(2)

  integer :: i, j, n
  double precision :: x, y, r, pi, tpi, fpi, fac, theta, beta, dbdrfac

  pi = 4.d0 * datan(1.d0)
  tpi = 2.0d0 * pi
  fpi = 4.0d0 * pi
  fac = 2.0d0 * tpi**2
  
  theta = 0.5d0*log(3.0) / w

  do n = 1, nc
     do j = lo(2), hi(2)
        y = (dble(j)+0.5d0)*dx(2)
        do i = lo(1), hi(1)
           x = (dble(i)+0.5d0)*dx(1)
           
           r = sqrt((x-0.5d0)**2+(y-0.5d0)**2)

           beta = (sigma-1.d0)/2.d0*tanh(theta*(r-0.25d0)) + (sigma+1.d0)/2.d0
           beta = beta * b
           dbdrfac = (sigma-1.d0)/2.d0/(cosh(theta*(r-0.25d0)))**2 * theta/r
           dbdrfac = dbdrfac * b
        
           rhs(i,j,n) = beta*fac*(sin(tpi*x) * sin(tpi*y)  &
                &               + sin(fpi*x) * sin(fpi*y)) &
                &    + dbdrfac*((x-0.5d0)*(-tpi*cos(tpi*x) * sin(tpi*y)   &
                &                          - pi*cos(fpi*x) * sin(fpi*y))  &
                &             + (y-0.5d0)*(-tpi*sin(tpi*x) * cos(tpi*y)   &
                &                          - pi*sin(fpi*x) * cos(fpi*y))) &
                &    + a * (                    sin(tpi*x) * sin(tpi*y)   &
                &                    + 0.25d0 * sin(fpi*x) * sin(fpi*y))
        end do
     end do
  end do
end subroutine fort_set_rhs

subroutine fort_set_coef(xcoef, xlo, xhi, ycoef, ylo, yhi, nc, dx, sigma, w) bind(c)
  integer, intent(in) :: xlo(2), xhi(2), ylo(2), yhi(2)
  integer, intent(in), value :: nc
  double precision, intent(in), value :: sigma, w
  double precision, intent(inout) :: xcoef(xlo(1):xhi(1),xlo(2):xhi(2),nc)
  double precision, intent(inout) :: ycoef(ylo(1):yhi(1),ylo(2):yhi(2),nc)
  double precision, intent(in) :: dx(2)

  integer :: i, j, n
  double precision :: theta, x, y, r

  theta = 0.5d0*log(3.0) / w

  do n = 1, nc
     do j = xlo(2), xhi(2)
        y = (dble(j)+0.5d0)*dx(2)
        do i = xlo(1), xhi(1)
           x = (dble(i))*dx(1)
           
           r = sqrt((x-0.5d0)**2 + (y-0.5d0)**2)

           xcoef(i,j,n) = (sigma-1.d0)/2.d0*tanh(theta*(r-0.25d0)) &
                + (sigma+1.d0)/2.d0
        end do
     end do

     do j = ylo(2), yhi(2)
        y = (dble(j))*dx(2)
        do i = ylo(1), yhi(1)
           x = (dble(i)+0.5d0)*dx(1)
           
           r = sqrt((x-0.5d0)**2 + (y-0.5d0)**2)

           ycoef(i,j,n) = (sigma-1.d0)/2.d0*tanh(theta*(r-0.25d0)) &
                + (sigma+1.d0)/2.d0
        end do
     end do
  end do
end subroutine fort_set_coef
