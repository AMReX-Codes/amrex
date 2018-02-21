
module init_prob_module
  use amrex_fort_module, only : amrex_real
  implicit none
  private
  public :: actual_init_poisson, actual_init_abeclap

contains

  subroutine actual_init_poisson (lo, hi, rhs, rlo, rhi, exact, elo, ehi, prob_lo, prob_hi, dx) &
       bind(c,name='actual_init_poisson')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, elo, ehi
    real(amrex_real), intent(inout) :: rhs  (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(inout) :: exact(elo(1):ehi(1),elo(2):ehi(2),elo(3):ehi(3))
    real(amrex_real), dimension(3), intent(in) :: prob_lo, prob_hi, dx

    integer :: i,j,k
    real(amrex_real) :: x, y, z
    real(amrex_real), parameter :: tpi =  8.d0*atan(1.0)
    real(amrex_real), parameter :: fpi = 16.d0*atan(1.0)
    real(amrex_real), parameter :: fac = tpi*tpi*3.d0

    do k = lo(3), hi(3)
       z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
       do j = lo(2), hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i = lo(1), hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
             
             exact(i,j,k) = 1.d0 * (sin(tpi*x) * sin(tpi*y) * sin(tpi*z))  &
                  &      + .25d0 * (sin(fpi*x) * sin(fpi*y) * sin(fpi*z))
                
             rhs(i,j,k) = -fac * (sin(tpi*x) * sin(tpi*y) * sin(tpi*z))  &
                  &       -fac * (sin(fpi*x) * sin(fpi*y) * sin(fpi*z))
          end do
       end do
    end do

  end subroutine actual_init_poisson


  subroutine actual_init_abeclap (lo, hi, glo, ghi, rhs, rlo, rhi, exact, elo, ehi, &
       alpha, alo, ahi, beta, blo, bhi, a, b, prob_lo, prob_hi, dx) &
       bind(c,name='actual_init_abeclap')
    integer, dimension(3), intent(in) :: lo, hi, glo, ghi, rlo, rhi, elo, ehi, alo, ahi, blo, bhi
    real(amrex_real), intent(inout) :: rhs  (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(inout) :: exact(elo(1):ehi(1),elo(2):ehi(2),elo(3):ehi(3))
    real(amrex_real), intent(inout) :: alpha(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
    real(amrex_real), intent(inout) :: beta (blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))
    real(amrex_real), intent(in), value :: a, b
    real(amrex_real), dimension(3), intent(in) :: prob_lo, prob_hi, dx

    integer :: i,j,k
    real(amrex_real) x, y, z, xc, yc, zc
    real(amrex_real) r, theta, dbdrfac
    real(amrex_real) pi, fpi, tpi, fac
    real(amrex_real), parameter :: w = 0.05d0
    real(amrex_real), parameter :: sigma = 10.d0

    pi = 4.d0 * atan(1.d0)
    tpi = 2.0d0 * pi
    fpi = 4.0d0 * pi
    fac = 12.d0 * pi**2

    xc = (prob_hi(1) + prob_lo(1))/2.d0
    yc = (prob_hi(2) + prob_lo(2))/2.d0
    zc = (prob_hi(3) + prob_lo(3))/2.d0

    theta = 0.5d0*log(3.d0) / (w + 1.d-50)
    
    do k = glo(3), ghi(3)
       z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
       do j = glo(2), ghi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i = glo(1), ghi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
             
             r = sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
             
             beta(i,j,k) = (sigma-1.d0)/2.d0*tanh(theta*(r-0.25d0)) + (sigma+1.d0)/2.d0
          end do
       end do
    end do
    
    do k = lo(3), hi(3)
       z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
       do j = lo(2), hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i = lo(1), hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
             
             r = sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
             
             dbdrfac = (sigma-1.d0)/2.d0/(cosh(theta*(r-0.25d0)))**2 * theta/r
             dbdrfac = dbdrfac * b
             
             alpha(i,j,k) = 1.d0

             exact(i,j,k) = 1.d0 * cos(tpi*x) * cos(tpi*y) * cos(tpi*z)   &
                  &      + .25d0 * cos(fpi*x) * cos(fpi*y) * cos(fpi*z)

             rhs(i,j,k) = beta(i,j,k)*b*fac*(cos(tpi*x) * cos(tpi*y) * cos(tpi*z)   &
                  &                        + cos(fpi*x) * cos(fpi*y) * cos(fpi*z))  &
                  &   + dbdrfac*((x-xc)*(tpi*sin(tpi*x) * cos(tpi*y) * cos(tpi*z)   &
                  &                     + pi*sin(fpi*x) * cos(fpi*y) * cos(fpi*z))  &
                  &            + (y-yc)*(tpi*cos(tpi*x) * sin(tpi*y) * cos(tpi*z)   &
                  &                     + pi*cos(fpi*x) * sin(fpi*y) * cos(fpi*z))  &
                  &            + (z-zc)*(tpi*cos(tpi*x) * cos(tpi*y) * sin(tpi*z)   &
                  &                     + pi*cos(fpi*x) * cos(fpi*y) * sin(fpi*z))) &
                  &                   + a * (cos(tpi*x) * cos(tpi*y) * cos(tpi*z)   &
                  &               + 0.25d0 * cos(fpi*x) * cos(fpi*y) * cos(fpi*z))
          end do
       end do
    end do

  end subroutine actual_init_abeclap

end module init_prob_module
