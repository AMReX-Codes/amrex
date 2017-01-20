
module stencil_test_module

  implicit none

  private

  public :: lapleb_MSD, lapl_MSD

contains

  subroutine lapleb_MSD(&
       lph,lph_l1,lph_l2,lph_h1,lph_h2, &
       phi,phi_l1,phi_l2,phi_h1,phi_h2, &
       fd0,fd0_l1,fd0_l2,fd0_h1,fd0_h2, &
       fd1,fd1_l1,fd1_l2,fd1_h1,fd1_h2, &
       lo, hi, dx) bind(C,name="lapleb_MSD")

    integer, intent(in) :: lph_l1,lph_l2,lph_h1,lph_h2
    integer, intent(in) :: phi_l1,phi_l2,phi_h1,phi_h2
    integer, intent(in) :: fd0_l1,fd0_l2,fd0_h1,fd0_h2
    integer, intent(in) :: fd1_l1,fd1_l2,fd1_h1,fd1_h2

    integer, intent(in) :: lo(2),hi(2)

    double precision, intent(inout) :: lph(lph_l1:lph_h1,lph_l2:lph_h2)
    double precision, intent(in   ) :: phi(phi_l1:phi_h1,phi_l2:phi_h2)
    double precision, intent(in   ) :: fd0(fd0_l1:fd0_h1,fd0_l2:fd0_h2,3)
    double precision, intent(in   ) :: fd1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,3)
    double precision, intent(in   ) :: dx(2)

    integer :: i,j,ip,im,jp,jm, in, jn, ii, jj, iface, jface
    double precision :: fac1, fac2, cent, f(-1:1), Fx(0:1), Fy(0:1)

    fac1 = -1.d0 / (dx(1)*dx(1))
    fac2 = -1.d0 / (dx(2)*dx(2))

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          do ii=0,1
             iface = i+ii
             jp = j+1
             jm = j-1
             f(-1) = fac1*(phi(iface,jm) - phi(iface-1,jm))
             f( 0) = fac1*(phi(iface,j ) - phi(iface-1,j ))
             f(+1) = fac1*(phi(iface,jp) - phi(iface-1,jp))
             cent = ABS(fd0(iface,j,3))
             jn = SIGN(1.d0,fd0(iface,j,3))
             Fx(ii) = (f(0)*(1.d0 - cent) + f(jn)*cent)*fd0(iface,j,1)
          enddo

          do jj=0,1
             jface = j+jj
             ip = i+1
             im = i-1
             f(-1) = fac2*(phi(im,jface) - phi(im,jface-1))
             f( 0) = fac2*(phi(i ,jface) - phi(i ,jface-1))
             f(+1) = fac2*(phi(ip,jface) - phi(ip,jface-1))
             cent = ABS(fd0(i,jface,2))
             in = SIGN(1.d0,fd0(i,jface,2))
             Fy(jj) = (f(0)*(1.d0 - cent) + f(in)*cent)*fd1(i,jface,1)
          enddo

          lph(i,j) = Fx(1) - Fx(0) + Fy(1) - Fy(0)

       enddo
    enddo

  end subroutine lapleb_MSD

  subroutine lapl_MSD(&
       lph,lph_l1,lph_l2,lph_h1,lph_h2, &
       phi,phi_l1,phi_l2,phi_h1,phi_h2, &
       lo, hi, dx) bind(C,name="lapl_MSD")

    integer, intent(in) :: lph_l1,lph_l2,lph_h1,lph_h2
    integer, intent(in) :: phi_l1,phi_l2,phi_h1,phi_h2

    integer, intent(in) :: lo(2),hi(2)

    double precision, intent(inout) :: lph(lph_l1:lph_h1,lph_l2:lph_h2)
    double precision, intent(in   ) :: phi(phi_l1:phi_h1,phi_l2:phi_h2)
    double precision, intent(in   ) :: dx(2)

    integer :: i,j
    double precision :: fac1, fac2

    fac1 = 1.d0 / (dx(1)*dx(1))
    fac2 = 1.d0 / (dx(2)*dx(2))

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          lph(i,j) = fac1*(phi(i+1,j) + phi(i-1,j) - 2*phi(i,j)) &
               +     fac2*(phi(i,j+1) + phi(i,j-1) - 2*phi(i,j))
       enddo
    enddo

  end subroutine lapl_MSD

end module stencil_test_module
