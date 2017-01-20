
module stencil_test_module

  implicit none

  private

  public :: lapl_MSD

contains

  subroutine lapleb_MSD(&
       lph,lph_l1,lph_l2,lph_l3,lph_h1,lph_h2,lph_h3, &
       phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
       fd0,fd0_l1,fd0_l2,fd0_l3,fd0_h1,fd0_h2,fd0_h3, &
       fd1,fd1_l1,fd1_l2,fd1_l3,fd1_h1,fd1_h2,fd1_h3, &
       fd2,fd2_l1,fd2_l2,fd2_l3,fd2_h1,fd2_h2,fd2_h3, &
       lo, hi, dx) bind(C,name="lapleb_MSD")

    integer, intent(in) :: lph_l1,lph_l2,lph_l3,lph_h1,lph_h2,lph_h3
    integer, intent(in) :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
    integer, intent(in) :: fd0_l1,fd0_l2,fd0_l3,fd0_h1,fd0_h2,fd0_h3
    integer, intent(in) :: fd1_l1,fd1_l2,fd1_l3,fd1_h1,fd1_h2,fd1_h3
    integer, intent(in) :: fd2_l1,fd2_l2,fd2_l3,fd2_h1,fd2_h2,fd2_h3

    integer, intent(in) :: lo(3),hi(3)
                                                         
    double precision, intent(inout) :: lph(lph_l1:lph_h1,lph_l2:lph_h2,lph_l3:lph_h3)
    double precision, intent(in   ) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
    double precision, intent(in   ) :: fd0(fd0_l1:fd0_h1,fd0_l2:fd0_h2,fd0_l3:fd0_h3,4)
    double precision, intent(in   ) :: fd1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,fd1_l3:fd1_h3,4)
    double precision, intent(in   ) :: fd2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,fd2_l3:fd2_h3,4)
    double precision, intent(in   ) :: dx(3)

    integer :: i,j,k,ip,im,jp,jm,kp,km, in, jn, kn, ii, jj, kk, iface, jface, kface, is, js, ks
    double precision :: fac0, fac1, fac2, c0, c1, f00, f01, f10, f11, Fx(0:1), Fy(0:1), Fz(0:1)

    fac0 = -1.d0 / (dx(1)*dx(1))
    fac1 = -1.d0 / (dx(2)*dx(2))
    fac2 = -1.d0 / (dx(3)*dx(3))

    do i=lo(1),hi(1)
       do j=lo(2),hi(2)
          do k=lo(3),hi(3)

             do ii=0,1
                iface = i+ii

                js = SIGN(1.d0, fd0(iface,j,k,3))
                ks = SIGN(1.d0, fd0(iface,j,k,4))

                jn = j + js
                kn = k + ks

                f00 = fac0*(phi(iface,j ,k ) - phi(iface-1,j ,k ))
                f10 = fac0*(phi(iface,jn,k ) - phi(iface-1,jn,k ))*js
                f01 = fac0*(phi(iface,j ,kn) - phi(iface-1,j ,kn))*ks
                f11 = fac0*(phi(iface,jn,kn) - phi(iface-1,jn,kn))*js*ks

                c0 = ABS(fd0(iface,j,k,3))
                c1 = ABS(fd0(iface,j,k,4))

                Fx(ii) = (f00*(1.d0-c0)*(1.d0-c1) + f10*c0*(1.d0-c1) + f01*(1.d0-c0)*c1 + f11*c0*c1)*fd0(iface,j,k,1)
             enddo

             do jj=0,1
                jface = j+jj

                is = SIGN(1.d0, fd1(i,jface,k,2))
                ks = SIGN(1.d0, fd1(i,jface,k,4))

                in = i + is
                kn = k + ks

                f00 = fac1*(phi(i ,jface,k ) - phi(i ,jface-1,k ))
                f10 = fac1*(phi(in,jface,k ) - phi(in,jface-1,k ))*is
                f01 = fac1*(phi(i ,jface,kn) - phi(i ,jface-1,kn))*ks
                f11 = fac1*(phi(in,jface,kn) - phi(in,jface-1,kn))*is*ks

                c0 = ABS(fd1(i,jface,k,2))
                c1 = ABS(fd1(i,jface,k,4))

                Fy(jj) = (f00*(1.d0-c0)*(1.d0-c1) + f10*c0*(1.d0-c1) + f01*(1.d0-c0)*c1 + f11*c0*c1)*fd1(i,jface,k,1)
             enddo

             do kk=0,1
                kface = k+kk

                is = SIGN(1.d0, fd2(i,j,kface,2))
                js = SIGN(1.d0, fd2(i,j,kface,3))

                in = i + is
                jn = j + js

                f00 = fac2*(phi(i ,j ,kface) - phi(i ,j ,kface-1))
                f10 = fac2*(phi(in,j ,kface) - phi(in,j ,kface-1))*is
                f01 = fac2*(phi(i ,jn,kface) - phi(i ,jn,kface-1))*js
                f11 = fac2*(phi(in,jn,kface) - phi(in,jn,kface-1))*is*js

                c0 = ABS(fd2(i,j,kface,2))
                c1 = ABS(fd2(i,j,kface,3))

                Fz(jj) = (f00*(1.d0-c0)*(1.d0-c1) + f10*c0*(1.d0-c1) + f01*(1.d0-c0)*c1 + f11*c0*c1)*fd2(i,j,kface,1)
             enddo

             lph(i,j,k) = Fx(1) - Fx(0) + Fy(1) - Fy(0) + Fz(1) - Fz(0)

          enddo
       enddo
    enddo

  end subroutine lapleb_MSD

  subroutine lapl_MSD(&
       lph,lph_l1,lph_l2,lph_l3,lph_h1,lph_h2,lph_h3, &
       phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
       lo, hi, dx) bind(C,name="lapl_MSD")

    integer, intent(in) :: lph_l1,lph_l2,lph_l3,lph_h1,lph_h2,lph_h3
    integer, intent(in) :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3

    integer, intent(in) :: lo(3),hi(3)
                                                         
    double precision, intent(inout) :: lph(lph_l1:lph_h1,lph_l2:lph_h2,lph_l3:lph_h3)
    double precision, intent(in   ) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
    double precision, intent(in   ) :: dx(3)

    integer :: i,j,k
    double precision :: fac1, fac2, fac3

    fac1 = 1.d0 / (dx(1)*dx(1))
    fac2 = 1.d0 / (dx(2)*dx(2))
    fac3 = 1.d0 / (dx(3)*dx(3))

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             lph(i,j,k) = fac1*(phi(i+1,j,k) + phi(i-1,j,k) - 2*phi(i,j,k)) &
                  +       fac2*(phi(i,j+1,k) + phi(i,j-1,k) - 2*phi(i,j,k)) &
                  +       fac3*(phi(i,j,k+1) + phi(i,j,k-1) - 2*phi(i,j,k))
          enddo
       enddo
    enddo

  end subroutine lapl_MSD

end module stencil_test_module
