
module stencil_test_module

  implicit none

  private

  public :: lapleb_MSD, lapl_MSD, init_phi

contains

  subroutine lapleb_MSD(&
       lph,lph_l1,lph_l2,lph_h1,lph_h2, &
       phi,phi_l1,phi_l2,phi_h1,phi_h2, &
       cd ,cd_l1, cd_l2, cd_h1, cd_h2, &
       fd0,fd0_l1,fd0_l2,fd0_h1,fd0_h2, &
       fd1,fd1_l1,fd1_l2,fd1_h1,fd1_h2, &
       lo, hi, dx) bind(C,name="lapleb_MSD")

    integer, intent(in) :: lph_l1,lph_l2,lph_h1,lph_h2
    integer, intent(in) :: phi_l1,phi_l2,phi_h1,phi_h2
    integer, intent(in) ::  cd_l1, cd_l2, cd_h1, cd_h2
    integer, intent(in) :: fd0_l1,fd0_l2,fd0_h1,fd0_h2
    integer, intent(in) :: fd1_l1,fd1_l2,fd1_h1,fd1_h2

    integer, intent(in) :: lo(2),hi(2)

    double precision, intent(inout) :: lph(lph_l1:lph_h1,lph_l2:lph_h2)
    double precision, intent(in   ) :: phi(phi_l1:phi_h1,phi_l2:phi_h2)
    double precision, intent(in   ) ::  cd( cd_l1:cd_h1 , cd_l2:cd_h2, 5)
    double precision, intent(in   ) :: fd0(fd0_l1:fd0_h1,fd0_l2:fd0_h2,3)
    double precision, intent(in   ) :: fd1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,3)
    double precision, intent(in   ) :: dx

    integer :: i,j,ip,im,jp,jm,in,jn,ii,jj,iface,jface,is,js
    double precision :: fac, c0, f0, f1, Fx(0:1), Fy(0:1)
    double precision :: n(2), x1, x2, y1, y2, v1, v2, xbc, ybc, vbc, dVdn, l1, l2
    double precision :: u(3), x, y, d(2), v(2)
    integer :: iv(2), c1, c2, L, M, s1, s2

    fac = 1.d0 / (dx*dx)

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          do ii=0,1
             iface = i+ii
             js = SIGN(1.d0, fd0(iface,j,3))
             jn = j + js
             f0 = fac*(phi(iface,j ) - phi(iface-1,j ))
             f1 = fac*(phi(iface,jn) - phi(iface-1,jn))
             c0 = ABS(fd0(iface,j,3))
             Fx(ii) = (f0*(1.d0 - c0) + f1*c0)*fd0(iface,j,1)
          enddo

          do jj=0,1
             jface = j+jj
             is = SIGN(1.d0, fd1(i,jface,2))
             in = i + is
             f0 = fac*(phi(i ,jface) - phi(i ,jface-1))
             f1 = fac*(phi(in,jface) - phi(in,jface-1))
             c0 = ABS(fd1(i,jface,2))
             Fy(jj) = (f0*(1.d0 - c0) + f1*c0)*fd1(i,jface,1)
          enddo

          lph(i,j) = Fx(1) - Fx(0) + Fy(1) - Fy(0)

       enddo
    enddo

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          n(1) = cd(i,j,2)
          n(2) = cd(i,j,3)

          if (ABS(n(1)).ne.0.d0 .and. ABS(n(2)).ne.0.d0) then

             c1 = 1 + INT( MIN(1.d0, ABS(n(2)/n(1)) ) ) ! Either 1 or 2 only
             c2 = 3 - c1 ! Either 1 or 2 only

             s1 = SIGN(1.d0,n(c1))
             s2 = SIGN(1.d0,n(c2))

             xbc = cd(i,j,3+c1) * s1
             ybc = cd(i,j,3+c2) * s2
             vbc = cd(i,j,1)
             
             do M=1,2
                iv(c1) = M * s1
                do L=1,3
                   iv(c2) = (L - 1)*s2
                   u(L) = phi(iv(1)+i,iv(2)+j)
                enddo
                x = M
                y = ybc + (x - xbc)*ABS(cd(i,j,1+c2)/cd(i,j,1+c1))

                v(M) =   0.5d0 * (y-1.d0) * (y-2.d0) * u(1) &
                     -               y    * (y-2.d0) * u(2) &
                     +      0.5d0 *  y    * (y-1.d0) * u(3)

                d(M) = SQRT( (x-xbc)**2 + (y-ybc)**2)
             enddo

             dVdn = - vbc  * (d(2)+d(1))/(d(2)*d(1)) &
                  +   v(1) * d(2)/(d(1)*(d(2)-d(1))) &
                  -   v(2) * d(1)/(d(2)*(d(2)-d(1)))
          endif

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
    double precision, intent(in   ) :: dx

    integer :: i,j
    double precision :: fac

    fac = 1.d0 / (dx*dx)

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          lph(i,j) = fac*(phi(i+1,j) + phi(i-1,j) - 2*phi(i,j)) &
               +     fac*(phi(i,j+1) + phi(i,j-1) - 2*phi(i,j))
       enddo
    enddo

  end subroutine lapl_MSD

  subroutine init_phi(&
       phi,phi_l1,phi_l2,phi_h1,phi_h2, &
       lo, hi, plo, dx) bind(C,name="init_phi")

    integer, intent(in) :: phi_l1,phi_l2,phi_h1,phi_h2

    integer, intent(in) :: lo(2),hi(2)

    double precision, intent(inout) :: phi(phi_l1:phi_h1,phi_l2:phi_h2)
    double precision, intent(in   ) :: plo(2), dx

    integer :: i,j
    double precision :: x,y

    do j=lo(2),hi(2)
       y = plo(2) + (j+0.5d0)*dx
       do i=lo(1),hi(1)
          x = plo(1) + (i+0.5d0)*dx
          phi(i,j) = 0.5d0 * (x-0.5d0)**2
       enddo
    enddo

  end subroutine init_phi

end module stencil_test_module
