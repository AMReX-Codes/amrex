module bcoef
  use amrex_fort_module, only : amrex_real
private :: f 
public :: build_b_2d, build_b_3d
  interface f 
    module procedure f2
    module procedure f3
  end interface 
contains

  subroutine build_b_2d(lo, hi, problo, probhi, &
                        bx, bxlo, bxhi, & 
                        by, bylo, byhi, & 
                        dx) bind(c, name='build_b_2d')
  implicit none
  integer, dimension(2), intent(in   ) :: lo, hi, bxlo, bxhi, bylo, byhi 
  real(amrex_real), dimension(2), intent(in  ) :: problo, probhi, dx 
  real(amrex_real), intent(inout) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2))
  real(amrex_real), intent(inout) :: by(bylo(1):byhi(1),bylo(2):byhi(2))

  integer :: i, j
  real(amrex_real) :: xc, yc, xf, yf, xl, yl, xh, yh 
  xl = problo(1)
  yl = problo(2)
  xh = probhi(1)
  yh = probhi(2)

    do j = lo(2)+1, hi(2)-1
      yf = yl + dble(j)*dx(2)
      yc = yl + (dble(j) + 0.5d0)*dx(2)
      do i = lo(1)+1, hi(1)-1
    !interior
        xf = xl + dble(i)*dx(1)
        xc = xl + (dble(i) + 0.5d0)*dx(1)
        bx(i,j) = f(xf, yc) 
        by(i,j) = f(xc, yf)
      enddo 
    !xfaces 
      xf = xl + dble(hi(1)+1)*dx(1)
      bx(lo(1),j)   = f(xl,yc)
      by(lo(1),j)   = f(xl,yf)
      bx(hi(1),j)   = f(xf,yc)
      by(hi(1),j)   = f(xh,yf)
      bx(hi(1)+1,j) = f(xh,yc)
    enddo

    yf = yl + dble(hi(2)+1)*dx(2)
    do i = lo(1)+1, hi(1)-1
    !yfaces
      xc = xl + (dble(i) + 0.5d0)*dx(1)
      xf = xl + dble(i)*dx(1)
      bx(i,lo(2))   = f(xf,yl)
      by(i,lo(2))   = f(xc,yl)
      bx(i,hi(2))   = f(xf,yh)
      by(i,hi(2))   = f(xc,yf)
      by(i,hi(2)+1) = f(xc,yh)
    enddo
    xf = xl + dble(hi(1))*dx(1)
    yf = yl + dble(hi(2))*dx(2)
  !corners
    bx(lo(1),lo(2))   = f(xl,yl)
    by(lo(1),lo(2))   = f(xl,yl)

    bx(lo(1),hi(2))   = f(xl,yh)
    by(lo(1),hi(2))   = f(xl,yf)
    by(lo(1),hi(2)+1) = f(xl,yh)

    bx(hi(1),lo(2))   = f(xf,yl)
    bx(hi(1)+1,lo(2)) = f(xh,yl)
    by(hi(1),lo(2))   = f(xh,yl)

    bx(hi(1),hi(2))   = f(xf,yh)
    by(hi(1),hi(2))   = f(xh,yf)
    bx(hi(1)+1,hi(2)) = f(xh,yh)
    by(hi(1),hi(2)+1) = f(xh,yh)


  end subroutine build_b_2d

  subroutine build_b_3d(lo, hi, problo, probhi, &
                        bx, bxlo, bxhi, & 
                        by, bylo, byhi, & 
                        bz, bzlo, bzhi, &
                        dx) bind(c, name='build_b_3d')
  implicit none 
  integer, dimension(3), intent(in) :: lo, hi, bxlo, bxhi, bylo, byhi, bzlo, bzhi 
  real(amrex_real), dimension(3), intent(in) :: problo, probhi, dx 
  real(amrex_real), intent(inout) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
  real(amrex_real), intent(inout) :: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
  real(amrex_real), intent(inout) :: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))

  integer :: i, j, k 
  real(amrex_real) :: xc, xf, xl, xh, yc, yf, yl, yh, zc, zf, zl, zh
    xl = problo(1)
    yl = problo(2)
    zl = problo(3)
    xh = probhi(1)
    yh = probhi(2)
    zh = probhi(3) 
    do  k = lo(3)+1,hi(3)-1
      zf = zl + dble(k)*dx(3)
      zc = zl + (dble(k) + 0.5d0)*dx(3) 
      do j = lo(2)+1,hi(2)-1
        yf = yl + dble(j)*dx(2)
        yc = yl + (dble(j) + 0.5d0)*dx(2)
        do i = lo(1)+1,hi(1)-1
! interior
          xf = xl + dble(i)*dx(1) 
          xc = xl + (dble(i) + 0.5d0)*dx(1)
          bx(i,j,k) = f(xf,yc,zc)
          by(i,j,k) = f(xc,yf,zc) 
          bz(i,j,k) = f(xc,yc,zf)
        enddo
        xf = xl + dble(hi(1))*dx(1)
!    x faces
        bx(lo(1),j,k) = f(xl,yc,zc)
        by(lo(1),j,k) = f(xl,yf,zc)
        bz(lo(1),j,k) = f(xl,yc,zf)

        bx(hi(1),j,k) = f(xf,yc,zc)
        by(hi(1),j,k) = f(xh,yf,zc)
        bz(hi(1),j,k) = f(xh,yc,zf)
        bx(hi(1)+1,j,k) = f(xh,yc,zc)
      enddo
      yf = yl + dble(hi(2))*dx(2)
      do i = lo(1)+1,hi(1)-1
        xc = xl + (dble(i) + 0.5d0)*dx(1)
        xf = xl + (dble(i))*dx(1)
!   y faces
        bx(i,lo(2),k)   = f(xf,yl,zc)
        by(i,lo(2),k)   = f(xc,yl,zc)
        bz(i,lo(2),k)   = f(xc,yl,zf)

        bx(i,hi(2),k)   = f(xf,yh,zc)
        by(i,hi(2),k)   = f(xc,yf,zc)
        bz(i,hi(2),k)   = f(xc,yh,zf)
        by(i,hi(2)+1,k) = f(xc,yh,zc) 
      enddo

!    xy edges 

      bx(lo(1),lo(2),k)   = f(xl,yl,zc)
      by(lo(1),lo(2),k)   = f(xl,yl,zc)
      bz(lo(1),lo(2),k)   = f(xl,yl,zc)

      bx(lo(1),hi(2),k)   = f(xl,yl,zc)
      by(lo(1),hi(2),k)   = f(xl,yf,zc)
      bz(lo(1),hi(2),k)   = f(xl,yh,zf)
      by(lo(1),hi(2)+1,k) = f(xl,yh,zc)

      bx(hi(1),lo(2),k)   = f(xf,yl,zc)
      bx(hi(1)+1,lo(2),k) = f(xh,yl,zc)
      by(hi(1),lo(2),k)   = f(xh,yl,zc)
      bz(hi(1),lo(2),k)   = f(xh,yl,zc)

      bx(hi(1),hi(2),k)   = f(xf,yh,zc)
      bx(hi(1)+1,hi(2),k) = f(xh,yh,zc)
      by(hi(1),hi(2),k)   = f(xh,yf,zc)
      bz(hi(1),hi(2),k)   = f(xh,yh,zf)
      by(hi(1),hi(2)+1,k) = f(xh,yh,zc)
    enddo

    zf = yl + dble(hi(3))*dx(3)
    do j = lo(2)+1,hi(2)-1
      yc = yl + (dble(j) + 0.5d0)*dx(2)
      yf = yl + dble(j)*dx(2)
      do i = lo(1)+1,hi(1)-1
! z face
        xc = xl + (dble(i) + 0.5d0)*dx(1) 
        xf = xl + dble(i)*dx(1)
        bx(i,j,lo(3)) = f(xf,yc,zl)
        by(i,j,lo(3)) = f(xc,yf,zl)
        bz(i,j,lo(3)) = f(xc,yc,zl)

        bx(i,j,hi(3)) = f(xf,yc,zh)
        by(i,j,hi(3)) = f(xc,yf,zh)
        bz(i,j,hi(3)) = f(xc,yc,zf)
        bz(i,j,hi(3)+1) = f(xc,yc,zh)
      enddo
! xz edges 
      xf = xl + (dble(hi(1)))*dx(1)
      bx(lo(1),j,lo(3))   = f(xl,yc,zl)
      by(lo(1),j,lo(3))   = f(xl,yf,zl)
      bz(lo(1),j,lo(3))   = f(xl,yc,zl)

      bx(lo(1),j,hi(3))   = f(xl,yc,zh)
      by(lo(1),j,hi(3))   = f(xl,yf,zh)
      bz(lo(1),j,hi(3))   = f(xl,yc,zf)
      bz(lo(1),j,hi(3)+1) = f(xl,yc,zh)

      bx(hi(1),j,lo(3))   = f(xf,yc,zl)
      bx(hi(1)+1,j,lo(3)) = f(xh,yc,zl)
      by(hi(1),j,lo(3))   = f(xh,yf,zl)
      bz(hi(1),j,lo(3))   = f(xh,yc,zl)

      bx(hi(1),j,hi(3))   = f(xf,yc,zh)
      bx(hi(1)+1,j,hi(3)) = f(xh,yc,zh)
      by(hi(1),j,hi(3))   = f(xh,yf,zh)
      bz(hi(1),j,hi(3))   = f(xh,yc,zf)
      bz(hi(1),j,hi(3)+1) = f(xh,yc,zh)
    enddo

!   yz edges
    do i = lo(1)+1,hi(1)-1
      xc = xl +(dble(i) + 0.5)*dx(1) 
      xf = xl +dble(i)*dx(1)
      bx(i,lo(2),lo(3))   = f(xf,yl,zl)
      by(i,lo(2),lo(3))   = f(xc,yl,zl)
      bz(i,lo(2),lo(3))   = f(xc,yl,zl)

      bx(i,lo(2),hi(3))   = f(xf,yl,zh)
      by(i,lo(2),hi(3))   = f(xc,yl,zh)
      bz(i,lo(2),hi(3))   = f(xc,yl,zf)
      bz(i,lo(2),hi(3)+1) = f(xc,yl,zh)

      bx(i,hi(2),lo(3))   = f(xf,yh,zl)
      by(i,hi(2),lo(3))   = f(xc,yf,zl)
      by(i,hi(2)+1,lo(3)) = f(xc,yh,zl)
      bz(i,hi(2),lo(3))   = f(xc,yh,zl)

      bx(i,hi(2),hi(3))   = f(xf,yh,zh)
      by(i,hi(2),hi(3))   = f(xc,yf,zh)
      by(i,hi(2)+1,hi(3)) = f(xh,yh,zh)
      bz(i,hi(2),hi(3))   = f(xc,yh,zf)
      bz(i,hi(2),hi(3)+1) = f(xc,yh,zh)
    enddo

!     Cornerns...
      xf = xl + dble(hi(1))*dx(1)
      yf = yl + dble(hi(2))*dx(2)
      zf = zl + dble(hi(3))*dx(3)

      bx(lo(1),lo(2),lo(3))   = f(xl,yl,zl)
      by(lo(1),lo(2),lo(3))   = f(xl,yl,zl)
      bz(lo(1),lo(2),lo(3))   = f(xl,yl,zl)

      bx(lo(1),lo(2),hi(3))   = f(xl,yl,zh)
      by(lo(1),lo(2),hi(3))   = f(xl,yl,zh)
      bz(lo(1),lo(2),hi(3))   = f(xl,yl,zf)
      bz(lo(1),lo(2),hi(3)+1) = f(xl,yl,zh)

      bx(lo(1),hi(2),lo(3))   = f(xl,yh,zl)
      by(lo(1),hi(2),lo(3))   = f(xl,yf,zl)
      by(lo(1),hi(2)+1,lo(3)) = f(xl,yh,zl)
      bz(lo(1),hi(2),lo(3))   = f(xl,yh,zl)

      bx(lo(1),hi(2),hi(3))   = f(xl,yh,zh)
      by(lo(1),hi(2),hi(3))   = f(xl,yf,zh)
      by(lo(1),hi(2)+1,hi(3)) = f(xl,yh,zh)
      bz(lo(1),hi(2),hi(3))   = f(xl,yh,zf)
      bz(lo(1),hi(2),hi(3)+1) = f(xl,yh,zh)
     
      bx(hi(1),lo(2),lo(3))   = f(xf,yl,zl)
      bx(hi(1)+1,lo(2),lo(3)) = f(xh,yl,zl)
      by(hi(1),lo(2),lo(3))   = f(xh,yl,zl)
      bz(hi(1),lo(2),lo(3))   = f(xh,yl,zl)

      bx(hi(1),lo(2),hi(3))   = f(xf,yl,zh)
      bx(hi(1)+1,lo(2),hi(3)) = f(xh,yl,zh)
      by(hi(1),lo(2),hi(3))   = f(xh,yl,zh)
      bz(hi(1),lo(2),hi(3))   = f(xh,yl,zf)
      bz(hi(1),lo(2),hi(3)+1) = f(xh,yl,zh)

      bx(hi(1),hi(2),lo(3))   = f(xf,yh,zl)
      bx(hi(1)+1,hi(2),lo(3)) = f(xh,yh,zl)
      by(hi(1),hi(2),lo(3))   = f(xh,yf,zl)
      by(hi(1),hi(2)+1,lo(3)) = f(xh,yh,zl)
      bz(hi(1),hi(2),lo(3))   = f(xh,yh,zl)

      bx(hi(1),hi(2),hi(3))   = f(xf,yh,zh)
      bx(hi(1)+1,hi(2),hi(3)) = f(xh,yh,zh)
      by(hi(1),hi(2),hi(3))   = f(xh,yf,zh)
      by(hi(1),hi(2)+1,hi(3)) = f(xh,yh,zh)
      bz(hi(1),hi(2),hi(3))   = f(xh,yh,zf)
      bz(hi(1),hi(2),hi(3)+1) = f(xh,yh,zh)

  end subroutine build_b_3d

  function f2(x,y)
  real(amrex_real), intent(in) :: x, y 
  real(amrex_real) :: f2 
    f2 = sin(x)*cos(y)  
  end function f2 
  
  function f3(x,y,z)
  real(amrex_real), intent(in) :: x, y, z
  real(amrex_real) :: f3
    f3 = sin(x)*sin(y)*cos(z)
  end function f3

end module bcoef
