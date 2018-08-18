module acoef
  use amrex_fort_module, only : amrex_real
private :: f 
public :: build_a_2d, build_a_3d
  interface f 
    module procedure f2
    module procedure f3
  end interface 
contains

  subroutine build_a_2d(lo, hi, problo, probhi, &
                        a, alo, ahi, & 
                        dx) bind(c, name='build_a_2d')
  implicit none
  integer, dimension(2), intent(in   ) :: lo, hi, alo, ahi
  real(amrex_real), dimension(2), intent(in  ) :: problo, probhi, dx 
  real(amrex_real), intent(inout) :: a(alo(1):ahi(1),alo(2):ahi(2))

  integer :: i, j
  real(amrex_real) :: x, y, xl, yl, xh, yh 
    xl = problo(1)
    yl = problo(2)
    xh = probhi(1)
    yh = probhi(2)
    
    do j = lo(2)+1, hi(2)-1
      y = yl + (dble(j) + 0.5d0)*dx(2)
      do i = lo(1)+1, hi(1)-1
       x = xl + (dble(i) + 0.5d0)*dx(1) 
       a(i,j) = f(x,y)
      enddo
      a(lo(1),j) = f(xl,y)
      a(hi(1),j) = f(xh,y)
    enddo
    do i = lo(1)+1,hi(1)-1
      x = xl + (dble(i) + 0.5d0)*dx(1)
      a(i,lo(2)) = f(x,yl)
      a(i,hi(2)) = f(x,yh)
    enddo
    a(lo(1),lo(2)) = f(xl,yl)
    a(hi(1),lo(2)) = f(xh,yl)
    a(lo(1),hi(2)) = f(xl,yh)
    a(hi(1),hi(2)) = f(xh,yh)

  end subroutine build_a_2d

  subroutine build_a_3d(lo, hi, problo, probhi, &
                        a, alo, ahi, & 
                        dx) bind(c, name='build_a_3d')
  implicit none 
  integer, dimension(3), intent(in) :: lo, hi, alo, ahi
  real(amrex_real), dimension(3), intent(in) :: problo, probhi, dx 
  real(amrex_real), intent(inout) :: a(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

  integer :: i, j, k 
  real(amrex_real) :: x, xl, xh, y, yl, yh, z, zl, zh
    xl = problo(1)
    yl = problo(2)
    zl = problo(3)
    xh = probhi(1)
    yh = probhi(2)
    zh = probhi(3) 
    do k = lo(3)+1, hi(3)-1
      z  = zl + (dble(k) + 0.5d0)*dx(3)
      do j = lo(2)+1, hi(2)-1
        y  = yl + (dble(j) + 0.5d0)*dx(2)
        do i = lo(1)+1, hi(1)-1
          x  = xl + (dble(i) + 0.5d0)*dx(1)
          a(i,j,k) = f(x,y,z)
        enddo
        a(lo(1),j,k) = f(xl,y,z)
        a(hi(1),j,k) = f(xh,y,z)
      enddo
      do i = lo(1)+1,hi(1)-1
        x = xl + (dble(i) + 0.5d0)*dx(1)
        a(i,lo(2),k) = f(x,yl,z)
        a(i,hi(2),k) = f(x,yh,z)
      enddo
      a(lo(1),lo(2),k) = f(xl,yl,z)
      a(hi(1),lo(2),k) = f(xh,yl,z)
      a(lo(1),hi(2),k) = f(xl,yh,z)
      a(hi(1),hi(2),k) = f(xh,yh,z)
    enddo
    do j = lo(2)+1, hi(2)-1
      y = yl + (dble(j) + 0.5d0)*dx(2)
      do i = lo(1)+1, hi(1)-1
        x = xl + (dble(i) + 0.5d0)*dx(1)
        a(i,j,lo(3)) = f(x,y,zl)
        a(i,j,hi(3)) = f(x,y,zh)
      enddo
      a(lo(1),j,lo(3)) = f(xl,y,zl) 
      a(hi(1),j,lo(3)) = f(xh,y,zl) 
      a(lo(1),j,hi(3)) = f(xl,y,zh)
      a(hi(1),j,hi(3)) = f(xh,y,zh)
    enddo
    a(lo(1),lo(2),lo(3)) = f(xl,yl,zl)
    a(hi(1),lo(2),lo(3)) = f(xh,yl,zl)
    a(lo(1),hi(2),lo(3)) = f(xl,yh,zl)
    a(hi(1),hi(2),lo(3)) = f(xh,yh,zl)
    a(lo(1),lo(2),hi(3)) = f(xl,yl,zh)
    a(hi(1),lo(2),hi(3)) = f(xh,yl,zh)
    a(lo(1),hi(2),hi(3)) = f(xh,yh,zh)
    a(hi(1),hi(2),hi(3)) = f(xh,yh,zh)

  end subroutine build_a_3d

  function f2(x,y)
  real(amrex_real), intent(in) :: x, y 
  real(amrex_real) :: f2 
    f2 = cos(x)*cos(y)  
  end function f2 
  
  function f3(x,y,z)
  real(amrex_real), intent(in) :: x, y, z
  real(amrex_real) :: f3
    f3 = cos(x)*sin(y)*cos(z)
  end function f3

end module acoef
