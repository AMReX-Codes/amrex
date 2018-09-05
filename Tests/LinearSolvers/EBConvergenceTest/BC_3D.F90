module bc
use amrex_fort_module, only: amrex_real

private :: f3
public :: apply_bc

contains
  subroutine apply_bc(lo, hi, phi, phlo, phhi, &
                      dx, problo, probhi)&
                      bind(c, name='apply_bc')
  implicit none
  integer, dimension(3), intent(in) :: lo, hi, phlo, phhi
    real(amrex_real), dimension(3), intent(in) :: dx, problo, probhi
    real(amrex_real), intent(inout) :: phi(phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3))
  
    real(amrex_real) :: xl, xh, x, yl, yh, y, zl, zh, z, denom
    integer          :: i, j, k

    xl = problo(1) - 0.5d0 
    xh = probhi(1) - 0.5d0
    yl = problo(2) - 0.5d0
    yh = probhi(2) - 0.5d0
    zl = problo(3) - 0.5d0
    zh = probhi(3) - 0.5d0 

    do k = lo(3)+1, hi(3)-1
      z = zl + (dble(k) + 0.5d0)*dx(3)
      do j = lo(2)+1, hi(2)-1
      !xface
        y = yl + (dble(j) + 0.5d0)*dx(2)
        phi(lo(1),j,k) = f3(xl,y,z)
        phi(hi(1),j,k) = f3(xh,y,z)
      end do
      do i = lo(1)+1, hi(1)-1
      !yface
        x = xl + (dble(i) + 0.5d0)*dx(1)
        phi(i,lo(2),k) = f3(x,yl,z)
        phi(i,hi(2),k) = f3(x,yh,z)
      end do
      !x-y edges 
      phi(lo(1),lo(2),k) = f3(xl,yl,z)
      phi(hi(1),hi(2),k) = f3(xh,yh,z)
      phi(lo(1),hi(2),k) = f3(xl,yh,z)
      phi(hi(1),lo(2),k) = f3(xh,yl,z)
    enddo
    
    !z-face 
    do j = lo(2)+1, hi(2)-1
      y = yl + (dble(j) + 0.5d0)*dx(2)
      do i = lo(1)+1,hi(2)-1
        x = xl +(dble(i) +0.5d0)*dx(1) 
        phi(i,j,lo(3)) = f3(x,y,zl)
        phi(i,j,hi(3)) = f3(x,y,zh)
      enddo
      !xz edge
      phi(lo(1),j,lo(3)) = f3(xl,y,zl)
      phi(lo(1),j,hi(3)) = f3(xl,y,zh) 
      phi(hi(1),j,hi(3)) = f3(xh,y,zh)
      phi(hi(1),j,lo(3)) = f3(xh,y,zl)
    enddo

    !yz edge
    do i = lo(1)+1,hi(1)-1
      x = xl +(dble(i) +0.5d0)*dx(1)
      phi(i,lo(2),lo(3)) = f3(x,yl,zl)
      phi(i,lo(2),hi(3)) = f3(x,yl,zh) 
      phi(i,hi(2),hi(3)) = f3(x,yh,zh)
      phi(i,hi(2),lo(3)) = f3(x,yh,zl)
    enddo

    !Corners! 
    phi(lo(1),lo(2),lo(3)) = f3(xl,yl,zl)
    phi(lo(1),lo(2),hi(3)) = f3(xl,yl,zh) 
    phi(lo(1),hi(2),hi(3)) = f3(xl,yh,zh)
    phi(lo(1),hi(2),lo(3)) = f3(xl,yh,zl)
    phi(hi(1),lo(2),lo(3)) = f3(xh,yl,zl)
    phi(hi(1),lo(2),hi(3)) = f3(xh,yl,zh) 
    phi(hi(1),hi(2),hi(3)) = f3(xh,yh,zh)
    phi(hi(1),hi(2),lo(3)) = f3(xh,yh,zl)

  end subroutine apply_bc
  
  function f3(x,y,z)
    real(amrex_real), intent(in) :: x, y, z
    real(amrex_real) :: f3
    
    f3 = (x+z)/sqrt(x*x + y*y + z*z)
  end function f3 

end module bc
