module rhs
 use amrex_fort_module, only : amrex_real
 
private :: RHS3, dbdx, dbdy, dbdz
public :: build_rhs_2D, build_rhs_3D
  interface dbdx 
    module procedure dbdx2
    module procedure dbdx3
  end interface dbdx
  interface dbdy
    module procedure dbdy2
    module procedure dbdy3
  end interface dbdy
  interface dbdz 
    module procedure dbdz3
  end interface dbdz

contains 
  subroutine build_rhs_2D(lo, hi, rhs, rlo, rhi, &
                               a  , alo, ahi, &
                             bx , bxlo, bxhi, &
                             by , bylo, byhi, & 
                          dx, problo, probhi) &
                          bind(c, name='build_rhs_2d')
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
        b     = 0.25d0*(bx(i+1,j) + bx(i,j) + by(i,j+1) + by(i,j))
        dxb   = (bx(i+1,j) - bx(i,j))*dxinv
        x     = xl + (dble(i) + 0.5d0)*dx(1) - 0.5d0; 
!        dxb   = dbdx(x,y)
!        dyb   = dbdy(x,y)
        denom = sqrt(x*x + y*y) 
        term1 = a(i,j)*x/denom
        term2 = (y*(y*dxb -x*dyb) - x*b)/(denom**3)
        rhs(i,j) = term1 - term2
      enddo
      !Boundary Case xlo and xhi, y inside
      b   = 0.25d0*(bx(lo(1)+1,j) + bx(lo(1),j) + by(lo(1),j+1) + by(lo(1),j))
      dyb = (by(lo(1), j+1) - by(lo(1),j))*(dyinv)
      x = xl - 0.5d0
      dxb = (bx(lo(1)+1,j)-bx(lo(1),j))*(dxinv)  
!      dxb   = dbdx(x,y)
!      dyb   = dbdy(x,y)
      denom = sqrt(x*x + y*y)
      term1 = a(lo(1),j)*x/denom
      term2 = (y*(y*dxb - x*dyb) - x*b)/(denom**3)
      rhs(lo(1),j) = term1 - term2 

      dyb = (by(hi(1), j+1) - by(hi(1), j))*dyinv
      x = xh - 0.5d0 
      b   = 0.25d0*(bx(hi(1)+1,j) + bx(hi(1),j) + by(hi(1),j+1) + by(hi(1),j))
      dxb = (bx(hi(1)+1, j) - bx(hi(1),j))*(dxinv*0.5)
!      dxb   = dbdx(x,y)
!      dyb   = dbdy(x,y)
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
!      dxb   = dbdx(x,y)
!      dyb   = dbdy(x,y)
      denom = sqrt(x**2 + y**2)
      term1 = a(i,lo(2))*x/denom
      term2 = (y*(y*dxb - x*dyb) - x*b)/(denom**3)
      rhs(i,lo(2)) = term1 - term2 

      y = yh -0.5d0
      b   = 0.25d0*(bx(i+1,hi(2)) + bx(i,hi(2)) + by(i,hi(2)+1) + by(i,hi(2)))
      dyb = (by(i,hi(2)+1) - by(i,hi(2)))*dyinv
      dxb = (bx(i+1,hi(2)) - bx(i,hi(2)))*dxinv
!      dxb   = dbdx(x,y)
!      dyb   = dbdy(x,y)
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
  
  end subroutine build_rhs_2D

  subroutine build_rhs_3D(lo, hi, rhs, rlo, rhi, &
                               a  , alo, ahi, &
                             bx , bxlo, bxhi, &
                             by , bylo, byhi, &
                             bz , bzlo, bzhi, &
                          dx, problo, probhi) &
                          bind(c, name='build_rhs_3d')
  implicit none
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, alo, ahi, bxlo, bxhi, bylo, byhi, bzlo, bzhi 
    real(amrex_real), dimension(3), intent(in) :: dx, problo, probhi
    real(amrex_real), intent(inout) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) :: a  ( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3))
    real(amrex_real), intent(in   ) :: bx (bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
    real(amrex_real), intent(in   ) :: by (bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
    real(amrex_real), intent(in   ) :: bz (bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
   
    real(amrex_real) :: dxb, dyb, dzb, denom, x, y, z, dxinv, dyinv, dzinv
    real(amrex_real) :: term1, term2, b, xl, xh, yl, yh, zl, zh
    integer :: i, j, k
   
    dxinv   = 1.d0/dx(1)
    dyinv   = 1.d0/dx(2)
    dzinv   = 1.d0/dx(3)
    zl = problo(3)
    yl = problo(2)
    xl = problo(1)
    zh = probhi(3)
    yh = probhi(2)
    xh = probhi(1) 
    
    do k = lo(3)+1, hi(3)-1
      z = zl + (dble(k) + 0.5d0)*dx(3) - 0.5 ! Shifted to center of sphere <0.5 0.5 0.5>
      do j = lo(2)+1, hi(2)-1
        y = yl + (dble(j) + 0.5d0)*dx(2) - 0.5
        do i = lo(1)+1, hi(1)-1 !Interior
          b = 1.d0/6.d0*(bx(i+1,j,k) + bx(i,j,k) + by(i,j+1,k) + by(i,j,k) + bz(i,j,k+1) + bz(i,j,k)) !Average to CC
          x = xl + (dble(i) + 0.5d0)*dx(1) - 0.5 
          dzb = (bz(i,j,k+1) - bz(i,j,k))*dzinv
          dyb = (by(i,j+1,k) - by(i,j,k))*dyinv
          dxb = (bx(i+1,j,k) - bx(i,j,k))*dxinv
          rhs(i,j,k) = RHS3(a(i,j,k), b, x, y, z, dxb, dyb, dzb)
        enddo
        !Xlo faces 
        x = xl - 0.5d0
        b = bx(lo(1), j, k)
        dxb = (bx(lo(1)+1,j,k) - b)*dxinv 
        dyb = (by(lo(1),j+1,k) - by(lo(1),j,k))*dyinv
        dzb = (bz(lo(1),j,k+1) - bz(lo(1),j,k))*dzinv 
        rhs(lo(1), j, k) = RHS3(a(lo(1),j,k), b, x, y, z, dxb, dyb, dzb)

        !Xhi faces 
        x = xh - 0.5d0
        b = (bx(hi(1),j,k) + bx(hi(1)+1,j,k))*0.5d0
        dxb = (bx(hi(1)+1,j,k) - b)*dxinv 
        dyb = (by(hi(1),j+1,k) - by(hi(1),j,k))*dyinv
        dzb = (bz(hi(1),j,k+1) - bz(hi(1),j,k))*dzinv 
        rhs(hi(1), j, k) = RHS3(a(hi(1),j,k), b, x, y, z, dxb, dyb, dzb)
      enddo

      do i = lo(1)+1, hi(1)-1
        x = xl + (dble(i) + 0.5d0)*dx(1) - 0.5d0
      !Ylo Faces
        y = yl - 0.5d0 
        b = by(i,lo(2),k)
        dxb = (bx(i+1,lo(2),k) - bx(i,lo(2),k))*dxinv
        dyb = (by(i,lo(2)+1,k) - by(i,lo(2),k))*dyinv
        dzb = (bz(i,lo(2),k+1) - bz(i,lo(2),k))*dzinv
        rhs(i,lo(2),k) = RHS3(a(i,lo(2),k), b, x, y, z, dxb, dyb, dzb) 


      !Yhi Faces
        y = yh - 0.5d0 
        b = (by(i,hi(2),k) + by(i,hi(2)+1,k))*0.5d0
        dxb = (bx(i+1,hi(2),k) - bx(i,hi(2),k))*dxinv
        dyb = (by(i,hi(2)+1,k) - by(i,hi(2),k))*dyinv
        dzb = (bz(i,hi(2),k+1) - bz(i,hi(2),k))*dzinv
        rhs(i,hi(2),k) = RHS3(a(i,hi(2),k), b, x, y, z, dxb, dyb, dzb) 
     enddo

     !xlo ylo edge
      y = yl - 0.5d0
      x = xl - 0.5d0
      b = 0.5d0*(bx(lo(1), lo(2), k) + by(lo(1),lo(2),k))
      dxb = (bx(lo(1)+1,lo(2),k) - bx(lo(1),lo(2),k))*dxinv 
      dyb = (by(lo(1),lo(2)+1,k) - by(lo(1),lo(2),k))*dyinv
      dzb = (bz(lo(1),lo(2),k+1) - bz(lo(1),lo(2),k))*dzinv 
      rhs(lo(1), lo(2), k) = RHS3(a(lo(1),lo(2),k), b, x, y, z, dxb, dyb, dzb)

      !xhi ylo edge
      x = xh - 0.5d0 
      b = 0.5d0*(bx(hi(1), lo(2), k) + by(hi(1),lo(2),k))
      dxb = (bx(hi(1)+1,lo(2),k) - bx(hi(1),lo(2),k))*dxinv 
      dyb = (by(hi(1),lo(2)+1,k) - by(hi(1),lo(2),k))*dyinv
      dzb = (bz(hi(1),lo(2),k+1) - bz(hi(1),lo(2),k))*dzinv 
      rhs(hi(1), lo(2), k) = RHS3(a(hi(1),lo(2),k), b, x, y, z, dxb, dyb, dzb)
      
      !xhi yhi edge 
      y = yh - 0.5d0
      b = 0.5d0*(bx(hi(1), hi(2), k) + by(hi(1),hi(2),k))
      dxb = (bx(hi(1)+1,hi(2),k) - bx(hi(1),hi(2),k))*dxinv 
      dyb = (by(hi(1),hi(2)+1,k) - by(hi(1),hi(2),k))*dyinv
      dzb = (bz(hi(1),hi(2),k+1) - bz(hi(1),hi(2),k))*dzinv 
      rhs(hi(1), hi(2), k) = RHS3(a(hi(1),hi(2),k), b, x, y, z, dxb, dyb, dzb)

      !xlo yhi edge 
      x = xl -0.5d0
      b = 0.5d0*(bx(lo(1), hi(2), k) + by(lo(1),hi(2),k))
      dxb = (bx(lo(1)+1,hi(2),k) - bx(lo(1),hi(2),k))*dxinv 
      dyb = (by(lo(1),hi(2)+1,k) - by(lo(1),hi(2),k))*dyinv
      dzb = (bz(lo(1),hi(2),k+1) - bz(lo(1),hi(2),k))*dzinv 
      rhs(lo(1), hi(2), k) = RHS3(a(lo(1),hi(2),k), b, x, y, z, dxb, dyb, dzb)
    enddo 

    
    do j = lo(2)+1,hi(2)-1
      y = yl + (dble(j) + 0.5d0)*dx(2) - 0.5d0
      do i = lo(1)+1,hi(1)-1
        x = xl + (dble(i) + 0.5d0)*dx(1) -0.5d0
        !zlo face 
        z = zl - 0.5d0
        b = bz(i,j,lo(3))
        dxb = (bx(i+1,j,lo(3))-bx(i,j,lo(3)))*dxinv
        dyb = (by(i,j+1,lo(3))-by(i,j,lo(3)))*dyinv
        dzb = (bz(i,j,lo(3)+1)-b)*dzinv
        rhs(i,j,lo(3)) = RHS3(a(i,j,lo(3)),b,x,y,z,dxb, dyb, dzb)

        !zhi face 
        z = zh - 0.5d0
        b = (bz(i,j,hi(3)) + bz(i,j,hi(3)+1))*0.5d0 
        dxb = (bx(i+1,j,hi(3))-bx(i,j,hi(3)))*dxinv
        dyb = (by(i,j+1,hi(3))-by(i,j,hi(3)))*dyinv
        dzb = (bz(i,j,hi(3)+1)-b)*dzinv
        rhs(i,j,hi(3)) = RHS3(a(i,j,hi(3)),b,x,y,z,dxb, dyb, dzb)
      enddo
      !zlo xlo edge
      z = zl - 0.5d0
      x = xl - 0.5d0
      b = 0.5d0*(bx(lo(1), j, lo(3)) + bz(lo(1),j,lo(3)))
      dxb = (bx(lo(1)+1,j,lo(3)) - bx(lo(1),j,lo(3)))*dxinv 
      dyb = (by(lo(1),j+1,lo(3)) - by(lo(1),j,lo(3)))*dyinv
      dzb = (bz(lo(1),j,lo(3)+1) - bz(lo(1),j,lo(3)))*dzinv 
      rhs(lo(1), j, lo(3)) = RHS3(a(lo(1),j,lo(3)), b, x, y, z, dxb, dyb, dzb)

      !zlo xhi edge
      x = xh - 0.5d0 
      b = 0.5d0*(bx(hi(1), j, lo(3)) + bz(hi(1),j,lo(3)))
      dxb = (bx(hi(1)+1,j,lo(3)) - bx(hi(1),j,lo(3)))*dxinv 
      dyb = (by(hi(1),j+1,lo(3)) - by(hi(1),j,lo(3)))*dyinv
      dzb = (bz(hi(1),j,lo(3)+1) - bz(hi(1),j,lo(3)))*dzinv 
      rhs(hi(1), j, lo(3)) = RHS3(a(hi(1),j,lo(3)), b, x, y, z, dxb, dyb, dzb)
      
      !zhi xhi edge 
      z = zh - 0.5d0
      b = 0.5d0*(bx(hi(1), j, hi(3)) + bz(hi(1),j,hi(3)))
      dxb = (bx(hi(1)+1,j,hi(3)) - bx(hi(1),j,hi(3)))*dxinv 
      dyb = (by(hi(1),j+1,hi(3)) - by(hi(1),j,hi(3)))*dyinv
      dzb = (bz(hi(1),j,hi(3)+1) - bz(hi(1),j,hi(3)))*dzinv 
      rhs(hi(1), j, hi(3)) = RHS3(a(hi(1),j,hi(3)), b, x, y, z, dxb, dyb, dzb)

      !zhi xlo edge 
      x = xl -0.5d0
      b = 0.5d0*(bx(lo(1), j, hi(3)) + bz(lo(1),j,hi(3)))
      dxb = (bx(lo(1)+1,j,hi(3)) - bx(lo(1),j,hi(3)))*dxinv 
      dyb = (by(lo(1),j+1,hi(3)) - by(lo(1),j,hi(3)))*dyinv
      dzb = (bz(lo(1),j,hi(3)+1) - bz(lo(1),j,hi(3)))*dzinv 
      rhs(lo(1), j, hi(3)) = RHS3(a(lo(1),j,hi(3)), b, x, y, z, dxb, dyb, dzb)
   enddo 
  
   do i = lo(1)+1,hi(1)-1 
       !zlo ylo edge
      z = zl - 0.5d0
      y = yl - 0.5d0
      b = 0.5d0*(by(i, lo(2), lo(3)) + bz(i,lo(2),lo(3)))
      dxb = (bx(i+1,lo(2),lo(3)) - bx(i,lo(2),lo(3)))*dxinv 
      dyb = (by(i,lo(2)+1,lo(3)) - by(i,lo(2),lo(3)))*dyinv
      dzb = (bz(i,lo(2),lo(3)+1) - bz(i,lo(2),lo(3)))*dzinv 
      rhs(i, lo(2), lo(3)) = RHS3(a(i,lo(2),lo(3)), b, x, y, z, dxb, dyb, dzb)

      !zlo yhi edge
      y = yh - 0.5d0 
      b = 0.5d0*(by(i, hi(2), lo(3)) + bz(i,hi(2),lo(3)))
      dxb = (bx(i+1,hi(2),lo(3)) - bx(i,hi(2),lo(3)))*dxinv 
      dyb = (by(i,hi(2)+1,lo(3)) - by(i,hi(2),lo(3)))*dyinv
      dzb = (bz(i,hi(2),lo(3)+1) - bz(i,hi(2),lo(3)))*dzinv 
      rhs(i, hi(2), lo(3)) = RHS3(a(i,hi(2),lo(3)), b, x, y, z, dxb, dyb, dzb)
      
      !zhi yhi edge 
      z = zh - 0.5d0
      b = 0.5d0*(by(i, hi(2), hi(3)) + bz(i,hi(2),hi(3)))
      dxb = (bx(i+1,hi(2),hi(3)) - bx(i,hi(2),hi(3)))*dxinv 
      dyb = (by(i,hi(2)+1,hi(3)) - by(i,hi(2),hi(3)))*dyinv
      dzb = (bz(i,hi(2),hi(3)+1) - bz(i,hi(2),hi(3)))*dzinv 
      rhs(i, hi(2), hi(3)) = RHS3(a(i,hi(2),hi(3)), b, x, y, z, dxb, dyb, dzb)

      !zhi ylo edge 
      y = yl -0.5d0
      b = 0.5d0*(by(i, lo(2), hi(3)) + bz(i,lo(2),hi(3)))
      dxb = (bx(i+1,lo(2),hi(3)) - bx(i,lo(2),hi(3)))*dxinv 
      dyb = (by(i,lo(2)+1,hi(3)) - by(i,lo(2),hi(3)))*dyinv
      dzb = (bz(i,lo(2),hi(3)+1) - bz(i,lo(2),hi(3)))*dzinv 
      rhs(i, lo(2), hi(3)) = RHS3(a(i,lo(2),hi(3)), b, x, y, z, dxb, dyb, dzb)
   enddo 

  !Corners 
  !xlo ylo zlo 
  z = zl - 0.5d0
  y = yl - 0.5d0
  x = xl - 0.5d0 
  b = 1.d0/3.d0*(bx(lo(1),lo(2),lo(3)) + by(lo(1), lo(2), lo(3)) + bz(lo(1),lo(2),lo(3)))
  dxb = (bx(lo(1)+1,lo(2),lo(3)) - bx(lo(1),lo(2),lo(3)))*dxinv 
  dyb = (by(lo(1),lo(2)+1,lo(3)) - by(lo(1),lo(2),lo(3)))*dyinv
  dzb = (bz(lo(1),lo(2),lo(3)+1) - bz(lo(1),lo(2),lo(3)))*dzinv 
  rhs(lo(1), lo(2), lo(3)) = RHS3(a(lo(1),lo(2),lo(3)), b, x, y, z, dxb, dyb, dzb)

  !xlo ylo zhi 
  z = zh - 0.5d0
  b = 1.d0/3.d0*(bx(lo(1),lo(2),hi(3)) + by(lo(1), lo(2), hi(3)) + bz(lo(1),lo(2),hi(3)))
  dxb = (bx(lo(1)+1,lo(2),hi(3)) - bx(lo(1),lo(2),hi(3)))*dxinv 
  dyb = (by(lo(1),lo(2)+1,hi(3)) - by(lo(1),lo(2),hi(3)))*dyinv
  dzb = (bz(lo(1),lo(2),hi(3)+1) - bz(lo(1),lo(2),hi(3)))*dzinv 
  rhs(lo(1), lo(2), hi(3)) = RHS3(a(lo(1),lo(2),hi(3)), b, x, y, z, dxb, dyb, dzb)

  !xlo yhi zlo 
  z = zl - 0.5d0
  y = yh - 0.5d0
  b = 1.d0/3.d0*(bx(lo(1),hi(2),lo(3)) + by(lo(1), hi(2), lo(3)) + bz(lo(1),hi(2),lo(3)))
  dxb = (bx(lo(1)+1,hi(2),lo(3)) - bx(lo(1),hi(2),lo(3)))*dxinv 
  dyb = (by(lo(1),hi(2)+1,lo(3)) - by(lo(1),hi(2),lo(3)))*dyinv
  dzb = (bz(lo(1),hi(2),lo(3)+1) - bz(lo(1),hi(2),lo(3)))*dzinv 
  rhs(lo(1), hi(2), lo(3)) = RHS3(a(lo(1),hi(2),lo(3)), b, x, y, z, dxb, dyb, dzb)

  !xlo yhi zhi 
  z = zh - 0.5d0
  b = 1.d0/3.d0*(bx(lo(1),hi(2),hi(3)) + by(lo(1), hi(2), hi(3)) + bz(lo(1),hi(2),hi(3)))
  dxb = (bx(lo(1)+1,hi(2),hi(3)) - bx(lo(1),hi(2),hi(3)))*dxinv 
  dyb = (by(lo(1),hi(2)+1,hi(3)) - by(lo(1),hi(2),hi(3)))*dyinv
  dzb = (bz(lo(1),hi(2),hi(3)+1) - bz(lo(1),hi(2),hi(3)))*dzinv 
  rhs(lo(1), hi(2), hi(3)) = RHS3(a(lo(1),hi(2),hi(3)), b, x, y, z, dxb, dyb, dzb)

  !xhi ylo zlo 
  z = zl - 0.5d0
  y = yl - 0.5d0
  x = xh - 0.5d0 
  b = 1.d0/3.d0*(bx(hi(1),lo(2),lo(3)) + by(hi(1), lo(2), lo(3)) + bz(hi(1),lo(2),lo(3)))
  dxb = (bx(hi(1)+1,lo(2),lo(3)) - bx(hi(1),lo(2),lo(3)))*dxinv 
  dyb = (by(hi(1),lo(2)+1,lo(3)) - by(hi(1),lo(2),lo(3)))*dyinv
  dzb = (bz(hi(1),lo(2),lo(3)+1) - bz(hi(1),lo(2),lo(3)))*dzinv 
  rhs(hi(1), lo(2), lo(3)) = RHS3(a(hi(1),lo(2),lo(3)), b, x, y, z, dxb, dyb, dzb)

  !xhi ylo zhi 
  z = zh - 0.5d0
  b = 1.d0/3.d0*(bx(hi(1),lo(2),hi(3)) + by(hi(1), lo(2), hi(3)) + bz(hi(1),lo(2),hi(3)))
  dxb = (bx(hi(1)+1,lo(2),hi(3)) - bx(hi(1),lo(2),hi(3)))*dxinv 
  dyb = (by(hi(1),lo(2)+1,hi(3)) - by(hi(1),lo(2),hi(3)))*dyinv
  dzb = (bz(hi(1),lo(2),hi(3)+1) - bz(hi(1),lo(2),hi(3)))*dzinv 
  rhs(hi(1), lo(2), hi(3)) = RHS3(a(hi(1),lo(2),hi(3)), b, x, y, z, dxb, dyb, dzb)

  !xhi yhi zlo 
  z = zl - 0.5d0
  y = yh - 0.5d0
  b = 1.d0/3.d0*(bx(hi(1),hi(2),lo(3)) + by(hi(1), hi(2), lo(3)) + bz(hi(1),hi(2),lo(3)))
  dxb = (bx(hi(1)+1,hi(2),lo(3)) - bx(hi(1),hi(2),lo(3)))*dxinv 
  dyb = (by(hi(1),hi(2)+1,lo(3)) - by(hi(1),hi(2),lo(3)))*dyinv
  dzb = (bz(hi(1),hi(2),lo(3)+1) - bz(hi(1),hi(2),lo(3)))*dzinv 
  rhs(hi(1), hi(2), lo(3)) = RHS3(a(hi(1),hi(2),lo(3)), b, x, y, z, dxb, dyb, dzb)

  !xhi yhi zhi 
  z = zh - 0.5d0
  b = 1.d0/3.d0*(bx(hi(1),hi(2),hi(3)) + by(hi(1), hi(2), hi(3)) + bz(hi(1),hi(2),hi(3)))
  dxb = (bx(hi(1)+1,hi(2),hi(3)) - bx(hi(1),hi(2),hi(3)))*dxinv 
  dyb = (by(hi(1),hi(2)+1,hi(3)) - by(hi(1),hi(2),hi(3)))*dyinv
  dzb = (bz(hi(1),hi(2),hi(3)+1) - bz(hi(1),hi(2),hi(3)))*dzinv 
  rhs(hi(1), hi(2), hi(3)) = RHS3(a(hi(1),hi(2),hi(3)), b, x, y, z, dxb, dyb, dzb)
end subroutine build_rhs_3D

  function RHS3(a, b, x, y, z, dxb, dyb, dzb)
    real(amrex_real), intent(in) :: a, b, x, y, z, dxb, dyb, dzb 
    real(amrex_real)             :: RHS3
    real(amrex_real)             :: denom, term1, term2 
      denom = sqrt(x*x + y*y + z*z)
      term1 = a*(x+z)/denom 
      term2 = (x*x - x*z + y*y)*dzb + (y*y - x*z + z*z)*dxb - 2.d0*(x+z)*b -(x*y + y*z)*dyb 
      term2 = term2/(denom**3)
      RHS3  = term1 - term2 
  end function 

  function dbdx2(x,y)
  real(amrex_real), intent(in) :: x, y 
  real(amrex_real) :: dbdx2
    dbdx2 = cos(x)*cos(y) 
  end function dbdx2

  function dbdy2(x,y)
  real(amrex_real), intent(in) :: x, y
  real(amrex_real) :: dbdy2
    dbdy2 = -sin(x)*sin(y) 
  end function dbdy2

  function dbdx3(x,y,z)
  real(amrex_real), intent(in) :: x,y,z
  real(amrex_real) :: dbdx3
    dbdx3 = cos(x)*sin(y)*cos(z)
  end function dbdx3 

  function dbdy3(x,y,z)
  real(amrex_real), intent(in) :: x,y,z
  real(amrex_real) :: dbdy3
    dbdy3 = sin(x)*cos(y)*cos(z)
  end function dbdy3

  function dbdz3(x,y,z)
  real(amrex_real), intent(in) :: x,y,z
  real(amrex_real) :: dbdz3
    dbdz3 = -sin(x)*sin(y)*sin(z)
  end function dbdz3

end module rhs
