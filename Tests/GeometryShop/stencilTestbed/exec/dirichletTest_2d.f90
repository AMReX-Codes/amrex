
module stencil_test_module

  use bl_fort_module, only : dim=>bl_spacedim
  implicit none

  private

  public :: get_bndry_grad_stencil, apply_bndry_grad_stencil, init_phi

contains

  subroutine get_bndry_grad_stencil(stencil,stencilB,baseiv,normals,centroids,ivs,num) bind(C,name="get_bndry_grad_stencil")
    integer, intent(in) :: num
    double precision, intent(inout) :: stencil(0:2,0:2,num)
    double precision, intent(inout) ::      stencilB(num)
    integer,          intent(inout) ::    baseiv(dim,num)
    double precision, intent(in   ) ::   normals(dim,num)
    double precision, intent(in   ) :: centroids(dim,num)
    integer,          intent(in   ) ::       ivs(dim,num)
    integer :: i, c(dim), s(dim), iv(dim)
    double precision :: n(dim), b(dim), x(2), y(2), d(2)

    ! Rotate work to 1st quadrant, where 1-component is "indepenent".  Find 2 quadratically 
    ! interpolated values at x1=1 and x1=2 along ray normal to boundary, and intersecting
    ! eb at centroid.  These two values, plus the Dirichlet value at the eb are then used
    ! to compute a normal gradient at the eb.  The resulting stencil contains 6+1 points, and 
    ! is stored in a 3x3 matrix plus an aux value.  The normal gradient will then be evaluated:
    !
    !      dudn(i) = u(baseiv(i):baseiv(i)+(2,2)) . stencil(0:2,0:2,i)
    !                    + bcval(ivs(i)) . stencilB(i)
    !
    stencil(:,:,1:num) = 0.d0
    baseiv(:,1:num) = 0    

    do i=1,num
       
       n(:) = normals(:,i)

       c(1) = 1 + INT( MIN(1.d0, ABS(n(2)/n(1)) ) ) ! Either 1 or 2 only
       c(2) = 3 - c(1) ! Either 1 or 2 only

       s(:) = SIGN(1.d0,n(c(:)))
       b(:) = centroids(c(:),i) * s(:)

       baseiv(:,i) = ivs(:,i) - 1 + s(c(:))
             
       x(1) = 1
       x(2) = 2
       y(:) = b(2) + (x(:) - b(1))*ABS(n(c(2))/n(c(1)))
       d(:) = SQRT( (x(:)-b(1))**2 + (y(:)-b(2))**2)

       iv(c(1)) = 1 * s(1) + ivs(c(1),i) - baseiv(c(1),i)
       iv(c(2)) = (1 - 1)*s(2) + ivs(c(2),i) - baseiv(c(2),i)
       stencil(iv(1),iv(2),i) = (0.5d0 * (y(1)-1.d0) * (y(1)-2.d0)) * d(2)/(d(1)*(d(2)-d(1)))

       iv(c(2)) = (2 - 1)*s(2) + ivs(c(2),i) - baseiv(c(2),i)
       stencil(iv(1),iv(2),i) = (- y(1) * (y(1)-2.d0))              * d(2)/(d(1)*(d(2)-d(1)))

       iv(c(2)) = (3 - 1)*s(2) + ivs(c(2),i) - baseiv(c(2),i)
       stencil(iv(1),iv(2),i) = (0.5d0 * y(1) * (y(1)-1.d0))        * d(2)/(d(1)*(d(2)-d(1)))


       iv(c(1)) = 2 * s(1) + ivs(c(1),i) - baseiv(c(1),i)  
       iv(c(2)) = (1 - 1)*s(2) + ivs(c(2),i) - baseiv(c(2),i)
       stencil(iv(1),iv(2),i) = (0.5d0 * (y(2)-1.d0) * (y(2)-2.d0)) * (-d(1)/(d(2)*(d(2)-d(1))))

       iv(c(2)) = (2 - 1)*s(2) + ivs(c(2),i) - baseiv(c(2),i)
       stencil(iv(1),iv(2),i) = (- y(2) * (y(2)-2.d0))              * (-d(1)/(d(2)*(d(2)-d(1))))

       iv(c(2)) = (3 - 1)*s(2) + ivs(c(2),i) - baseiv(c(2),i)
       stencil(iv(1),iv(2),i) = (0.5d0 * y(2) * (y(2)-1.d0))        * (-d(1)/(d(2)*(d(2)-d(1))))

       stencilB(i) = - (d(2)+d(1))/(d(2)*d(1))
       
    enddo

  end subroutine get_bndry_grad_stencil


  subroutine apply_bndry_grad_stencil(dudn, b,&
       u,u_l1,u_l2,u_h1,u_h2, &
       stencil,stencilB,baseiv,num,dx) bind(C,name="apply_bndry_grad_stencil")
    integer, intent(in) :: num
    double precision, intent(inout) :: dudn(num)
    double precision, intent(in   ) :: b(num)
    integer,          intent(in   ) :: u_l1,u_l2,u_h1,u_h2
    double precision, intent(in   ) :: u(u_l1:u_h1,u_l2:u_h2)
    double precision, intent(in   ) :: stencil(0:2,0:2,num)
    double precision, intent(in   ) :: stencilB(num)
    integer,          intent(in   ) :: baseiv(dim,num)
    double precision, intent(in   ) :: dx

    double precision :: fac
    integer :: i, j, L

    fac = 1.d0 / dx

    do L=1,num
       
       i = baseiv(1,L)
       j = baseiv(2,L)

       dudn(L) = fac * (b(L) * stencilB(L) &
            + u(i  ,j  ) * stencil(0,0,L) &
            + u(i+1,j  ) * stencil(1,0,L) &
            + u(i+2,j  ) * stencil(2,0,L) &
            + u(i  ,j+1) * stencil(0,1,L) & 
            + u(i+1,j+1) * stencil(1,1,L) &
            + u(i+2,j+1) * stencil(2,1,L) &
            + u(i  ,j+2) * stencil(0,2,L) &
            + u(i+1,j+2) * stencil(1,2,L) &
            + u(i+2,j+2) * stencil(2,2,L) )
    enddo

  end subroutine apply_bndry_grad_stencil


  subroutine init_phi(&
       phi,phi_l1,phi_l2,phi_h1,phi_h2, &
       lo, hi, plo, dx) bind(C,name="init_phi")

    integer, intent(in) :: phi_l1,phi_l2,phi_h1,phi_h2

    integer, intent(in) :: lo(2),hi(2)

    double precision, intent(inout) :: phi(phi_l1:phi_h1,phi_l2:phi_h2)
    double precision, intent(in   ) :: plo(2), dx(2)

    integer :: i,j
    double precision :: x,y

    do j=lo(2),hi(2)
       y = plo(2) + (j+0.5d0)*dx(2)
       do i=lo(1),hi(1)
          x = plo(1) + (i+0.5d0)*dx(1)
          phi(i,j) = 0.5d0 * (x-0.5d0)**2
       enddo
    enddo

  end subroutine init_phi

end module stencil_test_module
