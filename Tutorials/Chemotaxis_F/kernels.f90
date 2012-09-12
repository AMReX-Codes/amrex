module chemotaxis_kernels
  !
  ! These kernels use point based reconstructions from polynomials and
  ! Gauss quadrature to perform flux and volume integrals.
  !
  implicit none
contains


  !
  ! Cell motility
  !
  subroutine cell_motility(f, u, lo, hi, ng, dx, diff) &
       bind(c, name='cell_motility')
    use iso_c_binding
    integer(c_int), intent(in)        :: lo(2), hi(2)
    integer(c_int), intent(in), value :: ng
    real(c_double), intent(in), value :: dx, diff
    real(c_double), intent(in)        :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(inout)     :: f(lo(1):hi(1),lo(2):hi(2))

    integer :: i, j

    real(c_double) :: u_edge, u_x_edge, u_y_edge

    ! this kernel was generated

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
u_x_edge = (-u(i+0,j+0) + u(i+1,j+0))/dx
u_edge = u(i+0,j+0)/2 + u(i+1,j+0)/2
f(i, j) = f(i, j) + (diff*u_edge*u_x_edge)/dx
u_x_edge = (u(i+0,j+0) - u(i-1,j+0))/dx
u_edge = u(i+0,j+0)/2 + u(i-1,j+0)/2
f(i, j) = f(i, j) - (diff*u_edge*u_x_edge)/dx
u_y_edge = (-u(i+0,j+0) + u(i+0,j+1))/dx
u_edge = u(i+0,j+0)/2 + u(i+0,j+1)/2
f(i, j) = f(i, j) + (diff*u_edge*u_y_edge)/dx
u_y_edge = (u(i+0,j+0) - u(i+0,j-1))/dx
u_edge = u(i+0,j+0)/2 + u(i+0,j-1)/2
f(i, j) = f(i, j) - (diff*u_edge*u_y_edge)/dx      
       end do
    end do

  end subroutine cell_motility


  !
  ! Chemotactic sensitivity
  !
  subroutine chemotactic_sensitivity(f, u, v, lo, hi, ng, dx, chi, alpha, gamma) &
       bind(c, name='chemotactic_sensitivity')
    use iso_c_binding
    integer(c_int), intent(in)        :: lo(2), hi(2)
    integer(c_int), intent(in), value :: ng
    real(c_double), intent(in), value :: dx, chi, alpha, gamma
    real(c_double), intent(in)        :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(in)        :: v(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(inout)     :: f(lo(1):hi(1),lo(2):hi(2))

    integer :: i, j

    real(c_double) :: u_edge, v_edge, v_x_edge, v_y_edge

    ! this kernel was generated

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
v_edge = v(i+0,j+0)/2 + v(i+1,j+0)/2
u_edge = u(i+0,j+0)/2 + u(i+1,j+0)/2
v_x_edge = (-v(i+0,j+0) + v(i+1,j+0))/dx
f(i, j) = f(i, j) + (-chi*u_edge*v_x_edge*(1 - u_edge/gamma)/(alpha*v_edge + 1)**2)/dx
v_edge = v(i+0,j+0)/2 + v(i-1,j+0)/2
u_edge = u(i+0,j+0)/2 + u(i-1,j+0)/2
v_x_edge = (v(i+0,j+0) - v(i-1,j+0))/dx
f(i, j) = f(i, j) - (-chi*u_edge*v_x_edge*(1 - u_edge/gamma)/(alpha*v_edge + 1)**2)/dx
v_edge = v(i+0,j+0)/2 + v(i+0,j+1)/2
u_edge = u(i+0,j+0)/2 + u(i+0,j+1)/2
v_y_edge = (-v(i+0,j+0) + v(i+0,j+1))/dx
f(i, j) = f(i, j) + (-chi*u_edge*v_y_edge*(1 - u_edge/gamma)/(alpha*v_edge + 1)**2)/dx
v_edge = v(i+0,j+0)/2 + v(i+0,j-1)/2
u_edge = u(i+0,j+0)/2 + u(i+0,j-1)/2
v_y_edge = (v(i+0,j+0) - v(i+0,j-1))/dx
f(i, j) = f(i, j) - (-chi*u_edge*v_y_edge*(1 - u_edge/gamma)/(alpha*v_edge + 1)**2)/dx
       end do
    end do

  end subroutine chemotactic_sensitivity


  !
  ! Signal diffusion
  !
  subroutine signal_diffusion(f, v, lo, hi, ng, dx) &
       bind(c, name='signal_diffusion')
    use iso_c_binding
    integer(c_int), intent(in)        :: lo(2), hi(2)
    integer(c_int), intent(in), value :: ng
    real(c_double), intent(in), value :: dx
    real(c_double), intent(in)        :: v(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(inout)     :: f(lo(1):hi(1),lo(2):hi(2))

    integer :: i, j

    real(c_double) :: v_x_edge, v_y_edge

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
v_x_edge = (-v(i+0,j+0) + v(i+1,j+0))/dx
f(i, j) = f(i, j) + (v_x_edge)/dx
v_x_edge = (v(i+0,j+0) - v(i-1,j+0))/dx
f(i, j) = f(i, j) - (v_x_edge)/dx
v_y_edge = (-v(i+0,j+0) + v(i+0,j+1))/dx
f(i, j) = f(i, j) + (v_y_edge)/dx
v_y_edge = (v(i+0,j+0) - v(i+0,j-1))/dx
f(i, j) = f(i, j) - (v_y_edge)/dx
       end do
    end do

  end subroutine signal_diffusion


  !
  ! Signal production
  !
  subroutine signal_production(f, u, v, lo, hi, ng, dx, phi) &
       bind(c, name='signal_production')
    use iso_c_binding
    integer(c_int), intent(in)        :: lo(2), hi(2)
    integer(c_int), intent(in), value :: ng
    real(c_double), intent(in), value :: dx, phi
    real(c_double), intent(in)        :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(in)        :: v(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(inout)     :: f(lo(1):hi(1),lo(2):hi(2))

    integer :: i, j

    ! this kernel is trivial!

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          ! f(i, j) = f(i, j) + u(i, j) / (1.0d0 + phi * u(i, j))
          f(i, j) = f(i, j) + u(i, j)
       end do
    end do

  end subroutine signal_production


  !
  ! Signal degredation
  !
  subroutine signal_degradation(f, u, v, lo, hi, ng, dx) &
       bind(c, name='signal_degradation')
    use iso_c_binding
    integer(c_int), intent(in)        :: lo(2), hi(2)
    integer(c_int), intent(in), value :: ng
    real(c_double), intent(in), value :: dx
    real(c_double), intent(in)        :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(in)        :: v(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(inout)     :: f(lo(1):hi(1),lo(2):hi(2))

    integer :: i, j

    ! this kernel is trivial!

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          f(i, j) = f(i, j) - v(i, j)
       end do
    end do

  end subroutine signal_degradation

end module chemotaxis_kernels
