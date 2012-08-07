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

    real(c_double), dimension(-2:2) :: cof, vec, tmp
    real(c_double), dimension(0:2)  :: qwts, qvls, u_edge, u_x_edge, u_y_edge

    ! this kernel was generated

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
{motility_kernel}      
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

    real(c_double), dimension(-2:2) :: cof, vec, tmp
    real(c_double), dimension(0:2)  :: qwts, qvls, u_edge, v_edge, v_x_edge, v_y_edge

    ! this kernel was generated

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
{sensitivity_kernel}
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

    real(c_double), dimension(-2:2) :: cof, vec, tmp
    real(c_double), dimension(0:2)  :: qwts, qvls, v_x_edge, v_y_edge

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
{diffusion_kernel}
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

    real(c_double), dimension(-2:2)     :: cof, vec
    real(c_double), dimension(-2:2,0:2) :: tmp
    real(c_double), dimension(0:2)      :: qwts, qvls1
    real(c_double), dimension(0:2,0:2)  :: qvls2, u_interior

    ! this kernel was generated

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
{production_kernel}
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
