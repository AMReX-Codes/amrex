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
qwts = [ 0.555555555555556d0, 0.888888888888889d0, 0.555555555555556d0 ]
cof = [ 0.0d0, 0.0833333333333333d0/dx, -1.25d0/dx, 1.25d0/dx, -0.0833333333333333d0/dx ]
vec = [ u(i-2,j-2), u(i-1,j-2), u(i+0,j-2), u(i+1,j-2), u(i+2,j-2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ u(i-2,j-1), u(i-1,j-1), u(i+0,j-1), u(i+1,j-1), u(i+2,j-1) ]; tmp(-1) = dot_product(cof, vec)
vec = [ u(i-2,j+0), u(i-1,j+0), u(i+0,j+0), u(i+1,j+0), u(i+2,j+0) ]; tmp(+0) = dot_product(cof, vec)
vec = [ u(i-2,j+1), u(i-1,j+1), u(i+0,j+1), u(i+1,j+1), u(i+2,j+1) ]; tmp(+1) = dot_product(cof, vec)
vec = [ u(i-2,j+2), u(i-1,j+2), u(i+0,j+2), u(i+1,j+2), u(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_x_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_x_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_x_edge(2) = dot_product(cof, tmp)
cof = [ 0.0333333333333333d0, -0.216666666666667d0, 0.783333333333333d0, 0.450000000000000d0, -0.0500000000000000d0 ]
vec = [ u(i-2,j-2), u(i-1,j-2), u(i+0,j-2), u(i+1,j-2), u(i+2,j-2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ u(i-2,j-1), u(i-1,j-1), u(i+0,j-1), u(i+1,j-1), u(i+2,j-1) ]; tmp(-1) = dot_product(cof, vec)
vec = [ u(i-2,j+0), u(i-1,j+0), u(i+0,j+0), u(i+1,j+0), u(i+2,j+0) ]; tmp(+0) = dot_product(cof, vec)
vec = [ u(i-2,j+1), u(i-1,j+1), u(i+0,j+1), u(i+1,j+1), u(i+2,j+1) ]; tmp(+1) = dot_product(cof, vec)
vec = [ u(i-2,j+2), u(i-1,j+2), u(i+0,j+2), u(i+1,j+2), u(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_edge(2) = dot_product(cof, tmp)
qvls = diff*u_edge*u_x_edge
f(i, j) = f(i, j) + 0.5d0*dx*dot_product(qwts, qvls)
cof = [ 0.0833333333333333d0/dx, -1.25d0/dx, 1.25d0/dx, -0.0833333333333333d0/dx, 0.0d0 ]
vec = [ u(i-2,j-2), u(i-1,j-2), u(i+0,j-2), u(i+1,j-2), u(i+2,j-2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ u(i-2,j-1), u(i-1,j-1), u(i+0,j-1), u(i+1,j-1), u(i+2,j-1) ]; tmp(-1) = dot_product(cof, vec)
vec = [ u(i-2,j+0), u(i-1,j+0), u(i+0,j+0), u(i+1,j+0), u(i+2,j+0) ]; tmp(+0) = dot_product(cof, vec)
vec = [ u(i-2,j+1), u(i-1,j+1), u(i+0,j+1), u(i+1,j+1), u(i+2,j+1) ]; tmp(+1) = dot_product(cof, vec)
vec = [ u(i-2,j+2), u(i-1,j+2), u(i+0,j+2), u(i+1,j+2), u(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_x_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_x_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_x_edge(2) = dot_product(cof, tmp)
cof = [ -0.0500000000000000d0, 0.450000000000000d0, 0.783333333333333d0, -0.216666666666667d0, 0.0333333333333333d0 ]
vec = [ u(i-2,j-2), u(i-1,j-2), u(i+0,j-2), u(i+1,j-2), u(i+2,j-2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ u(i-2,j-1), u(i-1,j-1), u(i+0,j-1), u(i+1,j-1), u(i+2,j-1) ]; tmp(-1) = dot_product(cof, vec)
vec = [ u(i-2,j+0), u(i-1,j+0), u(i+0,j+0), u(i+1,j+0), u(i+2,j+0) ]; tmp(+0) = dot_product(cof, vec)
vec = [ u(i-2,j+1), u(i-1,j+1), u(i+0,j+1), u(i+1,j+1), u(i+2,j+1) ]; tmp(+1) = dot_product(cof, vec)
vec = [ u(i-2,j+2), u(i-1,j+2), u(i+0,j+2), u(i+1,j+2), u(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_edge(2) = dot_product(cof, tmp)
qvls = diff*u_edge*u_x_edge
f(i, j) = f(i, j) - 0.5d0*dx*dot_product(qwts, qvls)
cof = [ 0.0d0, 0.0833333333333333d0/dx, -1.25d0/dx, 1.25d0/dx, -0.0833333333333333d0/dx ]
vec = [ u(i-2,j-2), u(i-2,j-1), u(i-2,j+0), u(i-2,j+1), u(i-2,j+2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ u(i-1,j-2), u(i-1,j-1), u(i-1,j+0), u(i-1,j+1), u(i-1,j+2) ]; tmp(-1) = dot_product(cof, vec)
vec = [ u(i+0,j-2), u(i+0,j-1), u(i+0,j+0), u(i+0,j+1), u(i+0,j+2) ]; tmp(+0) = dot_product(cof, vec)
vec = [ u(i+1,j-2), u(i+1,j-1), u(i+1,j+0), u(i+1,j+1), u(i+1,j+2) ]; tmp(+1) = dot_product(cof, vec)
vec = [ u(i+2,j-2), u(i+2,j-1), u(i+2,j+0), u(i+2,j+1), u(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_y_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_y_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_y_edge(2) = dot_product(cof, tmp)
cof = [ 0.0333333333333333d0, -0.216666666666667d0, 0.783333333333333d0, 0.450000000000000d0, -0.0500000000000000d0 ]
vec = [ u(i-2,j-2), u(i-2,j-1), u(i-2,j+0), u(i-2,j+1), u(i-2,j+2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ u(i-1,j-2), u(i-1,j-1), u(i-1,j+0), u(i-1,j+1), u(i-1,j+2) ]; tmp(-1) = dot_product(cof, vec)
vec = [ u(i+0,j-2), u(i+0,j-1), u(i+0,j+0), u(i+0,j+1), u(i+0,j+2) ]; tmp(+0) = dot_product(cof, vec)
vec = [ u(i+1,j-2), u(i+1,j-1), u(i+1,j+0), u(i+1,j+1), u(i+1,j+2) ]; tmp(+1) = dot_product(cof, vec)
vec = [ u(i+2,j-2), u(i+2,j-1), u(i+2,j+0), u(i+2,j+1), u(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_edge(2) = dot_product(cof, tmp)
qvls = diff*u_edge*u_y_edge
f(i, j) = f(i, j) + 0.5d0*dx*dot_product(qwts, qvls)
cof = [ 0.0833333333333333d0/dx, -1.25d0/dx, 1.25d0/dx, -0.0833333333333333d0/dx, 0.0d0 ]
vec = [ u(i-2,j-2), u(i-2,j-1), u(i-2,j+0), u(i-2,j+1), u(i-2,j+2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ u(i-1,j-2), u(i-1,j-1), u(i-1,j+0), u(i-1,j+1), u(i-1,j+2) ]; tmp(-1) = dot_product(cof, vec)
vec = [ u(i+0,j-2), u(i+0,j-1), u(i+0,j+0), u(i+0,j+1), u(i+0,j+2) ]; tmp(+0) = dot_product(cof, vec)
vec = [ u(i+1,j-2), u(i+1,j-1), u(i+1,j+0), u(i+1,j+1), u(i+1,j+2) ]; tmp(+1) = dot_product(cof, vec)
vec = [ u(i+2,j-2), u(i+2,j-1), u(i+2,j+0), u(i+2,j+1), u(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_y_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_y_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_y_edge(2) = dot_product(cof, tmp)
cof = [ -0.0500000000000000d0, 0.450000000000000d0, 0.783333333333333d0, -0.216666666666667d0, 0.0333333333333333d0 ]
vec = [ u(i-2,j-2), u(i-2,j-1), u(i-2,j+0), u(i-2,j+1), u(i-2,j+2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ u(i-1,j-2), u(i-1,j-1), u(i-1,j+0), u(i-1,j+1), u(i-1,j+2) ]; tmp(-1) = dot_product(cof, vec)
vec = [ u(i+0,j-2), u(i+0,j-1), u(i+0,j+0), u(i+0,j+1), u(i+0,j+2) ]; tmp(+0) = dot_product(cof, vec)
vec = [ u(i+1,j-2), u(i+1,j-1), u(i+1,j+0), u(i+1,j+1), u(i+1,j+2) ]; tmp(+1) = dot_product(cof, vec)
vec = [ u(i+2,j-2), u(i+2,j-1), u(i+2,j+0), u(i+2,j+1), u(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_edge(2) = dot_product(cof, tmp)
qvls = diff*u_edge*u_y_edge
f(i, j) = f(i, j) - 0.5d0*dx*dot_product(qwts, qvls)      
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
qwts = [ 0.555555555555556d0, 0.888888888888889d0, 0.555555555555556d0 ]
cof = [ 0.0333333333333333d0, -0.216666666666667d0, 0.783333333333333d0, 0.450000000000000d0, -0.0500000000000000d0 ]
vec = [ v(i-2,j-2), v(i-1,j-2), v(i+0,j-2), v(i+1,j-2), v(i+2,j-2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ v(i-2,j-1), v(i-1,j-1), v(i+0,j-1), v(i+1,j-1), v(i+2,j-1) ]; tmp(-1) = dot_product(cof, vec)
vec = [ v(i-2,j+0), v(i-1,j+0), v(i+0,j+0), v(i+1,j+0), v(i+2,j+0) ]; tmp(+0) = dot_product(cof, vec)
vec = [ v(i-2,j+1), v(i-1,j+1), v(i+0,j+1), v(i+1,j+1), v(i+2,j+1) ]; tmp(+1) = dot_product(cof, vec)
vec = [ v(i-2,j+2), v(i-1,j+2), v(i+0,j+2), v(i+1,j+2), v(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
v_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
v_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
v_edge(2) = dot_product(cof, tmp)
cof = [ 0.0333333333333333d0, -0.216666666666667d0, 0.783333333333333d0, 0.450000000000000d0, -0.0500000000000000d0 ]
vec = [ u(i-2,j-2), u(i-1,j-2), u(i+0,j-2), u(i+1,j-2), u(i+2,j-2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ u(i-2,j-1), u(i-1,j-1), u(i+0,j-1), u(i+1,j-1), u(i+2,j-1) ]; tmp(-1) = dot_product(cof, vec)
vec = [ u(i-2,j+0), u(i-1,j+0), u(i+0,j+0), u(i+1,j+0), u(i+2,j+0) ]; tmp(+0) = dot_product(cof, vec)
vec = [ u(i-2,j+1), u(i-1,j+1), u(i+0,j+1), u(i+1,j+1), u(i+2,j+1) ]; tmp(+1) = dot_product(cof, vec)
vec = [ u(i-2,j+2), u(i-1,j+2), u(i+0,j+2), u(i+1,j+2), u(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_edge(2) = dot_product(cof, tmp)
cof = [ 0.0d0, 0.0833333333333333d0/dx, -1.25d0/dx, 1.25d0/dx, -0.0833333333333333d0/dx ]
vec = [ v(i-2,j-2), v(i-1,j-2), v(i+0,j-2), v(i+1,j-2), v(i+2,j-2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ v(i-2,j-1), v(i-1,j-1), v(i+0,j-1), v(i+1,j-1), v(i+2,j-1) ]; tmp(-1) = dot_product(cof, vec)
vec = [ v(i-2,j+0), v(i-1,j+0), v(i+0,j+0), v(i+1,j+0), v(i+2,j+0) ]; tmp(+0) = dot_product(cof, vec)
vec = [ v(i-2,j+1), v(i-1,j+1), v(i+0,j+1), v(i+1,j+1), v(i+2,j+1) ]; tmp(+1) = dot_product(cof, vec)
vec = [ v(i-2,j+2), v(i-1,j+2), v(i+0,j+2), v(i+1,j+2), v(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
v_x_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
v_x_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
v_x_edge(2) = dot_product(cof, tmp)
qvls = chi*u_edge*v_x_edge*(1 - u_edge/gamma)/(alpha*v_edge + 1)**2
f(i, j) = f(i, j) + 0.5d0*dx*dot_product(qwts, qvls)
cof = [ -0.0500000000000000d0, 0.450000000000000d0, 0.783333333333333d0, -0.216666666666667d0, 0.0333333333333333d0 ]
vec = [ v(i-2,j-2), v(i-1,j-2), v(i+0,j-2), v(i+1,j-2), v(i+2,j-2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ v(i-2,j-1), v(i-1,j-1), v(i+0,j-1), v(i+1,j-1), v(i+2,j-1) ]; tmp(-1) = dot_product(cof, vec)
vec = [ v(i-2,j+0), v(i-1,j+0), v(i+0,j+0), v(i+1,j+0), v(i+2,j+0) ]; tmp(+0) = dot_product(cof, vec)
vec = [ v(i-2,j+1), v(i-1,j+1), v(i+0,j+1), v(i+1,j+1), v(i+2,j+1) ]; tmp(+1) = dot_product(cof, vec)
vec = [ v(i-2,j+2), v(i-1,j+2), v(i+0,j+2), v(i+1,j+2), v(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
v_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
v_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
v_edge(2) = dot_product(cof, tmp)
cof = [ -0.0500000000000000d0, 0.450000000000000d0, 0.783333333333333d0, -0.216666666666667d0, 0.0333333333333333d0 ]
vec = [ u(i-2,j-2), u(i-1,j-2), u(i+0,j-2), u(i+1,j-2), u(i+2,j-2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ u(i-2,j-1), u(i-1,j-1), u(i+0,j-1), u(i+1,j-1), u(i+2,j-1) ]; tmp(-1) = dot_product(cof, vec)
vec = [ u(i-2,j+0), u(i-1,j+0), u(i+0,j+0), u(i+1,j+0), u(i+2,j+0) ]; tmp(+0) = dot_product(cof, vec)
vec = [ u(i-2,j+1), u(i-1,j+1), u(i+0,j+1), u(i+1,j+1), u(i+2,j+1) ]; tmp(+1) = dot_product(cof, vec)
vec = [ u(i-2,j+2), u(i-1,j+2), u(i+0,j+2), u(i+1,j+2), u(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_edge(2) = dot_product(cof, tmp)
cof = [ 0.0833333333333333d0/dx, -1.25d0/dx, 1.25d0/dx, -0.0833333333333333d0/dx, 0.0d0 ]
vec = [ v(i-2,j-2), v(i-1,j-2), v(i+0,j-2), v(i+1,j-2), v(i+2,j-2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ v(i-2,j-1), v(i-1,j-1), v(i+0,j-1), v(i+1,j-1), v(i+2,j-1) ]; tmp(-1) = dot_product(cof, vec)
vec = [ v(i-2,j+0), v(i-1,j+0), v(i+0,j+0), v(i+1,j+0), v(i+2,j+0) ]; tmp(+0) = dot_product(cof, vec)
vec = [ v(i-2,j+1), v(i-1,j+1), v(i+0,j+1), v(i+1,j+1), v(i+2,j+1) ]; tmp(+1) = dot_product(cof, vec)
vec = [ v(i-2,j+2), v(i-1,j+2), v(i+0,j+2), v(i+1,j+2), v(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
v_x_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
v_x_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
v_x_edge(2) = dot_product(cof, tmp)
qvls = chi*u_edge*v_x_edge*(1 - u_edge/gamma)/(alpha*v_edge + 1)**2
f(i, j) = f(i, j) - 0.5d0*dx*dot_product(qwts, qvls)
cof = [ 0.0333333333333333d0, -0.216666666666667d0, 0.783333333333333d0, 0.450000000000000d0, -0.0500000000000000d0 ]
vec = [ v(i-2,j-2), v(i-2,j-1), v(i-2,j+0), v(i-2,j+1), v(i-2,j+2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ v(i-1,j-2), v(i-1,j-1), v(i-1,j+0), v(i-1,j+1), v(i-1,j+2) ]; tmp(-1) = dot_product(cof, vec)
vec = [ v(i+0,j-2), v(i+0,j-1), v(i+0,j+0), v(i+0,j+1), v(i+0,j+2) ]; tmp(+0) = dot_product(cof, vec)
vec = [ v(i+1,j-2), v(i+1,j-1), v(i+1,j+0), v(i+1,j+1), v(i+1,j+2) ]; tmp(+1) = dot_product(cof, vec)
vec = [ v(i+2,j-2), v(i+2,j-1), v(i+2,j+0), v(i+2,j+1), v(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
v_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
v_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
v_edge(2) = dot_product(cof, tmp)
cof = [ 0.0333333333333333d0, -0.216666666666667d0, 0.783333333333333d0, 0.450000000000000d0, -0.0500000000000000d0 ]
vec = [ u(i-2,j-2), u(i-2,j-1), u(i-2,j+0), u(i-2,j+1), u(i-2,j+2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ u(i-1,j-2), u(i-1,j-1), u(i-1,j+0), u(i-1,j+1), u(i-1,j+2) ]; tmp(-1) = dot_product(cof, vec)
vec = [ u(i+0,j-2), u(i+0,j-1), u(i+0,j+0), u(i+0,j+1), u(i+0,j+2) ]; tmp(+0) = dot_product(cof, vec)
vec = [ u(i+1,j-2), u(i+1,j-1), u(i+1,j+0), u(i+1,j+1), u(i+1,j+2) ]; tmp(+1) = dot_product(cof, vec)
vec = [ u(i+2,j-2), u(i+2,j-1), u(i+2,j+0), u(i+2,j+1), u(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_edge(2) = dot_product(cof, tmp)
cof = [ 0.0d0, 0.0833333333333333d0/dx, -1.25d0/dx, 1.25d0/dx, -0.0833333333333333d0/dx ]
vec = [ v(i-2,j-2), v(i-2,j-1), v(i-2,j+0), v(i-2,j+1), v(i-2,j+2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ v(i-1,j-2), v(i-1,j-1), v(i-1,j+0), v(i-1,j+1), v(i-1,j+2) ]; tmp(-1) = dot_product(cof, vec)
vec = [ v(i+0,j-2), v(i+0,j-1), v(i+0,j+0), v(i+0,j+1), v(i+0,j+2) ]; tmp(+0) = dot_product(cof, vec)
vec = [ v(i+1,j-2), v(i+1,j-1), v(i+1,j+0), v(i+1,j+1), v(i+1,j+2) ]; tmp(+1) = dot_product(cof, vec)
vec = [ v(i+2,j-2), v(i+2,j-1), v(i+2,j+0), v(i+2,j+1), v(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
v_y_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
v_y_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
v_y_edge(2) = dot_product(cof, tmp)
qvls = chi*u_edge*v_y_edge*(1 - u_edge/gamma)/(alpha*v_edge + 1)**2
f(i, j) = f(i, j) + 0.5d0*dx*dot_product(qwts, qvls)
cof = [ -0.0500000000000000d0, 0.450000000000000d0, 0.783333333333333d0, -0.216666666666667d0, 0.0333333333333333d0 ]
vec = [ v(i-2,j-2), v(i-2,j-1), v(i-2,j+0), v(i-2,j+1), v(i-2,j+2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ v(i-1,j-2), v(i-1,j-1), v(i-1,j+0), v(i-1,j+1), v(i-1,j+2) ]; tmp(-1) = dot_product(cof, vec)
vec = [ v(i+0,j-2), v(i+0,j-1), v(i+0,j+0), v(i+0,j+1), v(i+0,j+2) ]; tmp(+0) = dot_product(cof, vec)
vec = [ v(i+1,j-2), v(i+1,j-1), v(i+1,j+0), v(i+1,j+1), v(i+1,j+2) ]; tmp(+1) = dot_product(cof, vec)
vec = [ v(i+2,j-2), v(i+2,j-1), v(i+2,j+0), v(i+2,j+1), v(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
v_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
v_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
v_edge(2) = dot_product(cof, tmp)
cof = [ -0.0500000000000000d0, 0.450000000000000d0, 0.783333333333333d0, -0.216666666666667d0, 0.0333333333333333d0 ]
vec = [ u(i-2,j-2), u(i-2,j-1), u(i-2,j+0), u(i-2,j+1), u(i-2,j+2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ u(i-1,j-2), u(i-1,j-1), u(i-1,j+0), u(i-1,j+1), u(i-1,j+2) ]; tmp(-1) = dot_product(cof, vec)
vec = [ u(i+0,j-2), u(i+0,j-1), u(i+0,j+0), u(i+0,j+1), u(i+0,j+2) ]; tmp(+0) = dot_product(cof, vec)
vec = [ u(i+1,j-2), u(i+1,j-1), u(i+1,j+0), u(i+1,j+1), u(i+1,j+2) ]; tmp(+1) = dot_product(cof, vec)
vec = [ u(i+2,j-2), u(i+2,j-1), u(i+2,j+0), u(i+2,j+1), u(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_edge(2) = dot_product(cof, tmp)
cof = [ 0.0833333333333333d0/dx, -1.25d0/dx, 1.25d0/dx, -0.0833333333333333d0/dx, 0.0d0 ]
vec = [ v(i-2,j-2), v(i-2,j-1), v(i-2,j+0), v(i-2,j+1), v(i-2,j+2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ v(i-1,j-2), v(i-1,j-1), v(i-1,j+0), v(i-1,j+1), v(i-1,j+2) ]; tmp(-1) = dot_product(cof, vec)
vec = [ v(i+0,j-2), v(i+0,j-1), v(i+0,j+0), v(i+0,j+1), v(i+0,j+2) ]; tmp(+0) = dot_product(cof, vec)
vec = [ v(i+1,j-2), v(i+1,j-1), v(i+1,j+0), v(i+1,j+1), v(i+1,j+2) ]; tmp(+1) = dot_product(cof, vec)
vec = [ v(i+2,j-2), v(i+2,j-1), v(i+2,j+0), v(i+2,j+1), v(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
v_y_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
v_y_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
v_y_edge(2) = dot_product(cof, tmp)
qvls = chi*u_edge*v_y_edge*(1 - u_edge/gamma)/(alpha*v_edge + 1)**2
f(i, j) = f(i, j) - 0.5d0*dx*dot_product(qwts, qvls)
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
qwts = [ 0.555555555555556d0, 0.888888888888889d0, 0.555555555555556d0 ]
cof = [ 0.0d0, 0.0833333333333333d0/dx, -1.25d0/dx, 1.25d0/dx, -0.0833333333333333d0/dx ]
vec = [ v(i-2,j-2), v(i-1,j-2), v(i+0,j-2), v(i+1,j-2), v(i+2,j-2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ v(i-2,j-1), v(i-1,j-1), v(i+0,j-1), v(i+1,j-1), v(i+2,j-1) ]; tmp(-1) = dot_product(cof, vec)
vec = [ v(i-2,j+0), v(i-1,j+0), v(i+0,j+0), v(i+1,j+0), v(i+2,j+0) ]; tmp(+0) = dot_product(cof, vec)
vec = [ v(i-2,j+1), v(i-1,j+1), v(i+0,j+1), v(i+1,j+1), v(i+2,j+1) ]; tmp(+1) = dot_product(cof, vec)
vec = [ v(i-2,j+2), v(i-1,j+2), v(i+0,j+2), v(i+1,j+2), v(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
v_x_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
v_x_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
v_x_edge(2) = dot_product(cof, tmp)
qvls = v_x_edge
f(i, j) = f(i, j) + 0.5d0*dx*dot_product(qwts, qvls)
cof = [ 0.0833333333333333d0/dx, -1.25d0/dx, 1.25d0/dx, -0.0833333333333333d0/dx, 0.0d0 ]
vec = [ v(i-2,j-2), v(i-1,j-2), v(i+0,j-2), v(i+1,j-2), v(i+2,j-2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ v(i-2,j-1), v(i-1,j-1), v(i+0,j-1), v(i+1,j-1), v(i+2,j-1) ]; tmp(-1) = dot_product(cof, vec)
vec = [ v(i-2,j+0), v(i-1,j+0), v(i+0,j+0), v(i+1,j+0), v(i+2,j+0) ]; tmp(+0) = dot_product(cof, vec)
vec = [ v(i-2,j+1), v(i-1,j+1), v(i+0,j+1), v(i+1,j+1), v(i+2,j+1) ]; tmp(+1) = dot_product(cof, vec)
vec = [ v(i-2,j+2), v(i-1,j+2), v(i+0,j+2), v(i+1,j+2), v(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
v_x_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
v_x_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
v_x_edge(2) = dot_product(cof, tmp)
qvls = v_x_edge
f(i, j) = f(i, j) - 0.5d0*dx*dot_product(qwts, qvls)
cof = [ 0.0d0, 0.0833333333333333d0/dx, -1.25d0/dx, 1.25d0/dx, -0.0833333333333333d0/dx ]
vec = [ v(i-2,j-2), v(i-2,j-1), v(i-2,j+0), v(i-2,j+1), v(i-2,j+2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ v(i-1,j-2), v(i-1,j-1), v(i-1,j+0), v(i-1,j+1), v(i-1,j+2) ]; tmp(-1) = dot_product(cof, vec)
vec = [ v(i+0,j-2), v(i+0,j-1), v(i+0,j+0), v(i+0,j+1), v(i+0,j+2) ]; tmp(+0) = dot_product(cof, vec)
vec = [ v(i+1,j-2), v(i+1,j-1), v(i+1,j+0), v(i+1,j+1), v(i+1,j+2) ]; tmp(+1) = dot_product(cof, vec)
vec = [ v(i+2,j-2), v(i+2,j-1), v(i+2,j+0), v(i+2,j+1), v(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
v_y_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
v_y_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
v_y_edge(2) = dot_product(cof, tmp)
qvls = v_y_edge
f(i, j) = f(i, j) + 0.5d0*dx*dot_product(qwts, qvls)
cof = [ 0.0833333333333333d0/dx, -1.25d0/dx, 1.25d0/dx, -0.0833333333333333d0/dx, 0.0d0 ]
vec = [ v(i-2,j-2), v(i-2,j-1), v(i-2,j+0), v(i-2,j+1), v(i-2,j+2) ]; tmp(-2) = dot_product(cof, vec)
vec = [ v(i-1,j-2), v(i-1,j-1), v(i-1,j+0), v(i-1,j+1), v(i-1,j+2) ]; tmp(-1) = dot_product(cof, vec)
vec = [ v(i+0,j-2), v(i+0,j-1), v(i+0,j+0), v(i+0,j+1), v(i+0,j+2) ]; tmp(+0) = dot_product(cof, vec)
vec = [ v(i+1,j-2), v(i+1,j-1), v(i+1,j+0), v(i+1,j+1), v(i+1,j+2) ]; tmp(+1) = dot_product(cof, vec)
vec = [ v(i+2,j-2), v(i+2,j-1), v(i+2,j+0), v(i+2,j+1), v(i+2,j+2) ]; tmp(+2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
v_y_edge(0) = dot_product(cof, tmp)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
v_y_edge(1) = dot_product(cof, tmp)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
v_y_edge(2) = dot_product(cof, tmp)
qvls = v_y_edge
f(i, j) = f(i, j) - 0.5d0*dx*dot_product(qwts, qvls)
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
qwts = [ 0.555555555555556d0, 0.888888888888889d0, 0.555555555555556d0 ]
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
vec = [ u(i-2,j-2), u(i-1,j-2), u(i+0,j-2), u(i+1,j-2), u(i+2,j-2) ]
tmp(-2,0) = dot_product(cof, vec)
vec = [ u(i-2,j-1), u(i-1,j-1), u(i+0,j-1), u(i+1,j-1), u(i+2,j-1) ]
tmp(-1,0) = dot_product(cof, vec)
vec = [ u(i-2,j+0), u(i-1,j+0), u(i+0,j+0), u(i+1,j+0), u(i+2,j+0) ]
tmp(+0,0) = dot_product(cof, vec)
vec = [ u(i-2,j+1), u(i-1,j+1), u(i+0,j+1), u(i+1,j+1), u(i+2,j+1) ]
tmp(+1,0) = dot_product(cof, vec)
vec = [ u(i-2,j+2), u(i-1,j+2), u(i+0,j+2), u(i+1,j+2), u(i+2,j+2) ]
tmp(+2,0) = dot_product(cof, vec)
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
vec = [ u(i-2,j-2), u(i-1,j-2), u(i+0,j-2), u(i+1,j-2), u(i+2,j-2) ]
tmp(-2,1) = dot_product(cof, vec)
vec = [ u(i-2,j-1), u(i-1,j-1), u(i+0,j-1), u(i+1,j-1), u(i+2,j-1) ]
tmp(-1,1) = dot_product(cof, vec)
vec = [ u(i-2,j+0), u(i-1,j+0), u(i+0,j+0), u(i+1,j+0), u(i+2,j+0) ]
tmp(+0,1) = dot_product(cof, vec)
vec = [ u(i-2,j+1), u(i-1,j+1), u(i+0,j+1), u(i+1,j+1), u(i+2,j+1) ]
tmp(+1,1) = dot_product(cof, vec)
vec = [ u(i-2,j+2), u(i-1,j+2), u(i+0,j+2), u(i+1,j+2), u(i+2,j+2) ]
tmp(+2,1) = dot_product(cof, vec)
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
vec = [ u(i-2,j-2), u(i-1,j-2), u(i+0,j-2), u(i+1,j-2), u(i+2,j-2) ]
tmp(-2,2) = dot_product(cof, vec)
vec = [ u(i-2,j-1), u(i-1,j-1), u(i+0,j-1), u(i+1,j-1), u(i+2,j-1) ]
tmp(-1,2) = dot_product(cof, vec)
vec = [ u(i-2,j+0), u(i-1,j+0), u(i+0,j+0), u(i+1,j+0), u(i+2,j+0) ]
tmp(+0,2) = dot_product(cof, vec)
vec = [ u(i-2,j+1), u(i-1,j+1), u(i+0,j+1), u(i+1,j+1), u(i+2,j+1) ]
tmp(+1,2) = dot_product(cof, vec)
vec = [ u(i-2,j+2), u(i-1,j+2), u(i+0,j+2), u(i+1,j+2), u(i+2,j+2) ]
tmp(+2,2) = dot_product(cof, vec)
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_interior(0,0) = dot_product(cof, tmp(:,0))
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_interior(0,1) = dot_product(cof, tmp(:,0))
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_interior(0,2) = dot_product(cof, tmp(:,0))
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_interior(1,0) = dot_product(cof, tmp(:,1))
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_interior(1,1) = dot_product(cof, tmp(:,1))
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_interior(1,2) = dot_product(cof, tmp(:,1))
cof = [ -0.0147591832102836d0, 0.0837542401323444d0, 1.07043333333333d0, -0.159887573465678d0, 0.0204591832102836d0 ]
u_interior(2,0) = dot_product(cof, tmp(:,2))
cof = [ 0.00468750000000000d0, -0.0604166666666667d0, 1.11145833333333d0, -0.0604166666666667d0, 0.00468750000000000d0 ]
u_interior(2,1) = dot_product(cof, tmp(:,2))
cof = [ 0.0204591832102836d0, -0.159887573465678d0, 1.07043333333333d0, 0.0837542401323443d0, -0.0147591832102836d0 ]
u_interior(2,2) = dot_product(cof, tmp(:,2))
qvls2 = u_interior/(phi*u_interior + 1)
qvls1(0) = 0.5d0*dx*dot_product(qwts, qvls2(:, 0))
qvls1(1) = 0.5d0*dx*dot_product(qwts, qvls2(:, 1))
qvls1(2) = 0.5d0*dx*dot_product(qwts, qvls2(:, 2))
f(i, j) = f(i, j) + 0.5d0*dot_product(qwts, qvls1)/dx
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
