module chemotaxis_kernels
  !
  ! These kernels use the 4th order stencils of KKM.
  !

  implicit none
contains


  !
  ! Cell motility: diff * grad(u)
  !
  subroutine cell_motility(f, u, lo, hi, ng, dx, diff) bind(c, name='cell_motility')
    use iso_c_binding
    integer(c_int), intent(in)        :: lo(2), hi(2)
    integer(c_int), intent(in), value :: ng
    real(c_double), intent(in), value :: dx, diff
    real(c_double), intent(in)        :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(inout)     :: f(lo(1):hi(1),lo(2):hi(2))

    integer :: i, j

    ! this kernel was generated

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          f(i, j) = f(i, j) + diff * ( &
             0.0833333333333333d0*(-15.0d0*u(i+0,j+0) + 15.0d0*u(i+0,j+1) - u(i+0,j+2) + u(i+0,j-1)) &
           + 0.0833333333333333d0*(-15.0d0*u(i+0,j+0) + 15.0d0*u(i+1,j+0) - u(i+2,j+0) + u(i-1,j+0)) &
           - 0.0833333333333333d0*( 15.0d0*u(i+0,j+0) - u(i+0,j+1) - 15.0d0*u(i+0,j-1) + u(i+0,j-2)) &
           - 0.0833333333333333d0*( 15.0d0*u(i+0,j+0) - u(i+1,j+0) - 15.0d0*u(i-1,j+0) + u(i-2,j+0)) ) / dx
       end do
    end do

  end subroutine cell_motility


  !
  ! Chemotactic sensitivity: chi * u * grad(v)
  !
  subroutine chemotactic_sensitivity(f, u, v, lo, hi, ng, dx, chi) bind(c, name='chemotactic_sensitivity')
    use iso_c_binding
    integer(c_int), intent(in)        :: lo(2), hi(2)
    integer(c_int), intent(in), value :: ng
    real(c_double), intent(in), value :: dx, chi
    real(c_double), intent(in)        :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(in)        :: v(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(inout)     :: f(lo(1):hi(1),lo(2):hi(2))

    real(c_double), dimension(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng) :: phi_x, phi_y, rho_x, rho_y

    integer :: i, j

    ! this kernel was generated

    do j = lo(2)-ng/2, hi(2)+ng/2
       do i = lo(1)-ng/2, hi(1)+ng/2
          rho_y(i,j) = 0.583333333333333d0*u(i+0,j+0) + 0.583333333333333d0*u(i+0,j+1) - &
               0.0833333333333333d0*u(i+0,j+2) - 0.0833333333333333d0*u(i+0,j-1)
          phi_x(i,j) = -1.25d0*v(i+0,j+0)/dx + 1.25d0*v(i+1,j+0)/dx - 0.0833333333333333d0*v( &
               i+2,j+0)/dx + 0.0833333333333333d0*v(i-1,j+0)/dx
          phi_y(i,j) = -1.25d0*v(i+0,j+0)/dx + 1.25d0*v(i+0,j+1)/dx - 0.0833333333333333d0*v( &
               i+0,j+2)/dx + 0.0833333333333333d0*v(i+0,j-1)/dx
          rho_x(i,j) = 0.583333333333333d0*u(i+0,j+0) + 0.583333333333333d0*u(i+1,j+0) - &
               0.0833333333333333d0*u(i+2,j+0) - 0.0833333333333333d0*u(i-1,j+0)
       end do
    end do

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          f(i, j) = f(i, j) - chi * (phi_x(i+0,j+0)*rho_x(i+0,j+0) - phi_x(i-1,j+0)*rho_x(i-1,j+0) + phi_y( &
               i+0,j+0)*rho_y(i+0,j+0) - phi_y(i+0,j-1)*rho_y(i+0,j-1) + &
               3.61689814814815d-5*(34.0d0*phi_x(i+0,j+1) - 5.0d0*phi_x(i+0,j+2 &
               ) - 34.0d0*phi_x(i+0,j-1) + 5.0d0*phi_x(i+0,j-2))*(34.0d0*rho_x( &
               i+0,j+1) - 5.0d0*rho_x(i+0,j+2) - 34.0d0*rho_x(i+0,j-1) + 5.0d0* &
               rho_x(i+0,j-2)) - 3.61689814814815d-5*(34.0d0*phi_x(i-1,j+1) - &
               5.0d0*phi_x(i-1,j+2) - 34.0d0*phi_x(i-1,j-1) + 5.0d0*phi_x(i-1, &
               j-2))*(34.0d0*rho_x(i-1,j+1) - 5.0d0*rho_x(i-1,j+2) - 34.0d0* &
               rho_x(i-1,j-1) + 5.0d0*rho_x(i-1,j-2)) + 3.61689814814815d-5*( &
               34.0d0*phi_y(i+1,j+0) - 5.0d0*phi_y(i+2,j+0) - 34.0d0*phi_y(i-1, &
               j+0) + 5.0d0*phi_y(i-2,j+0))*(34.0d0*rho_y(i+1,j+0) - 5.0d0*rho_y &
               (i+2,j+0) - 34.0d0*rho_y(i-1,j+0) + 5.0d0*rho_y(i-2,j+0)) - &
               3.61689814814815d-5*(34.0d0*phi_y(i+1,j-1) - 5.0d0*phi_y(i+2,j-1 &
               ) - 34.0d0*phi_y(i-1,j-1) + 5.0d0*phi_y(i-2,j-1))*(34.0d0*rho_y( &
               i+1,j-1) - 5.0d0*rho_y(i+2,j-1) - 34.0d0*rho_y(i-1,j-1) + 5.0d0* &
               rho_y(i-2,j-1)))
       end do
    end do

  end subroutine chemotactic_sensitivity


  !
  ! Signal diffusion: grad(v)
  !
  subroutine signal_diffusion(f, v, lo, hi, ng, dx) bind(c, name='signal_diffusion')
    use iso_c_binding
    integer(c_int), intent(in)        :: lo(2), hi(2)
    integer(c_int), intent(in), value :: ng
    real(c_double), intent(in), value :: dx
    real(c_double), intent(in)        :: v(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(inout)     :: f(lo(1):hi(1),lo(2):hi(2))

    integer :: i, j

    ! this kernel was generated

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          f(i, j) = f(i, j) + ( &
             0.0833333333333333d0*(-15.0d0*v(i+0,j+0) + 15.0d0*v(i+0,j+1) - v(i+0,j+2) + v(i+0,j-1))/dx &
           + 0.0833333333333333d0*(-15.0d0*v(i+0,j+0) + 15.0d0*v(i+1,j+0) - v(i+2,j+0) + v(i-1,j+0))/dx &
           - 0.0833333333333333d0*( 15.0d0*v(i+0,j+0) - v(i+0,j+1) - 15.0d0*v(i+0,j-1) + v(i+0,j-2))/dx &
           - 0.0833333333333333d0*( 15.0d0*v(i+0,j+0) - v(i+1,j+0) - 15.0d0*v(i-1,j+0) + v(i-2,j+0))/dx )
       end do
    end do

  end subroutine signal_diffusion


  !
  ! Signal production: u
  !
  subroutine signal_production(f, u, v, lo, hi, ng, dx) bind(c, name='signal_production')
    use iso_c_binding
    integer(c_int), intent(in)        :: lo(2), hi(2)
    integer(c_int), intent(in), value :: ng
    real(c_double), intent(in), value :: dx
    real(c_double), intent(in)        :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(in)        :: v(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(inout)     :: f(lo(1):hi(1),lo(2):hi(2))

    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          f(i, j) = f(i, j) + u(i, j)
       end do
    end do

  end subroutine signal_production


  !
  ! Signal degredation: -v
  !
  subroutine signal_degradation(f, u, v, lo, hi, ng, dx) bind(c, name='signal_degradation')
    use iso_c_binding
    integer(c_int), intent(in)        :: lo(2), hi(2)
    integer(c_int), intent(in), value :: ng
    real(c_double), intent(in), value :: dx
    real(c_double), intent(in)        :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(in)        :: v(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(c_double), intent(inout)     :: f(lo(1):hi(1),lo(2):hi(2))

    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          f(i, j) = f(i, j) - v(i, j)
       end do
    end do

  end subroutine signal_degradation

end module chemotaxis_kernels
