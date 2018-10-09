module react_system_module

  implicit none

  integer, parameter :: neqs = 3

contains

  subroutine cv_f_rhs(y, ydot) bind(C, name="cv_f_rhs")

    use amrex_fort_module, only: rt=>amrex_real

    implicit none

    real(rt), intent(in) :: y(neqs)
    real(rt), intent(inout) :: ydot(neqs)

    !$gpu

    ydot(1) = -.04e0*y(1) + 1.e4*y(2)*y(3)
    ydot(3) = 3.e7*y(2)*y(2)
    ydot(2) = -ydot(1)-ydot(3)

  end subroutine cv_f_rhs


  subroutine cv_f_jtv(v, Jv, y, ydot) bind(C, name="cv_f_jtv")

    use amrex_fort_module, only: rt=>amrex_real

    implicit none

    real(rt), intent(in) :: v(neqs), y(neqs), ydot(neqs)
    real(rt), intent(inout) :: Jv(neqs)

    !$gpu

    Jv(1) = -0.04e0*v(1) + 1.e4*y(3)*v(2) + 1.e4*y(2)*v(3)
    Jv(3) = 6.0e7*y(2)*v(2)
    Jv(2) = 0.04e0*v(1) + (-1.e4*y(3)-6.0e7*y(2))*v(2) + (-1.e4*y(2))*v(3)

  end subroutine cv_f_jtv

end module react_system_module
