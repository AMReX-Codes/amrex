
subroutine err_c_2_f_4 (i, j, k, n) bind(c)
  implicit none
  integer :: i,j,k,n
end subroutine err_c_2_f_4

subroutine err_c_value_f_pointer (x) bind(c)
  use iso_c_binding
  implicit none
  real(c_double) :: x
end subroutine err_c_value_f_pointer

subroutine err_c_float_f_double (ifab, lo, hi, x) bind(c)
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(inout) :: ifab(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision, intent(out) :: x
end subroutine err_c_float_f_double

subroutine pass_c_real_f_real (x, y, z) bind(c)
  use iso_c_binding
  real(c_float), value :: x
  double precision, value :: y
  real(c_double), value :: z
end subroutine pass_c_real_f_real

module x_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real

  implicit none

  type, bind(c) :: my_c_struct
     double precision :: x(3)
     integer :: i(100)
  end type my_c_struct

contains

  function err_c_ret_int_f_ret_real () result(r)
    double precision :: r
  end function err_c_ret_int_f_ret_real

  function pass_c_ret_int_f_ret_int () result(r)
    integer :: r
  end function pass_c_ret_int_f_ret_int

  subroutine err_c_poiner_f_value (x) bind(c)
    use amrex_fort_module, only : amrex_real
    real(amrex_real), intent(in), value :: x
  end subroutine err_c_poiner_f_value

  subroutine err_c_reference_f_value (x) bind(c)
    use amrex_fort_module, only : amrex_real
    real(amrex_real), intent(in), value :: x
  end subroutine err_c_reference_f_value

  subroutine pass_c_bl_fort_fab_arg (x, xlo1, xlo2, xlo3, xhi1, xhi2, xhi3) bind(c)
    integer, intent(in) :: xlo1, xlo2, xlo3, xhi1, xhi2, xhi3
    real(amrex_real), intent(inout) :: x(xlo1:xhi1,xlo2:xhi2,xlo3:xhi3)
  end subroutine pass_c_bl_fort_fab_arg

  subroutine pass_c_bl_fort_fab_arg_3d (x, lo, hi) bind(c)
    integer, intent(in) :: lo(3), hi(3)
    real(amrex_real), intent(inout) :: x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  end subroutine pass_c_bl_fort_fab_arg_3d

  subroutine pass_c_pointer_f_c_ptr (p) bind(c)
    use iso_c_binding
    type(c_ptr), value :: p
  end subroutine pass_c_pointer_f_c_ptr

  subroutine err_c_pointer_f_c_ptr (p) bind(c)
    use iso_c_binding
    type(c_ptr) :: p
  end subroutine err_c_pointer_f_c_ptr

  subroutine pass_c_void_pointer_f_any (x, y ,z) bind(c)
    integer :: x
    double precision :: y
    type(my_c_struct) :: z
  end subroutine pass_c_void_pointer_f_any

end module x_module

