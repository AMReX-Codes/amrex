module constants_module

  use iso_c_binding, only : c_float, c_double, c_size_t

  implicit none
    
  integer, parameter :: type_real = c_double
  ! We could/should use Fortran 2008 c_sizeof here.
  integer (kind=c_size_t), parameter :: type_real_size = 8_c_size_t
  
  real(kind = type_real), parameter :: M_PI    = &
       3.141592653589793238462643383279502884197_type_real
  
 end module constants_module
