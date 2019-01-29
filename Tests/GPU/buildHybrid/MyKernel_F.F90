module my_kernel_module
  use iso_c_binding
  implicit none

contains

  attributes(global) subroutine plusone_fort (num) &
       bind(c,name='plusone_fort')
    integer, intent(inout) :: num 

    num = num + 1

  end subroutine plusone_fort

end module my_kernel_module
