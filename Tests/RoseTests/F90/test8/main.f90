program main

  use init_data_module, only : multifab, init_data

  implicit none

  integer, parameter :: dm = 3
  type(multifab) :: data
  double precision, allocatable :: dx(:)
  
  data%dim = dm
  
  allocate(dx(data%dim))
  dx = 1.d-1
  
  call init_data(data,dx)

end program main
