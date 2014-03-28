module init_data_module

  implicit none

  type multifab
     integer :: dim
  end type multifab

contains

  subroutine init_data(data, dx)
    type(multifab)  , intent(inout) :: data
    double precision, intent(in   ) :: dx(data%dim)

  end subroutine init_data

end module init_data_module
