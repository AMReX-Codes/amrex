module nbrs_test_module

  use amrex_fort_module, only : amrex_real
  implicit none

  type, bind(c) :: nbr_sten
     real(amrex_real) :: val(-1:1,-1:1,-1:1)
     integer :: iv(0:2)
  end type nbr_sten

end module nbrs_test_module
