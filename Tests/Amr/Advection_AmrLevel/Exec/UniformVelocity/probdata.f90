module probdata_module

  use amrex_fort_module
  implicit none

  real(amrex_real), save, bind(c) :: adv_vel(3)

end module probdata_module
