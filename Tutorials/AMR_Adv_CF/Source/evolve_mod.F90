module evolve_module

  use amrex_amr_module

  use my_amr_module, only : phi_old, phi_new

  implicit none
  private

  public :: evolve

contains

  subroutine evolve ()

  end subroutine evolve

end module evolve_module
