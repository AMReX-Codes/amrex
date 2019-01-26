module reion_aux_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! Global variables (re)set on inputs
  real(rt), save :: zhi_flash=-1.0, zheii_flash=-1.0, T_zhi=0.0, T_zheii=0.0
  logical, save  :: flash_h=.false., flash_he=.false., inhomogeneous_on=.false.

end module reion_aux_module
