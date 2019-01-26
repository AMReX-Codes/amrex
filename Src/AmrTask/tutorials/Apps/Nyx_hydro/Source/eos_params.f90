
! This module stores the runtime EOS species IF they are defined to be constants.  
! These parameter are initialized in set_eos_params().

module eos_params_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), save ::  h_species
  real(rt), save :: he_species

end module eos_params_module
