
! This module stores the runtime EOS species IF they are defined to be constants.  
! These parameter are initialized in set_eos_params().

module eos_params_module

  use constants_module, only : rt => type_real, M_PI

  implicit none

  real(rt), save ::  h_species
  real(rt), save :: he_species

end module eos_params_module
