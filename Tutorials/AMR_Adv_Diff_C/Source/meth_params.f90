
! This module stores the runtime parameters.  
! These parameter are initialized in set_method_params().

module meth_params_module

  implicit none

  integer, parameter     :: NHYP    = 4
  integer, parameter     :: MAXADV  = 2

  integer         , save :: NTHERM, NVAR
  integer         , save :: URHO, UX, UY, UZ, UFA
  integer         , save :: nadv

end module meth_params_module
