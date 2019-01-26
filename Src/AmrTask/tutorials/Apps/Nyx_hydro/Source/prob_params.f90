
! This module stores the runtime parameters that define the problem domain.  
! These parameter are initialized in set_problem_params().

module prob_params_module

  implicit none
  integer         , save :: physbc_lo(3)
  integer         , save :: physbc_hi(3)
  integer         , save :: Outflow, Symmetry
  integer         , save :: coord_type

end module prob_params_module
