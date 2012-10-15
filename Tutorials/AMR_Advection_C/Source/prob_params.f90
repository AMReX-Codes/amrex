
! This module stores the runtime parameters that define the problem domain.  
! These parameter are initialized in set_problem_params().

module prob_params_module

  implicit none
  integer         , save, allocatable :: physbc_lo(:)
  integer         , save, allocatable :: physbc_hi(:)
  integer         , save :: Outflow, Symmetry
  integer         , save :: coord_type

end module prob_params_module
