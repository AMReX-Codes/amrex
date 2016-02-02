! make time available in a module

module time_module

  use bl_types

  implicit none

  real(dp_t), save :: time

  private

  public :: time

end module time_module
