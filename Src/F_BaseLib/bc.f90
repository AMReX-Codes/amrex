module bc_module

  implicit none

  integer, parameter, public :: BC_UNDEF     = -HUGE(1)
  integer, parameter, public :: BC_PER       = -1
  integer, parameter, public :: BC_INT       = 0
  integer, parameter, public :: BC_DIR       = 1
  integer, parameter, public :: BC_NEU       = 2
! integer, parameter, public :: BC_ROB       = 3 ! Not ready yet

  integer, parameter, public :: UNDEFINED    = -111

  integer, parameter, public :: PERIODIC     = -1
  integer, parameter, public :: INTERIOR     =  0

  integer, parameter, public :: INLET        = 11
  integer, parameter, public :: OUTLET       = 12
  integer, parameter, public :: SYMMETRY     = 13
  integer, parameter, public ::    SLIP_WALL = 14
  integer, parameter, public :: NO_SLIP_WALL = 15

  integer, parameter, public :: REFLECT_ODD  =  20
  integer, parameter, public :: REFLECT_EVEN =  21
  integer, parameter, public :: FOEXTRAP     =  22
  integer, parameter, public :: EXT_DIR      =  23
  integer, parameter, public :: HOEXTRAP     =  24

end module bc_module
