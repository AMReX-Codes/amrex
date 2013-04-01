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
  integer, parameter, public ::         SLIP = 16
  integer, parameter, public ::      NO_SLIP = 17

  integer, parameter, public :: REFLECT_ODD  =  20
  integer, parameter, public :: REFLECT_EVEN =  21
  integer, parameter, public :: HOM_DIR      =  22
  integer, parameter, public :: EXT_DIR      =  23
  integer, parameter, public :: FOEXTRAP     =  24
  integer, parameter, public :: HOEXTRAP     =  25
  integer, parameter, public :: DIR_VEL      =  26
  integer, parameter, public :: DIR_TRACT    =  27


contains

  function bc_string_to_integer(str) result (bc_int)

    character (len=*), intent(in) :: str

    integer :: bc_int

    select case (str)

    case ("periodic")
       bc_int = PERIODIC

    case ("inlet")
       bc_int = INLET

    case ("outlet")
       bc_int = OUTLET

    case ("symmetry")
       bc_int = SYMMETRY

    case ("slip wall")
       bc_int = SLIP_WALL

    case ("no slip wall")
       bc_int = NO_SLIP_WALL
      
    case default 
       bc_int = UNDEFINED

    end select

  end function bc_string_to_integer

  function bc_integer_to_string(int) result (str)

    integer, intent(in) :: int

    character (len=20) :: str

    select case (int)

    case (PERIODIC)
       str = "periodic"

    case (INLET)
       str = "inlet"

    case (OUTLET)
       str = "outlet"

    case (SYMMETRY)
       str = "symmetry"

    case (SLIP_WALL)
       str = "slip wall"
    
    case (NO_SLIP_WALL)
       str = "no slip wall"
      
    case default 
       str = "undefined"

    end select

  end function bc_integer_to_string

end module bc_module
