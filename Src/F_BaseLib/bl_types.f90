!! Defines basic types
!!
!! Most data in _BoxLib_ applications can be represented
!! with default types.  The exception is that for almost all
!! cases we need to use double precision floats.
module bl_types

  implicit none

  integer, parameter :: dp_t = kind(0.0d0)
  integer, parameter :: sp_t = kind(0.0)
  
  integer, parameter :: ll_t = selected_int_kind(15)

end module bl_types
