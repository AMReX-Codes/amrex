!! Defines basic types
!!
!! Most data in _BoxLib_ applications can be represented
!! with default types.  The exception is that for almost all
!! cases we need to use double precision floats.

module bl_types

  implicit none

  !! Double precision floating point data are declared as
  !! REAL(KIND=dp_t), this is akin to the old non-standard 
  !! REAL*8, or DOUBLE PRECISION

  integer, parameter :: dp_t = kind(0.0d0)

  !! Single precision floating point data are declared as
  !! REAL(kind=sp_t), or simply as REAL.

  integer, parameter :: sp_t = kind(0.0)
  
  !! A long long data type used to attempt to get sufficient
  !! decimal digits to compute box volues.

  integer, parameter :: ll_t = selected_int_kind(15)

end module bl_types
