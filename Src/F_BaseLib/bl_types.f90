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

  ! integer, parameter :: dp_t = kind(0.0d0)
  integer, parameter, public :: dp_t = selected_real_kind(15,307)

  !! Single precision floating point data are declared as
  !! REAL(kind=sp_t), or simply as REAL.

  ! integer, parameter :: sp_t = kind(0.0)
  integer, parameter, public :: sp_t = selected_real_kind(6, 37)

  !! Quad precision floating point data
  !! This is tricky, by definition a selected_real_kind returns -1 if 
  !! type with required precision and range can't be found, so if we don't
  !! have extended precision reals on this processor, just use dp
  integer, parameter, private :: qp_t_0 = selected_real_kind(30,1000)
  integer, parameter, public :: qp_t = (1 + sign(1,qp_t_0))/2*qp_t_0 + (1 - sign(1,qp_t_0))/2*dp_t

  !! A long long data type used to attempt to get sufficient
  !! decimal digits to compute box volues.
  integer, parameter, private :: i4_t = selected_int_kind(9)
  integer, parameter, private :: i8_t = selected_int_kind(15)
  ! integer, parameter, public :: ll_t = selected_int_kind(15)
  integer, parameter, public :: ll_t = (1 + sign(1,i8_t))/2*i8_t + (1 - sign(1,i8_t))/2*i4_t

  !! Prints useful type information for some of the basic type
  public :: bl_types_info
   
!  logical, private, parameter :: bigendian = IACHAR(TRANSFER(1,"a")) == 0 

contains
  
  subroutine bl_types_info(unit)
    integer, intent(in) :: unit
    write(unit=unit,fmt="(A)") "BL_TYPES_INFO"
    write(unit=unit,fmt="(2X,A)") "Real Types"

    write(unit=unit,fmt="(2x,A)") "  NAME   KIND RANGE PRECISION  DIGITS"
    write(unit=unit,fmt=100) &
         "  DP_T = ", dp_t, range(0.0_dp_t), precision(0.0_dp_t), digits(0.0_dp_t)
    write(unit=unit,fmt=100) &
         "  SP_T = ", sp_t, range(0.0_sp_t), precision(0.0_sp_t), digits(0.0_sp_t)
    if ( qp_t_0 > 0 ) then
       write(unit=unit,fmt="(2x,A,I4)") "QP_T_0 = ", qp_t_0
    else
       write(unit=unit,fmt="(2x,A)") &
            "QP_T_0 = Not Defined: see bl_types.f90 for possible adjustments"
    end if
    write(unit=unit,fmt=100) &
         "QP_T   = ", qp_t, range(0.0_qp_t), precision(0.0_qp_t), digits(0.0_qp_t)

    write(unit=unit,fmt="(2X,A)") "Integer Types"
    write(unit=unit,fmt="(2x,A)")    "  NAME   KIND RANGE  DIGITS"
    write(unit=unit,fmt=200) "  i4_T = ", &
         i4_t, range(0_i4_t), digits(0_i4_t)
    if ( i8_t > 0 ) then
       write(unit=unit,fmt="(2x,A,i4)") "  i8_T = ", i8_t
    else
       write(unit=unit,fmt="(2x,A)") &
            "  i8_T = Not Defined: see bl_types.f90 for possible adjustments"
    end if
    write(unit=unit,fmt=200) "  LL_T = ", ll_t, range(0_ll_t), digits(0_ll_t)
100 format(2x,A,I4,2X,I4,I10, I8)
200 format(2x,A,I4,2X,I4,I8)
  end subroutine bl_types_info

!  function bl_is_bigendian() result(r)
!    logical :: r
!    r = bigendian
!  end function bl_is_bigendian

end module bl_types
