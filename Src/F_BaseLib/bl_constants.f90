!! Convenient definitions of numerical constants.  This should not
!! be considered an exhaustive list of such constants.
module bl_constants_module

  use bl_types

  implicit none

  real(kind = dp_t), parameter :: ZERO    =  0.0_dp_t
  real(kind = dp_t), parameter :: ONE     =  1.0_dp_t
  real(kind = dp_t), parameter :: TWO     =  2.0_dp_t
  real(kind = dp_t), parameter :: THREE   =  3.0_dp_t
  real(kind = dp_t), parameter :: FOUR    =  4.0_dp_t
  real(kind = dp_t), parameter :: FIVE    =  5.0_dp_t
  real(kind = dp_t), parameter :: SIX     =  6.0_dp_t
  real(kind = dp_t), parameter :: SEVEN   =  7.0_dp_t
  real(kind = dp_t), parameter :: EIGHT   =  8.0_dp_t
  real(kind = dp_t), parameter :: NINE    =  9.0_dp_t
  real(kind = dp_t), parameter :: TEN     = 10.0_dp_t

  real(kind = dp_t), parameter :: ELEVEN  = 11.0_dp_t
  real(kind = dp_t), parameter :: TWELVE  = 12.0_dp_t
  real(kind = dp_t), parameter :: FIFTEEN = 15.0_dp_t
  real(kind = dp_t), parameter :: SIXTEEN = 16.0_dp_t

  real(kind = dp_t), parameter :: HALF    = 0.5_dp_t
  real(kind = dp_t), parameter :: THIRD   = ONE/THREE
  real(kind = dp_t), parameter :: FOURTH  = 0.25_dp_t
  real(kind = dp_t), parameter :: FIFTH   = ONE/FIVE
  real(kind = dp_t), parameter :: SIXTH   = ONE/SIX
  real(kind = dp_t), parameter :: SEVENTH = ONE/SEVEN
  real(kind = dp_t), parameter :: EIGHTH  = 0.125_dp_t
  real(kind = dp_t), parameter :: NINETH  = ONE/NINE
  real(kind = dp_t), parameter :: TENTH   = 0.10_dp_t
  real(kind = dp_t), parameter :: TWELFTH = ONE/TWELVE

  real(kind = dp_t), parameter :: TWO3RD    = TWO/THREE
  real(kind = dp_t), parameter :: FOUR3RD   = FOUR/THREE
  real(kind = dp_t), parameter :: FIVE3RD   = FIVE/THREE

  real(kind = dp_t), parameter :: SEVEN12TH = SEVEN/TWELVE

  !! Pi
  real(kind = dp_t), parameter :: M_PI    = &
       3.141592653589793238462643383279502884197_dp_t
  real(kind = dp_t), parameter :: M_SQRT_PI  = &
       1.772453850905516027298167483341145182798_dp_t

  !! Roots
  real(kind = dp_t), parameter :: M_SQRT_2  = &
       1.414213562373095048801688724209698078570_dp_t

end module bl_constants_module
