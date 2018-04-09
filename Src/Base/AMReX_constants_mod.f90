
module amrex_constants_module
  use amrex_fort_module, only : amrex_real
  implicit none

  real(kind = amrex_real), parameter :: ZERO    =  0.0_amrex_real
  real(kind = amrex_real), parameter :: ONE     =  1.0_amrex_real
  real(kind = amrex_real), parameter :: TWO     =  2.0_amrex_real
  real(kind = amrex_real), parameter :: THREE   =  3.0_amrex_real
  real(kind = amrex_real), parameter :: FOUR    =  4.0_amrex_real
  real(kind = amrex_real), parameter :: FIVE    =  5.0_amrex_real
  real(kind = amrex_real), parameter :: SIX     =  6.0_amrex_real
  real(kind = amrex_real), parameter :: SEVEN   =  7.0_amrex_real
  real(kind = amrex_real), parameter :: EIGHT   =  8.0_amrex_real
  real(kind = amrex_real), parameter :: NINE    =  9.0_amrex_real
  real(kind = amrex_real), parameter :: TEN     = 10.0_amrex_real

  real(kind = amrex_real), parameter :: ELEVEN  = 11.0_amrex_real
  real(kind = amrex_real), parameter :: TWELVE  = 12.0_amrex_real
  real(kind = amrex_real), parameter :: FIFTEEN = 15.0_amrex_real
  real(kind = amrex_real), parameter :: SIXTEEN = 16.0_amrex_real

  real(kind = amrex_real), parameter :: HALF    = 0.5_amrex_real
  real(kind = amrex_real), parameter :: THIRD   = ONE/THREE
  real(kind = amrex_real), parameter :: FOURTH  = 0.25_amrex_real
  real(kind = amrex_real), parameter :: FIFTH   = ONE/FIVE
  real(kind = amrex_real), parameter :: SIXTH   = ONE/SIX
  real(kind = amrex_real), parameter :: SEVENTH = ONE/SEVEN
  real(kind = amrex_real), parameter :: EIGHTH  = 0.125_amrex_real
  real(kind = amrex_real), parameter :: NINETH  = ONE/NINE
  real(kind = amrex_real), parameter :: TENTH   = 0.10_amrex_real
  real(kind = amrex_real), parameter :: TWELFTH = ONE/TWELVE

  real(kind = amrex_real), parameter :: TWO3RD    = TWO/THREE
  real(kind = amrex_real), parameter :: FOUR3RD   = FOUR/THREE
  real(kind = amrex_real), parameter :: FIVE3RD   = FIVE/THREE
  real(kind = amrex_real), parameter :: FIVE6TH   = FIVE/SIX

  real(kind = amrex_real), parameter :: THREE4TH  = 0.75_amrex_real

  real(kind = amrex_real), parameter :: FIVE12TH = FIVE/TWELVE
  real(kind = amrex_real), parameter :: SEVEN12TH = SEVEN/TWELVE

  real(kind = amrex_real), parameter :: FIVE32ND = FIVE/32.0_amrex_real

  !! Pi
  real(kind = amrex_real), parameter :: M_PI    = &
       3.141592653589793238462643383279502884197_amrex_real
  real(kind = amrex_real), parameter :: M_SQRT_PI  = &
       1.772453850905516027298167483341145182798_amrex_real

  !! Roots
  real(kind = amrex_real), parameter :: M_SQRT_2  = &
       1.414213562373095048801688724209698078570_amrex_real

end module amrex_constants_module
