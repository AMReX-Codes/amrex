module cns_physics_module
  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private

  real(rt), parameter, public :: gamma = 1.4d0     ! gamma-law EOS
  real(rt), parameter, public :: mu    = 28.97d0   ! mean molecular weight
  real(rt), parameter, public :: Pr    = 0.72d0    ! Prandtl number
  real(rt), parameter, public :: Sc    = 0.72d0    ! Schmidt number
  real(rt), parameter, public :: C_S   = 1.458d-5  ! constant in Sutherland's law
  real(rt), parameter, public :: T_S   = 110.4d0   ! Sutherland temperature

  real(rt), parameter, public :: Ru = 8.31451d+07

  real(rt), parameter, public :: R_over_mu = Ru/mu
  real(rt), parameter, public :: cv = Ru/(mu*(gamma-1.d0))

end module cns_physics_module
