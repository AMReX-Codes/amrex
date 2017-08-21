module cns_physics_module
  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private

  real(rt), parameter, public :: Ru = 8.31451d+07

  real(rt), save, public :: gamma = 1.4d0     ! gamma-law EOS
  real(rt), save, public :: mu    = 28.97d0   ! mean molecular weight

  real(rt), save, public :: Pr    = 0.72d0    ! Prandtl number
!  real(rt), save, public :: Sc    = 0.72d0    ! Schmidt number
  real(rt), save, public :: C_S   = 1.458d-5  ! constant in Sutherland's law
  real(rt), save, public :: T_S   = 110.4d0   ! Sutherland temperature

  real(rt), save, public :: cv                ! Ru/(mu*(gamma-1.d0))
  real(rt), save, public :: cp

  public :: physics_init

contains

  subroutine physics_init ()
    use amrex_parmparse_module
    type(amrex_parmparse) :: pp

    call amrex_parmparse_build(pp,"physics")
    call pp%query("gamma",gamma)
    call pp%query("mu", mu)
    call amrex_parmparse_destroy(pp)

    cv = Ru / (mu * (gamma-1.d0))
    cp = gamma * Ru / (mu * (gamma-1.d0))
  end subroutine physics_init

end module cns_physics_module
