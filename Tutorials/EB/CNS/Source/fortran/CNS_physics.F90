module cns_physics_module
  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private

  real(rt), parameter, public :: Ru = 8.31451d+07

  real(rt), save, public :: gamma  = 1.4d0     ! gamma-law EOS
  real(rt), save, public :: eos_mu = 28.97d0   ! mean molecular weight

  real(rt), save, public :: Pr    = 0.72d0    ! Prandtl number
!  real(rt), save, public :: Sc    = 0.72d0    ! Schmidt number
  real(rt), save, public :: C_S   = 1.458d-5  ! constant in Sutherland's law
  real(rt), save, public :: T_S   = 110.4d0   ! Sutherland temperature

  real(rt), save, public :: cv                ! Ru/(eos_mu*(gamma-1.d0))
  real(rt), save, public :: cp

  logical, save, public :: use_const_visc = .false.
  real(rt), save, public :: const_visc_mu = -1.d0
  real(rt), save, public :: const_visc_ki = -1.d0
  real(rt), save, public :: const_lambda  = -1.d0

  public :: physics_init

contains

  subroutine physics_init ()
    use amrex_parmparse_module
    type(amrex_parmparse) :: pp

    integer :: i_use_const_visc
    integer :: do_diffusion

    call amrex_parmparse_build(pp,"physics")

    call pp%query("gamma",gamma)
    call pp%query("eos_mu", eos_mu)

    call pp%query("Pr", Pr)
    call pp%query("C_S", C_S)
    call pp%query("T_S", T_S)

    i_use_const_visc = 0
    call pp%query("use_const_visc", i_use_const_visc)
    use_const_visc = i_use_const_visc .ne. 0
    if (use_const_visc) then
       call pp%get("const_visc_mu", const_visc_mu)
       call pp%get("const_visc_ki", const_visc_ki)
       call pp%get("const_lambda", const_lambda)
    end if

    do_diffusion = 1
    call pp%query("do_diffusion", do_diffusion)
    if (do_diffusion .eq. 0) then
       use_const_visc = .true.
       const_visc_mu = 0.d0
       const_visc_ki = 0.d0
       const_lambda = 0.d0
    end if

    call amrex_parmparse_destroy(pp)

    cv = Ru / (eos_mu * (gamma-1.d0))
    cp = gamma * Ru / (eos_mu * (gamma-1.d0))
  end subroutine physics_init

end module cns_physics_module
