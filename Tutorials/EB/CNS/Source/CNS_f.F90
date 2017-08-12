module cns_module

  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private

  integer, parameter, public :: URHO  = 1
  integer, parameter, public :: UMX   = 2
  integer, parameter, public :: UMY   = 3
  integer, parameter, public :: UMZ   = 4
  integer, parameter, public :: UEDEN = 5
  integer, parameter, public :: UEINT = 6
  integer, parameter, public :: UTEMP = 7
  integer, parameter, public :: NVAR  = 7

  ! physics of the gas
  real(rt), save, public :: gamma
  real(rt), save, public :: mu
  real(rt), save, public :: Pr
  real(rt), save, public :: Sc
  real(rt), save, public :: C_S
  real(rt), save, public :: T_S

  ! boundary condition information
  integer, save, public :: physbc_lo(3)
  integer, save, public :: physbc_hi(3)
  integer, save, public :: Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall

  ! problem domain
  real(rt), save, public :: problo(3), probhi(3)

  public :: cns_init_fort

contains

  subroutine cns_init_fort (gamma_in, mu_in, Pr_in, Sc_in, C_S_in, T_S_in, &
       physbc_lo_in, physbc_hi_in, &
       Interior_in, Inflow_in, Outflow_in, Symmetry_in, SlipWall_in, NoSlipWall_in, &
       problo_in, probhi_in) &
       bind(c,name='cns_init_fort')
    real(rt), value, intent(in) :: gamma_in, mu_in, Pr_in, Sc_in, C_S_in, T_S_in
    integer, intent(in) :: physbc_lo_in(3), physbc_hi_in(3)
    integer, value, intent(in) :: Interior_in, Inflow_in, Outflow_in, Symmetry_in, SlipWall_in, NoSlipWall_in
    real(rt), intent(in) :: problo_in(3), probhi_in(3)
    gamma = gamma_in
    mu    = mu_in
    Pr    = Pr_in
    Sc    = Sc_in
    C_S   = C_S_in
    T_S   = T_S_in

    physbc_lo = physbc_lo_in
    physbc_hi = physbc_hi_in

    Interior   = Interior_in
    Inflow     = Inflow_in
    Outflow    = Outflow_in
    Symmetry   = Symmetry_in
    SlipWall   = SlipWall_in
    NoSlipWall = NoSlipWall_in

    problo = problo_in
    probhi = probhi_in

  end subroutine cns_init_fort

end module cns_module
