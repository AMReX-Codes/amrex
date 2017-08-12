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

  ! boundary condition information
  integer, save, public :: physbc_lo(3)
  integer, save, public :: physbc_hi(3)
  integer, save, public :: Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall

  ! problem domain
  real(rt), save, public :: problo(3), probhi(3), center(3)

  public :: cns_init_fort

contains

  subroutine cns_init_fort (physbc_lo_in, physbc_hi_in, &
       Interior_in, Inflow_in, Outflow_in, Symmetry_in, SlipWall_in, NoSlipWall_in, &
       problo_in, probhi_in) &
       bind(c,name='cns_init_fort')
    integer, intent(in) :: physbc_lo_in(3), physbc_hi_in(3)
    integer, value, intent(in) :: Interior_in, Inflow_in, Outflow_in, Symmetry_in, SlipWall_in, NoSlipWall_in
    real(rt), intent(in) :: problo_in(3), probhi_in(3)

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
    center = 0.5_rt*(problo+probhi)

  end subroutine cns_init_fort

end module cns_module
