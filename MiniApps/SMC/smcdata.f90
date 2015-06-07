module smcdata_module

  use multifab_module

  implicit none

  type(multifab),save :: Uprime, Unew
  type(multifab),save :: Q
  type(multifab),save :: mu, xi ! viscosity
  type(multifab),save :: lam ! partial thermal conductivity
  type(multifab),save :: Ddiag ! diagonal components of rho * Y_k * D

  private

  public :: Uprime, Unew, Q, mu, xi, lam, Ddiag
  public :: build_smcdata, destroy_smcdata

contains

  subroutine build_smcdata(la)

    use variables_module, only : ncons, nprim
    use chemistry_module, only : nspecies
    use derivative_stencil_module, only : stencil_ng

    implicit none

    type(layout), intent(in) :: la

    call multifab_build(Uprime, la, ncons, 0)
    call multifab_build(Unew,   la, ncons, stencil_ng)

    call multifab_build(Q, la, nprim, stencil_ng)
    call multifab_setval(Q, 0.d0, .true.)

    call multifab_build(mu , la, 1, stencil_ng)
    call multifab_build(xi , la, 1, stencil_ng)
    call multifab_build(lam, la, 1, stencil_ng)
    call multifab_build(Ddiag, la, nspecies, stencil_ng)

  end subroutine build_smcdata


  subroutine destroy_smcdata()
    call destroy(Unew)
    call destroy(Uprime)
    call destroy(Q)
    call destroy(mu)
    call destroy(xi)
    call destroy(lam)
    call destroy(Ddiag)
  end subroutine destroy_smcdata

end module smcdata_module

