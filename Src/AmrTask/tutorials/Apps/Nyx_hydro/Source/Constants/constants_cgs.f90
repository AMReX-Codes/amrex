! Fundamental constants taken from NIST's 2010 CODATA recommended values

module fundamental_constants_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  ! Pi
  real(rt), parameter :: pi = 3.141592653589793238d0

  ! newton's gravitational constant
  real(rt), parameter :: Gconst = 6.67428e-8_rt      ! cm^3/g/s^2
  ! new value; if uncommented initial models will need to be re-HSE'ed
  !real(rt), parameter :: Gconst = 6.67384e-8_rt      ! cm^3/g/s^2

  ! boltzmann's constant
  real(rt), parameter :: k_B = 1.3806488e-16_rt   ! erg/K

  ! planck's constant over 2pi
  real(rt), parameter :: hbar = 1.054571726e-27_rt ! erg s

  ! planck's constant
  real(rt), parameter :: hplanck = 6.62606957e-27_rt ! erg s

  ! avogradro's Number
  real(rt), parameter :: n_A = 6.02214129e23_rt   ! mol^-1

  ! Thomson cross section [cm^2]
  real(rt), parameter :: sigma_T =  6.6524587158e-25_rt

  ! convert eV to erg
  real(rt), parameter :: ev2erg = 1.602176487e-12_rt

  ! convert MeV to eV
  real(rt), parameter :: MeV2eV = 1.0e6_rt

  ! mass of proton
  real(rt), parameter :: m_proton = 1.672621777e-24_rt ! g

  ! mass of neutron
  real(rt), parameter :: m_n = 1.674927351e-24_rt ! g

  ! mass of electron
  real(rt), parameter :: m_e = 9.10938291e-28_rt  ! g

  ! speed of light in vacuum
  real(rt), parameter :: c_light = 2.99792458e10_rt   ! cm/s

  real(rt), parameter :: mp_over_kb = m_proton/k_B

  ! electron charge
  ! NIST: q_e = 1.602176565e-19 C
  !
  ! C is the SI unit Coulomb; in cgs we have the definition:
  !     1 C = 0.1 * |c_light| * 1 statC
  ! where statC is the cgs unit statCoulomb; 1 statC = 1 erg^1/2 cm^1/2
  ! and |c_light| is the speed of light in cgs (but without units)
  real(rt), parameter :: q_e     = 4.80320451e-10_rt  ! erg^1/2 cm^1/2

  ! stefan-boltzmann constant
  real(rt), parameter :: sigma_SB = 5.670373e-5_rt    ! erg/s/cm^2/K^4

  ! radiation constant
  real(rt), parameter :: a_rad = 4.0_rt*sigma_SB/c_light

  ! Hubble constant (in s^{-1}, converted from 100 (km/s)/Mpc by dividing by 3.08568025e19km/Mpc)
  real(rt), parameter :: Hubble_const = 32.407764868e-19

  !!! conversion factors

  ! conversion factor for density from cosmo to cgs units.
  real(rt), parameter :: density_to_cgs = 1

  ! converstion factor for "heat" from cgs to cosmo units.
  real(rt), parameter :: heat_from_cgs = 1

  ! conversion factor for energy from cosmo to cgs units.
  real(rt), parameter :: e_to_cgs = 1

  ! For AGN accretion rate
  real(rt), parameter :: eddington_const = 4.0d0*pi * Gconst * m_proton / (sigma_T * c_light)

end module fundamental_constants_module
