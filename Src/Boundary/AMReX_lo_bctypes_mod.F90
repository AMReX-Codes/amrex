#include "AMReX_LO_BCTYPES.H"

module amrex_lo_bctypes_module

  implicit none
  integer, parameter :: amrex_lo_dirichlet         = AMREX_LO_DIRICHLET
  integer, parameter :: amrex_lo_neumann           = AMREX_LO_NEUMANN
  integer, parameter :: amrex_lo_reflect_odd       = AMREX_LO_REFLECT_ODD
  integer, parameter :: amrex_lo_marshak           = AMREX_LO_MARSHAK
  integer, parameter :: amrex_lo_sanchez_pomraning = AMREX_LO_SANCHEZ_POMRANING
  integer, parameter :: amrex_lo_inflow            = AMREX_LO_INFLOW
  integer, parameter :: amrex_lo_periodic          = AMREX_LO_PERIODIC
  integer, parameter :: amrex_lo_bogus             = AMREX_LO_BOGUS

end module amrex_lo_bctypes_module
