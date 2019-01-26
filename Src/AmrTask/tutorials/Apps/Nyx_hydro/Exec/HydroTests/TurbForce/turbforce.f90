! This module stores the coefficients for the turbulent forcing

module turbforce_module

  use bl_types
! use bl_space

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,save  :: nmodes, nxmodes, nymodes, nzmodes, mode_start
  integer,save  :: spectrum_type, forcing_type, use_rho_in_forcing
  integer,save  :: do_mode_division

  real(rt),save :: turb_scale
  real(rt),save ::  forcing_xlength, forcing_ylength, forcing_zlength, forcing_epsilon
  real(rt),save ::  forcing_time_scale_min, forcing_time_scale_max, force_scale
  real(rt),save ::  stop_forcing
  real(rt),save ::  AXY(32), BXY(32), CXY(32), DXY(32), PXY(32), QXY(32), RXY(32)
  real(rt),save ::  AZX(32), BZX(32), CZX(32), DZX(32), PZX(32), QZX(32), RZX(32)
  real(rt),save ::  AYZ(32), BYZ(32), CYZ(32), DYZ(32), PYZ(32), QYZ(32), RYZ(32)
  real(rt),save ::  FTX(0:32,0:32,0:32), FTY(0:32,0:32,0:32), FTZ(0:32,0:32,0:32) 
  real(rt),save ::  TAT(0:32,0:32,0:32), TAP(0:32,0:32,0:32)
  real(rt),save ::  FPXX(0:32,0:32,0:32), FPYX(0:32,0:32,0:32), FPZX(0:32,0:32,0:32)
  real(rt),save ::  FPXY(0:32,0:32,0:32), FPYY(0:32,0:32,0:32), FPZY(0:32,0:32,0:32)
  real(rt),save ::  FPXZ(0:32,0:32,0:32), FPYZ(0:32,0:32,0:32), FPZZ(0:32,0:32,0:32)
  real(rt),save ::  FAX(0:32,0:32,0:32), FAY(0:32,0:32,0:32), FAZ(0:32,0:32,0:32)

end module turbforce_module
