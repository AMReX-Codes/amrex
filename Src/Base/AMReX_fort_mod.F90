module amrex_fort_module

  use iso_c_binding, only : c_float, c_double

  implicit none

  integer, parameter ::    bl_spacedim = BL_SPACEDIM
  integer, parameter :: amrex_spacedim = BL_SPACEDIM

#ifdef BL_USE_FLOAT
  integer, parameter :: amrex_real = c_float
#else
  integer, parameter :: amrex_real = c_double
#endif

end module amrex_fort_module
