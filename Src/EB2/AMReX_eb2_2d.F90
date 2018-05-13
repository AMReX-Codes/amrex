
module amrex_eb2_2d_moudle

  use amrex_fort_module
  implicit none

  private
  public :: amrex_eb2_gfab_build_types

contains

  subroutine amrex_eb2_gfab_build_types (lo, hi, s, slo, shi, cell, clo, chi, &
       fx, fxlo, fxhi, fy, fylo, fyhi) bind(c,name='amrex_eb2_gfab_build_types')
    integer, dimension(2), intent(in) :: lo, hi, slo, shi, clo, chi, fxlo, fxhi, fylo, fyhi
    real(amrex_real), intent(in   ) ::    s( slo(1): shi(1), slo(2): shi(2))
    integer(c_int)  , intent(inout) :: cell( clo(1): chi(1), clo(2): chi(2))
    integer(c_int)  , intent(inout) ::   fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    integer(c_int)  , intent(inout) ::   fy(fylo(1):fyhi(1),fylo(2):fyhi(2))

    
  end subroutine amrex_eb2_gfab_build_types

end module amrex_eb2_2d_moudle
