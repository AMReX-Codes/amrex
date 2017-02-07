module amrex_multifabutil_module

  use iso_c_binding
  use amrex_multifab_module
  use amrex_geometry_module

  implicit none
  private

  public :: amrex_average_down

  interface
     subroutine amrex_fi_average_down (fmf, cmf, fgeom, cgeom, scomp, ncomp, rr) bind(c)
       import
       implicit none
       type(c_ptr), value :: fmf, cmf, fgeom, cgeom
       integer(c_int), value :: scomp, ncomp, rr
     end subroutine amrex_fi_average_down
  end interface

contains

  subroutine amrex_average_down (fmf, cmf, fgeom, cgeom, scomp, ncomp, rr)
    type(amrex_multifab), intent(in   ) :: fmf
    type(amrex_multifab), intent(inout) :: cmf
    type(amrex_geometry), intent(in) :: fgeom, cgeom
    integer, intent(in) :: scomp, ncomp, rr
    call amrex_fi_average_down(fmf%p, cmf%p, fgeom%p, cgeom%p, scomp-1, ncomp, rr)
  end subroutine amrex_average_down

end module amrex_multifabutil_module
