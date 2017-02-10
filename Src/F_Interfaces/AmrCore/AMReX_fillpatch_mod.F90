
module amrex_fillpatch_module

  use amrex_base_module

  implicit none
  private

  public :: amrex_fillpatch

  interface amrex_fillpatch
     module procedure amrex_fillpatch_single
     module procedure amrex_fillpatch_two
  end interface amrex_fillpatch

  interface
     subroutine amrex_fi_fillpatch_single(mf, time, smf, stime, ns, scomp, dcomp, ncomp, &
          geom, pbc) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf, geom, pbc
       type(c_ptr), intent(in) :: smf(*)
       real(amrex_real), value :: time
       real(amrex_real), intent(in) :: stime(*)
       integer(c_int), value :: scomp, dcomp, ncomp, ns
     end subroutine amrex_fi_fillpatch_single
  end interface

  type amrex_interpolater
  end type amrex_interpolater

  type amrex_bcrec
  end type amrex_bcrec

contains

  subroutine amrex_fillpatch_single (mf, time, smf, stime, scomp, dcomp, ncomp, geom, pbc)
    type(amrex_multifab), intent(inout) :: mf, smf(:)
    real(amrex_real) :: time, stime(*)
    integer, intent(in) :: scomp, dcomp, ncomp
    type(amrex_geometry), intent(in) :: geom
    type(amrex_physbc), intent(in) :: pbc

    integer :: i, n
    type(c_ptr) :: smf_c(size(smf))
    
    n = size(smf)
    do i = 1, n
       smf_c(i) = smf(i)%p
    end do
    call amrex_fi_fillpatch_single(mf%p, time, smf_c, stime, n, scomp, dcomp, ncomp, geom%p, pbc%p)
  end subroutine amrex_fillpatch_single


  subroutine amrex_fillpatch_two (mf, time, cmf, ctime, fmf, ftime, scomp, dcomp, ncomp, &
       cgeom, fgeom, cpbc, fpbc, rr, interp, bcs)
    type(amrex_multifab), intent(inout) :: mf, cmf(:), fmf(:)
    real(amrex_real), intent(in) :: time, ctime(*), ftime(*)
    integer, intent(in) :: scomp, dcomp, ncomp, rr
    type(amrex_geometry), intent(in) :: cgeom, fgeom
    type(amrex_physbc), intent(in) :: cpbc, fpbc
    type(amrex_interpolater), intent(in) :: interp
    type(amrex_bcrec), intent(in) :: bcs(*)

  end subroutine amrex_fillpatch_two

end module amrex_fillpatch_module
