
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

     subroutine amrex_fi_fillpatch_two(mf, time, &
          cmf, ctime, nc, fmf, ftime, nf, scomp, dcomp, ncomp, &
          cgeom, fgeom, cbc, fbc, rr, interp, lo_bc, hi_bc) bind(c)
       import
       implicit none
       type(c_ptr), value :: mf, cgeom, fgeom, cbc, fbc
       type(c_ptr), intent(in) :: cmf(*), fmf(*), lo_bc(*), hi_bc(*)
       real(amrex_real), value :: time
       real(amrex_real), intent(in) :: ctime(*), ftime(*)
       integer, value :: nc, nf, scomp, dcomp, ncomp, rr, interp
     end subroutine amrex_fi_fillpatch_two
  end interface

contains

  subroutine amrex_fillpatch_single (mf, told, mfold, tnew, mfnew, geom, fill_physbc, &
       time, scomp, dcomp, ncomp)
    type(amrex_multifab), intent(inout) :: mf
    type(amrex_multifab), intent(in   ) :: mfold, mfnew

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
    call amrex_fi_fillpatch_single(mf%p, time, smf_c, stime, n, &
         scomp-1, dcomp-1, ncomp, geom%p, pbc%p)
  end subroutine amrex_fillpatch_single


  subroutine amrex_fillpatch_two (mf, time, &
       cmf, ctime, fmf, ftime, scomp, dcomp, ncomp, &
       cgeom, fgeom, cpbc, fpbc, rr, interp, lo_bc, hi_bc)
    type(amrex_multifab), intent(inout) :: mf, cmf(:), fmf(:)
    real(amrex_real), intent(in) :: time, ctime(*), ftime(*)
    integer, intent(in) :: scomp, dcomp, ncomp, rr
    type(amrex_geometry), intent(in) :: cgeom, fgeom
    type(amrex_physbc), intent(in) :: cpbc, fpbc
    integer, intent(in) :: interp
    integer, target, intent(in) :: lo_bc(amrex_spacedim,scomp+ncomp-1)
    integer, target, intent(in) :: hi_bc(amrex_spacedim,scomp+ncomp-1)
    integer :: i, nc, nf
    type(c_ptr) :: cmf_c(size(cmf)), fmf_c(size(fmf)), &
         lo_bc_c(scomp+ncomp-1), hi_bc_c(scomp+ncomp-1)
    nc = size(cmf)
    do i = 1, nc
       cmf_c(i) = cmf(i)%p
    end do
    nf = size(fmf)
    do i = 1, nf
       fmf_c(i) = fmf(i)%p
    end do
    do i = 1, scomp-1
       lo_bc_c(i) = c_null_ptr
       hi_bc_c(i) = c_null_ptr
    end do
    do i = scomp, scomp+ncomp-1
       lo_bc_c(i) = c_loc(lo_bc(1,i))
       hi_bc_c(i) = c_loc(hi_bc(1,i))
    end do
    call amrex_fi_fillpatch_two(mf%p, time, &
         cmf_c, ctime, nc, fmf_c, ftime, nf, scomp-1, dcomp-1, ncomp, &
         cgeom%p, fgeom%p, cpbc%p, fpbc%p, rr, interp, lo_bc_c, hi_bc_c)
  end subroutine amrex_fillpatch_two

end module amrex_fillpatch_module
