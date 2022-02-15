
module amrex_abeclaplacian_module

  use amrex_base_module
  use amrex_linop_module
  implicit none

  private
  public :: amrex_abeclaplacian_build, amrex_abeclaplacian_destroy

  type, extends(amrex_linop), public :: amrex_abeclaplacian
   contains
     generic :: assignment(=) => amrex_abeclaplacian_assign   ! shallow copy
     procedure :: set_scalars => amrex_abeclaplacian_set_scalars
     procedure :: set_acoeffs => amrex_abeclaplacian_set_acoeffs
     procedure :: set_bcoeffs => amrex_abeclaplacian_set_bcoeffs
     procedure, private :: amrex_abeclaplacian_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_abeclaplacian_destroy
#endif
  end type amrex_abeclaplacian

  interface
     subroutine amrex_fi_new_abeclaplacian (abeclap,n,gm,ba,dm,mt,agg,con,mcl) bind(c)
       import
       implicit none
       type(c_ptr), intent(inout) :: abeclap
       integer(c_int), value :: n, mt, agg, con, mcl
       type(c_ptr), intent(in) :: gm(*), ba(*), dm(*)
     end subroutine amrex_fi_new_abeclaplacian

     subroutine amrex_fi_delete_linop (linop) bind(c)
       import
       implicit none
       type(c_ptr), value :: linop
     end subroutine amrex_fi_delete_linop

     subroutine amrex_fi_abeclap_set_scalars (abeclap,a,b) bind(c)
       import
       implicit none
       type(c_ptr), value :: abeclap
       real(amrex_real), value :: a, b
     end subroutine amrex_fi_abeclap_set_scalars

     subroutine amrex_fi_abeclap_set_acoeffs (abeclap, amrlev, alpha) bind(c)
       import
       implicit none
       type(c_ptr), value :: abeclap, alpha
       integer(c_int), value :: amrlev
     end subroutine amrex_fi_abeclap_set_acoeffs

     subroutine amrex_fi_abeclap_set_bcoeffs (abeclap, amrlev, beta) bind(c)
       import
       implicit none
       type(c_ptr), value :: abeclap
       integer(c_int), value :: amrlev
       type(c_ptr), intent(in) :: beta(*)
     end subroutine amrex_fi_abeclap_set_bcoeffs
  end interface

contains

  subroutine amrex_abeclaplacian_assign (dst, src)
    class(amrex_abeclaplacian), intent(inout) :: dst
    type (amrex_abeclaplacian), intent(in   ) :: src
    call amrex_abeclaplacian_destroy(dst)
    dst%owner = .false.
    dst%p = src%p
  end subroutine amrex_abeclaplacian_assign


  subroutine amrex_abeclaplacian_build (abeclap, geom, ba, dm, metric_term, &
       agglomeration, consolidation, max_coarsening_level)
    type(amrex_abeclaplacian), intent(inout) :: abeclap
    type(amrex_geometry), intent(in) :: geom(0:)
    type(amrex_boxarray), intent(in) :: ba(0:)
    type(amrex_distromap), intent(in) :: dm(0:)
    logical, intent(in), optional :: metric_term, agglomeration, consolidation
    integer, intent(in), optional :: max_coarsening_level

    integer(c_int) :: imt, iagg, icon, imcl
    integer(c_int) :: ilev, nlevs
    type(c_ptr) :: gm_cp(0:size(geom)-1)
    type(c_ptr) :: ba_cp(0:size(geom)-1)
    type(c_ptr) :: dm_cp(0:size(geom)-1)

    imt = -1
    if (present(metric_term)) then
       if (metric_term) then
          imt = 1
       else
          imt = 0
       end if
    end if

    iagg = -1
    if (present(agglomeration)) then
       if (agglomeration) then
          iagg = 1
       else
          iagg = 0
       end if
    end if

    icon = -1
    if (present(consolidation)) then
       if (consolidation) then
          icon = 1
       else
          icon = 0
       end if
    end if

    imcl = 30
    if (present(max_coarsening_level)) then
       imcl = max_coarsening_level
    end if

    nlevs = size(geom)
    do ilev = 0, nlevs-1
       gm_cp(ilev) = geom(ilev)%p
       ba_cp(ilev) = ba(ilev)%p
       dm_cp(ilev) = dm(ilev)%p
    end do

    abeclap%owner = .true.
    call amrex_fi_new_abeclaplacian(abeclap%p, nlevs, gm_cp, ba_cp, dm_cp, imt, iagg, icon, imcl)
  end subroutine amrex_abeclaplacian_build


  subroutine amrex_abeclaplacian_destroy (this)
    type(amrex_abeclaplacian), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_linop(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
  end subroutine amrex_abeclaplacian_destroy


  subroutine amrex_abeclaplacian_set_scalars (this, ascalar, bscalar)
    class(amrex_abeclaplacian), intent(inout) :: this
    real(amrex_real) :: ascalar, bscalar
    call amrex_fi_abeclap_set_scalars(this%p, ascalar, bscalar)
  end subroutine amrex_abeclaplacian_set_scalars


  subroutine amrex_abeclaplacian_set_acoeffs (this, lev, acoef)
    class(amrex_abeclaplacian), intent(inout) :: this
    integer, intent(in) :: lev
    type(amrex_multifab), intent(in) :: acoef
    call amrex_fi_abeclap_set_acoeffs(this%p, lev, acoef%p)
  end subroutine amrex_abeclaplacian_set_acoeffs


  subroutine amrex_abeclaplacian_set_bcoeffs (this, lev, bcoef)
    class(amrex_abeclaplacian), intent(inout) :: this
    integer, intent(in) :: lev
    type(amrex_multifab), intent(in) :: bcoef(amrex_spacedim)
    type(c_ptr) :: cp(amrex_spacedim)
    integer :: idim
    do idim = 1, amrex_spacedim
       cp(idim) = bcoef(idim)%p
    end do
    call amrex_fi_abeclap_set_bcoeffs(this%p, lev, cp)
  end subroutine amrex_abeclaplacian_set_bcoeffs

end module amrex_abeclaplacian_module
