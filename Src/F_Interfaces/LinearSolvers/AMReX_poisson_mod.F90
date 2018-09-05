
module amrex_poisson_module

  use amrex_base_module
  use amrex_linop_module
  implicit none

  private
  public :: amrex_poisson_build, amrex_poisson_destroy

  type, extends(amrex_linop), public :: amrex_poisson
   contains
     generic :: assignment(=) => amrex_poisson_assign   ! shallow copy
     procedure, private :: amrex_poisson_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_poisson_destroy
#endif
  end type amrex_poisson

  ! interfaces to C++ functions
  
  interface
     subroutine amrex_fi_new_poisson (poisson,n,gm,ba,dm,mt,agg,con,mcl) bind(c)
       import
       implicit none
       type(c_ptr), intent(inout) :: poisson
       integer(c_int), value :: n, mt, agg, con, mcl
       type(c_ptr), intent(in) :: gm(*), ba(*), dm(*)
     end subroutine amrex_fi_new_poisson

     subroutine amrex_fi_delete_linop (linop) bind(c)
       import
       implicit none
       type(c_ptr), value :: linop
     end subroutine amrex_fi_delete_linop
  end interface

contains

  subroutine amrex_poisson_assign (dst, src)
    class(amrex_poisson), intent(inout) :: dst
    type (amrex_poisson), intent(in   ) :: src
    call amrex_poisson_destroy(dst)
    dst%owner = .false.
    dst%p = src%p
  end subroutine amrex_poisson_assign


  subroutine amrex_poisson_build (poisson, geom, ba, dm, metric_term, &
       agglomeration, consolidation, max_coarsening_level)
    type(amrex_poisson), intent(inout) :: poisson
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

    poisson%owner = .true.
    call amrex_fi_new_poisson(poisson%p, nlevs, gm_cp, ba_cp, dm_cp, imt, iagg, icon, imcl)
  end subroutine amrex_poisson_build


  subroutine amrex_poisson_destroy (this)
    type(amrex_poisson), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_linop(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
  end subroutine amrex_poisson_destroy

end module amrex_poisson_module
