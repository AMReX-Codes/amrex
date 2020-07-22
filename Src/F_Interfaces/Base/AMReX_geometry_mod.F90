
module amrex_geometry_module

  use iso_c_binding
  use amrex_fort_module, only : ndims => amrex_spacedim, amrex_real
  use amrex_box_module
  use amrex_parmparse_module

  implicit none

  private

  public :: amrex_geometry_set_coord_sys, amrex_geometry_set_prob_domain, amrex_geometry_set_periodic
  public :: amrex_pmask, amrex_problo, amrex_probhi
  public :: amrex_geometry_build, amrex_geometry_destroy, amrex_geometry_init_data
  public :: amrex_is_periodic, amrex_is_any_periodic, amrex_is_all_periodic
  public :: amrex_geometry_finalize

  logical, save :: amrex_pmask(3)  = .false.  
  real(amrex_real), save :: amrex_problo(3) = 0.0_amrex_real
  real(amrex_real), save :: amrex_probhi(3) = 1.0_amrex_real

  logical, save :: amrex_geometry_initialized = .false.

!$omp threadprivate(amrex_geometry_initialized,amrex_pmask,amrex_problo,amrex_probhi)

  type, public :: amrex_geometry
     logical          :: owner     = .false.
     type(c_ptr)      :: p         = c_null_ptr
     !
     real(amrex_real) :: dx(3)     = 0.0_amrex_real
     type(amrex_box)  :: domain
   contains
     generic :: assignment(=)           => amrex_geometry_assign, amrex_geometry_install  ! shallow copy
     procedure :: get_physical_location => amrex_geometry_get_ploc
     procedure, private :: amrex_geometry_assign
     procedure, private :: amrex_geometry_install
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_geometry_destroy
#endif
  end type amrex_geometry

  ! interfaces to c++ functions

  interface
     subroutine amrex_fi_new_geometry (geom,lo,hi) bind(c)
       import
       implicit none
       type(c_ptr) :: geom
       integer, intent(in) :: lo(3), hi(3)
     end subroutine amrex_fi_new_geometry

     subroutine amrex_fi_delete_geometry (geom) bind(c)
       import
       implicit none
       type(c_ptr), value :: geom
     end subroutine amrex_fi_delete_geometry

     subroutine amrex_fi_geometry_get_pmask (pmask) bind(c)
       import
       implicit none
       integer(c_int) :: pmask(3)
     end subroutine amrex_fi_geometry_get_pmask

     subroutine amrex_fi_geometry_get_probdomain (problo,probhi) bind(c)
       import
       implicit none
       real(amrex_real) :: problo(3), probhi(3)
     end subroutine amrex_fi_geometry_get_probdomain

     subroutine amrex_fi_geometry_get_intdomain (geom,lo,hi) bind(c)
       import
       implicit none
       type(c_ptr), value :: geom
       integer(c_int), intent(out) :: lo(3), hi(3)
     end subroutine amrex_fi_geometry_get_intdomain
  end interface

contains

  subroutine amrex_geometry_finalize ()
    amrex_pmask(3)  = .false.
    amrex_problo(3) = 0.0_amrex_real
    amrex_probhi(3) = 1.0_amrex_real
    amrex_geometry_initialized = .false.
  end subroutine amrex_geometry_finalize

  subroutine amrex_geometry_set_coord_sys (csys)
    integer, intent(in) :: csys
    type(amrex_parmparse) :: pp
    call amrex_parmparse_build(pp,"geometry")
    call pp%add("coord_sys", csys)
    call amrex_parmparse_destroy(pp)
  end subroutine amrex_geometry_set_coord_sys

  subroutine amrex_geometry_set_prob_domain (pblo, pbhi)
    real(amrex_real), intent(in) :: pblo(3), pbhi(3)
    type(amrex_parmparse) :: pp
    call amrex_parmparse_build(pp,"geometry")
    call pp%addarr("prob_lo", pblo)
    call pp%addarr("prob_hi", pbhi)
    call amrex_parmparse_destroy(pp)
  end subroutine amrex_geometry_set_prob_domain

  subroutine amrex_geometry_set_periodic (periodic)
    logical, intent(in) :: periodic(3)
    integer :: imask(3)
    type(amrex_parmparse) :: pp
    call amrex_parmparse_build(pp,"geometry")
    imask = 0
    where(periodic) imask = 1
    call pp%addarr("is_periodic", imask)
    call amrex_parmparse_destroy(pp)
  end subroutine amrex_geometry_set_periodic

  subroutine amrex_geometry_init ()
    integer :: imask(3)
    imask = 0
    call amrex_fi_geometry_get_pmask(imask)
    where (imask .eq. 1) amrex_pmask = .true.
    call amrex_fi_geometry_get_probdomain(amrex_problo, amrex_probhi)
  end subroutine amrex_geometry_init

  subroutine amrex_geometry_build (geom, domain)
    type(amrex_geometry) :: geom
    type(amrex_box), intent(in) :: domain
    geom%owner = .true.
    call amrex_fi_new_geometry(geom%p, domain%lo, domain%hi)
    call amrex_geometry_init_data(geom)
  end subroutine amrex_geometry_build

  subroutine amrex_geometry_init_data (geom)  ! geom%p must be valid!
    type(amrex_geometry), intent(inout) :: geom
    integer :: i, lo(3), hi(3)
    logical, external :: omp_in_parallel
    !$omp parallel if(.not.omp_in_parallel())
    if (.not.amrex_geometry_initialized) then
       call amrex_geometry_init()
       amrex_geometry_initialized = .true.
    end if
    !$omp end parallel
    call amrex_fi_geometry_get_intdomain(geom%p, lo, hi)
    geom%domain = amrex_box(lo, hi)
    do i = 1, ndims
       geom%dx(i) = (amrex_probhi(i)-amrex_problo(i)) / dble(hi(i)-lo(i)+1)
    end do
  end subroutine amrex_geometry_init_data

  impure elemental subroutine amrex_geometry_destroy (geom)
    type(amrex_geometry), intent(inout) :: geom
    if (geom%owner) then
       if (c_associated(geom%p)) then
          call amrex_fi_delete_geometry(geom%p)
       end if
    end if
    geom%owner = .false.
    geom%p = c_null_ptr
  end subroutine amrex_geometry_destroy

  subroutine amrex_geometry_assign (dst, src)
    class(amrex_geometry), intent(inout) :: dst
    type (amrex_geometry), intent(in   ) :: src
    call amrex_geometry_destroy(dst)
    dst%owner  = .false.
    dst%p      = src%p
    dst%dx     = src%dx
    dst%domain = src%domain
  end subroutine amrex_geometry_assign

  subroutine amrex_geometry_install (this, p)
    class(amrex_geometry), intent(inout) :: this
    type(c_ptr), intent(in) :: p
    this%owner  = .false.
    this%p      = p
    call amrex_geometry_init_data(this)
  end subroutine amrex_geometry_install

  logical function amrex_is_periodic (dir)
    integer, intent(in) :: dir
    amrex_is_periodic = amrex_pmask(dir)
  end function amrex_is_periodic

  logical function amrex_is_any_periodic ()
    amrex_is_any_periodic = any(amrex_pmask(1:ndims))
  end function amrex_is_any_periodic

  logical function amrex_is_all_periodic ()
    amrex_is_all_periodic = all(amrex_pmask(1:ndims))
  end function amrex_is_all_periodic

  ! give integer location, compute real location
  function amrex_geometry_get_ploc (this, iloc) result (ploc)
    class(amrex_geometry), intent(in) :: this
    integer, intent(in) :: iloc(ndims)
    real(amrex_real) :: ploc(ndims)
    integer :: i
    do i = 1, ndims
       ploc(i) = amrex_problo(i) + (iloc(i) - this%domain%lo(i)) * this%dx(i)
    end do
  end function amrex_geometry_get_ploc
    
end module amrex_geometry_module
