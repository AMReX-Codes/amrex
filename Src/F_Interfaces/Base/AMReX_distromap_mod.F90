
module amrex_distromap_module

  use amrex_boxarray_module
  use iso_c_binding

  implicit none

  private

  public :: amrex_distromap_build, amrex_distromap_destroy, amrex_print

  type, public :: amrex_distromap
     logical     :: owner = .false.
     type(c_ptr) :: p = c_null_ptr
   contains
     generic   :: assignment(=) => amrex_distromap_assign, amrex_distromap_install   ! shallow copy
     procedure :: clone         => amrex_distromap_clone    ! deep copy
     procedure :: move          => amrex_distromap_move     ! transfer ownership
     procedure :: get_pmap      => amrex_distromap_get_pmap ! fill caller-owned array of PEs
     procedure, private :: amrex_distromap_assign
     procedure, private :: amrex_distromap_install 
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_distromap_destroy
#endif
  end type amrex_distromap

  interface amrex_distromap_build
     module procedure amrex_distromap_build_ba
     module procedure amrex_distromap_build_pmap
  end interface amrex_distromap_build

  interface amrex_print
     module procedure amrex_distromap_print
  end interface amrex_print

  ! interfaces to cpp functions

  interface
     subroutine amrex_fi_new_distromap (dm,ba) bind(c)
       import
       implicit none
       type(c_ptr) :: dm
       type(c_ptr), value :: ba
     end subroutine amrex_fi_new_distromap

     subroutine amrex_fi_new_distromap_from_pmap (dm,pmap,plen) bind(c)
       import
       implicit none
       type(c_ptr) :: dm
       integer(c_int), intent(in) :: pmap(*)
       integer(c_int), value :: plen
     end subroutine amrex_fi_new_distromap_from_pmap

     subroutine amrex_fi_delete_distromap (dm) bind(c)
       import
       implicit none
       type(c_ptr), value :: dm
     end subroutine amrex_fi_delete_distromap

     subroutine amrex_fi_clone_distromap (dmo, dmi) bind(c)
       import
       implicit none
       type(c_ptr) :: dmo
       type(c_ptr), value :: dmi
     end subroutine amrex_fi_clone_distromap

     subroutine amrex_fi_distromap_maxsize (dm,n) bind(c)
       import
       implicit none
       type(c_ptr), value :: dm
       integer(c_int), value :: n
     end subroutine amrex_fi_distromap_maxsize

     subroutine amrex_fi_distromap_get_pmap (dm,pmap,plen) bind(c)
       import
       implicit none
       type(c_ptr), value :: dm
       integer(c_int), intent(out) :: pmap(*)
       integer(c_int), value :: plen
     end subroutine amrex_fi_distromap_get_pmap

     subroutine amrex_fi_print_distromap (dm) bind(c)
       import
       implicit none
       type(c_ptr), value :: dm
     end subroutine amrex_fi_print_distromap
  end interface

contains

  subroutine amrex_distromap_build_ba (dm, ba)
    type(amrex_distromap) :: dm
    type(amrex_boxarray), intent(in) :: ba
    dm%owner = .true.
    call amrex_fi_new_distromap(dm%p, ba%p)
  end subroutine amrex_distromap_build_ba

  subroutine amrex_distromap_build_pmap (dm, pmap)
    type(amrex_distromap) :: dm
    integer, intent(in) :: pmap(:)
    dm%owner = .true.
    call amrex_fi_new_distromap_from_pmap(dm%p,pmap, size(pmap))
  end subroutine amrex_distromap_build_pmap

  impure elemental subroutine amrex_distromap_destroy (this)
    type(amrex_distromap), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_distromap(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
  end subroutine amrex_distromap_destroy

  subroutine amrex_distromap_assign (dst, src)
    class(amrex_distromap), intent(inout) :: dst
    type (amrex_distromap), intent(in   ) :: src
    dst%owner = .false.
    dst%p = src%p
  end subroutine amrex_distromap_assign

  subroutine amrex_distromap_install (this, p)
    class(amrex_distromap), intent(inout) :: this
    type(c_ptr), intent(in) :: p
    this%owner = .false.
    this%p     = p
  end subroutine amrex_distromap_install

  subroutine amrex_distromap_clone (dst, src)
    class(amrex_distromap), intent(inout) :: dst
    type (amrex_distromap), intent(in   ) :: src
    dst%owner = .true.
    call amrex_fi_clone_distromap(dst%p, src%p)
  end subroutine amrex_distromap_clone

  subroutine amrex_distromap_move (dst, src)
    class(amrex_distromap) :: dst, src
    dst%owner = src%owner
    dst%p = src%p
    src%owner = .false.
    src%p = c_null_ptr
  end subroutine amrex_distromap_move

  subroutine amrex_distromap_get_pmap (dm, pmap)
    class(amrex_distromap) :: dm
    integer, intent(out) :: pmap(:)
    call amrex_fi_distromap_get_pmap(dm%p,pmap, size(pmap))
  end subroutine amrex_distromap_get_pmap

  subroutine amrex_distromap_print (dm)
    type(amrex_distromap), intent(in) :: dm
    call amrex_fi_print_distromap(dm%p)
  end subroutine amrex_distromap_print

end module amrex_distromap_module
