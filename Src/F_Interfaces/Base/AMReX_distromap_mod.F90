
module amrex_distromap_module

  use amrex_boxarray_module
  use iso_c_binding

  implicit none

  private

  public :: amrex_distromap_build, amrex_distromap_destroy

  type, public :: amrex_distromap
     logical     :: owner = .false.
     type(c_ptr) :: p = c_null_ptr
   contains
     generic   :: assignment(=) => amrex_distromap_assign   ! shallow copy
     procedure :: clone         => amrex_distromap_clone    ! deep copy
     procedure :: move          => amrex_distromap_move     ! transfer ownership
     procedure, private :: amrex_distromap_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_distromap_destroy
#endif
  end type amrex_distromap

  interface amrex_distromap_build
     module procedure amrex_distromap_build_ba
  end interface amrex_distromap_build

  ! interfaces to cpp functions

  interface
     subroutine amrex_fi_new_distromap (dm,ba) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: dm
       type(c_ptr), value :: ba
     end subroutine amrex_fi_new_distromap

     subroutine amrex_fi_delete_distromap (dm) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: dm
     end subroutine amrex_fi_delete_distromap

     subroutine amrex_fi_clone_distromap (dmo, dmi) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: dmo
       type(c_ptr), value :: dmi
     end subroutine amrex_fi_clone_distromap
  end interface

contains

  subroutine amrex_distromap_build_ba (dm, ba)
    type(amrex_distromap) :: dm
    type(amrex_boxarray), intent(in) :: ba
    dm%owner = .true.
    call amrex_fi_new_distromap(dm%p, ba%p)
  end subroutine amrex_distromap_build_ba

  subroutine amrex_distromap_destroy (this)
    type(amrex_distromap) :: this
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

end module amrex_distromap_module
