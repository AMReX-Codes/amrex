
module amrex_boxarray_module

  use iso_c_binding

  use amrex_box_module

  implicit none

  private

  public :: amrex_boxarray_build, amrex_boxarray_destroy

  type, public :: amrex_boxarray
     logical     :: owner = .false.
     type(c_ptr) :: p = c_null_ptr
   contains
     procedure :: move          => amrex_boxarray_move     ! take ownship
     procedure :: clone         => amrex_boxarray_clone    ! deep copy
     procedure ::                  amrex_boxarray_assign   ! shallow copy
     generic   :: assignment(=) => amrex_boxarray_assign   ! shallow copy
     procedure :: maxSize       => amrex_boxarray_maxSize  ! make the boxes smaller
#if defined(__gfortran__) && (__GNUC__ <= 4)
     final :: amrex_boxarray_destroy
#endif
  end type amrex_boxarray

  interface amrex_boxarray_build
     module procedure amrex_boxarray_build_bx
  end interface amrex_boxarray_build

  ! interfaces to cpp functions

  interface
     subroutine amrex_fi_new_boxarray (ba,lo,hi) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ba
       integer(c_int), intent(in) :: lo(3), hi(3)
     end subroutine amrex_fi_new_boxarray

     subroutine amrex_fi_delete_boxarray (ba) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: ba
     end subroutine amrex_fi_delete_boxarray

     subroutine amrex_fi_clone_boxarray (bao, bai) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: bao
       type(c_ptr), value, intent(in) :: bai
     end subroutine amrex_fi_clone_boxarray

     subroutine amrex_fi_boxarray_maxsize (ba,n) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: ba
       integer(c_int), value, intent(in) :: n
     end subroutine amrex_fi_boxarray_maxsize
  end interface

contains

  subroutine amrex_boxarray_build_bx (ba, bx)
    type(amrex_boxarray) :: ba
    type(amrex_box), intent(in ) :: bx
    ba%owner = .true.
    call amrex_fi_new_boxarray(ba%p, bx%lo, bx%hi)
  end subroutine amrex_boxarray_build_bx

  subroutine amrex_boxarray_destroy (this)
    type(amrex_boxarray) :: this
    if (this%owner) then
       this%owner = .false.
       if (c_associated(this%p)) then
          call amrex_fi_delete_boxarray(this%p)
          this%p = c_null_ptr
       end if
    end if
  end subroutine amrex_boxarray_destroy

  subroutine amrex_boxarray_assign (dst, src)
    class(amrex_boxarray), intent(inout) :: dst
    type (amrex_boxarray), intent(in   ) :: src
    call amrex_boxarray_destroy(dst)
    dst%owner = .false.
    dst%p = src%p
  end subroutine amrex_boxarray_assign

  subroutine amrex_boxarray_clone (dst, src)
    class(amrex_boxarray), intent(inout) :: dst
    type (amrex_boxarray), intent(in   ) :: src
    dst%owner = .true.
    call amrex_fi_clone_boxarray(dst%p, src%p)
  end subroutine amrex_boxarray_clone

  subroutine amrex_boxarray_move (dst, src)
    class(amrex_boxarray) :: dst, src
    call amrex_boxarray_destroy(dst)
    dst%owner = src%owner
    dst%p = src%p
    src%owner = .false.
    src%p = c_null_ptr
  end subroutine amrex_boxarray_move

  subroutine amrex_boxarray_maxsize (this, sz)
    class(amrex_boxarray) this
    integer, intent(in)    :: sz
    call amrex_fi_boxarray_maxsize(this%p, sz)
  end subroutine amrex_boxarray_maxsize

end module amrex_boxarray_module
