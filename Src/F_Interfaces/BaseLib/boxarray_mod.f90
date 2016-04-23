
module boxarray_module

  use iso_c_binding

  use box_module

  implicit none

  private

  public :: boxarray_build

  type, public :: BoxArray
     logical     :: owner = .false.
     type(c_ptr) :: p = c_null_ptr
   contains
     procedure :: move => boxarray_move 
     procedure :: boxarray_assign
     generic   :: assignment(=) => boxarray_assign
     final :: boxarray_destroy
  end type BoxArray

  interface boxarray_build
     module procedure boxarray_build_bx
  end interface boxarray_build

  ! interfaces to cpp functions

  interface
     subroutine fi_new_boxarray (p,lo,hi) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), intent(out) :: p
       integer(c_int), intent(in) :: lo(3), hi(3)
     end subroutine fi_new_boxarray

     subroutine fi_delete_boxarray (p) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: p
     end subroutine fi_delete_boxarray
  end interface

contains

  subroutine boxarray_build_bx (ba, bx)
    type(BoxArray), intent(out) :: ba
    type(Box)     , intent(in ) :: bx
    ba%owner = .true.
    call fi_new_boxarray(ba%p, bx%lo, bx%hi)
  end subroutine boxarray_build_bx

  subroutine boxarray_destroy (this)
    type(BoxArray) :: this
    if (this%owner) then
       this%owner = .false.
       if (c_associated(this%p)) then
          call fi_delete_boxarray(this%p)
          this%p = c_null_ptr
       end if
    end if
  end subroutine boxarray_destroy

  subroutine boxarray_assign (dst, src)
    class(BoxArray), intent(out) :: dst
    type (BoxArray), intent(in ) :: src
    dst%owner = .false.
    dst%p = src%p
  end subroutine boxarray_assign

  subroutine boxarray_move (dst, src)
    class(BoxArray), intent(out)   :: dst
    type (BoxArray), intent(inout) :: src
    dst%owner = src%owner
    dst%p = src%p
    src%owner = .false.
    src%p = c_null_ptr
  end subroutine boxarray_move


end module boxarray_module
