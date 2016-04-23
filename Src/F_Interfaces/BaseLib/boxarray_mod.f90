
module boxarray_module

  use iso_c_binding

  use box_module

  implicit none

  private

  type, public :: BoxArray
     type(c_ptr) :: p = c_null_ptr
   contains
     final :: destroy_boxarray
  end type BoxArray

  interface BoxArray
     module procedure build_boxarray
  end interface BoxArray

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

  function build_boxarray (bx) result(ba)
    type(box), target, intent(in) :: bx
    type(BoxArray) :: ba
    call fi_new_boxarray(ba%p, bx%lo, bx%hi)
  end function build_boxarray

  subroutine destroy_boxarray (this)
    type(BoxArray) :: this
    if (c_associated(this%p)) then
       call fi_delete_boxarray(this%p)
       this%p = c_null_ptr
    end if
  end subroutine destroy_boxarray

end module boxarray_module
