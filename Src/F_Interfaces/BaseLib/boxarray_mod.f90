
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
     procedure :: move          => boxarray_move 
     procedure ::                  boxarray_assign
     procedure :: maxSize       => boxarray_maxSize
     generic   :: assignment(=) => boxarray_assign
     final :: boxarray_destroy
  end type BoxArray

  interface boxarray_build
     module procedure boxarray_build_bx
  end interface boxarray_build

  ! interfaces to cpp functions

  interface
     subroutine fi_new_boxarray (ba,lo,hi) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ba
       integer(c_int), intent(in) :: lo(3), hi(3)
     end subroutine fi_new_boxarray

     subroutine fi_delete_boxarray (ba) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: ba
     end subroutine fi_delete_boxarray

     subroutine fi_boxarray_maxsize (ba,n) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: ba
       integer(c_int), value, intent(in) :: n
     end subroutine fi_boxarray_maxsize
  end interface

contains

  subroutine boxarray_build_bx (ba, bx)
    type(BoxArray) :: ba
    type(Box), intent(in ) :: bx
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
    class(BoxArray), intent(inout) :: dst
    type (BoxArray), intent(in   ) :: src
    dst%owner = .false.
    dst%p = src%p
  end subroutine boxarray_assign

  subroutine boxarray_move (dst, src)
    class(BoxArray) :: dst, src
    dst%owner = src%owner
    dst%p = src%p
    src%owner = .false.
    src%p = c_null_ptr
  end subroutine boxarray_move

  subroutine boxarray_maxsize (this, sz)
    class(BoxArray) this
    integer, intent(in)    :: sz
    call fi_boxarray_maxsize(this%p, sz)
  end subroutine boxarray_maxsize

end module boxarray_module
