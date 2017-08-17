
module amrex_tagbox_module

  use iso_c_binding
  use amrex_base_module

  implicit none

  private

  type, public :: amrex_tagboxarray
     type(c_ptr) :: p = c_null_ptr
   contains
     generic   :: assignment(=) => amrex_tagboxarray_assign, amrex_tagboxarray_install
     procedure :: dataPtr       => amrex_tagboxarray_dataptr
     procedure, private :: amrex_tagboxarray_assign
     procedure, private :: amrex_tagboxarray_install
  end type amrex_tagboxarray

  interface
     subroutine amrex_fi_tagboxarray_dataptr (tag, mfi, dp, lo, hi) bind(c)
       import
       implicit none
       type(c_ptr), value :: tag, mfi
       type(c_ptr) :: dp
       integer(c_int) :: lo(3), hi(3)
     end subroutine amrex_fi_tagboxarray_dataptr
  end interface

contains

  subroutine amrex_tagboxarray_assign (dst, src)
    class(amrex_tagboxarray), intent(inout) :: dst
    type(amrex_tagboxarray), intent(in) :: src
    dst%p = src%p
  end subroutine amrex_tagboxarray_assign

  subroutine amrex_tagboxarray_install (this, p)
    class(amrex_tagboxarray), intent(inout) :: this
    type(c_ptr), intent(in), value :: p
    this%p = p
  end subroutine amrex_tagboxarray_install

  function amrex_tagboxarray_dataPtr (this, mfi) result(dp)
    class(amrex_tagboxarray) :: this
    type(amrex_mfiter), intent(in) :: mfi
    character(kind=c_char), contiguous, pointer, dimension(:,:,:,:) :: dp
    type(c_ptr) :: cp
    character(kind=c_char), contiguous, pointer :: fp(:,:,:,:)
    integer(c_int) :: n(4)
    type(amrex_box) :: bx
    call amrex_fi_tagboxarray_dataptr(this%p, mfi%p, cp, bx%lo, bx%hi)
    n(1:3) = bx%hi - bx%lo + 1
    n(4)   = 1
    call c_f_pointer(cp, fp, shape=n)
    dp(bx%lo(1):,bx%lo(2):,bx%lo(3):,1:) => fp
  end function amrex_tagboxarray_dataPtr
  
end module amrex_tagbox_module
