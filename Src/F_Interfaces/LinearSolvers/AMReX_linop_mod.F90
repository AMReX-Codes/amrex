module amrex_linop_module

  use amrex_base_module
  implicit none

  private
  public :: amrex_linop, amrex_linop_build, amrex_linop_destroy

  type, public :: amrex_linop
     logical     :: owner =.false.
     type(c_ptr) :: p = c_null_ptr
   contains
     generic :: assignment(=) => amrex_linop_assign   ! shallow copy
     procedure, private :: amrex_linop_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_linop_destroy
#endif
  end type amrex_linop

contains

  subroutine amrex_linop_assign (dst, src)
    class(amrex_linop), intent(inout) :: dst
    type (amrex_linop), intent(in   ) :: src
    call amrex_linop_destroy(dst)
    dst%owner = .false.
    dst%p     = src%p
  end subroutine amrex_linop_assign


  subroutine amrex_linop_build (linop)
    type(amrex_linop), intent(inout) :: linop
  end subroutine amrex_linop_build


  subroutine amrex_linop_destroy (this)
    type(amrex_linop), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          ! call amrex_fi_delete_linop
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
  end subroutine amrex_linop_destroy

end module amrex_linop_module
