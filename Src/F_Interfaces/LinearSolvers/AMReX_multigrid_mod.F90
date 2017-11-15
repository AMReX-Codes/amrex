
module amrex_multigrid_module

  use amrex_base_module
  use amrex_linop_module
  implicit none

  private
  public :: amrex_multigrid_build, amrex_multigrid_destroy

  type, public :: amrex_multigrid
     logical     :: owner = .false.
     type(c_ptr) :: p = c_null_ptr
   contains
     generic :: assignment(=) => amrex_multigrid_assign  ! shallow copy
     procedure, private :: amrex_multigrid_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_multigrid_assign
#endif
  end type amrex_multigrid

  ! interfaces to C++ functions

  interface
     subroutine amrex_fi_new_multigrid (multigrid, linop) bind(c)
       import
       implicit none
       type(c_ptr), intent(inout) :: multigrid
       type(c_ptr), value :: linop
     end subroutine amrex_fi_new_multigrid
  
     subroutine amrex_fi_delete_multigrid (multigrid) bind(c)
       import
       implicit none
       type(c_ptr), value :: multigrid
     end subroutine amrex_fi_delete_multigrid
  end interface

contains

  subroutine amrex_multigrid_assign (dst, src)
    class(amrex_multigrid), intent(inout) :: dst
    type(amrex_multigrid), intent(in) :: src
    call amrex_multigrid_destroy(dst)
    ! .....
  end subroutine amrex_multigrid_assign


  subroutine amrex_multigrid_build (multigrid, linop)
    type(amrex_multigrid), intent(inout) :: multigrid
    class(amrex_linop), intent(in) :: linop
    multigrid%owner = .true.
    call amrex_fi_new_multigrid(multigrid%p, linop%p)
  end subroutine amrex_multigrid_build


  subroutine amrex_multigrid_destroy (this)
    type(amrex_multigrid), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_multigrid(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
  end subroutine amrex_multigrid_destroy

end module amrex_multigrid_module
