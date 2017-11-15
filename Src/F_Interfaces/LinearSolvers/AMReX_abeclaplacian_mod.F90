
module amrex_abeclapacian_module

  use amrex_base_module
  use amrex_linop_module
  implicit none

  private
  public :: amrex_abeclapacian, amrex_abeclapacian_build, amrex_abeclapacian_destroy

  type, extends(amrex_linop), public :: amrex_abeclapacian
   contains
     generic :: assignment(=) => amrex_abeclapacian_assign   ! shallow copy
     procedure, private :: amrex_abeclapacian_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_abeclapacian_destroy
#endif
  end type amrex_abeclapacian

contains

  subroutine amrex_abeclapacian_assign (dst, src)
    class(amrex_abeclapacian), intent(inout) :: dst
    type (amrex_abeclapacian), intent(in   ) :: src
    call amrex_abeclapacian_destroy(dst)
    !
  end subroutine amrex_abeclapacian_assign


  subroutine amrex_abeclapacian_build (linop)
    type(amrex_abeclapacian), intent(inout) :: linop
  end subroutine amrex_abeclapacian_build


  subroutine amrex_abeclapacian_destroy (this)
    type(amrex_abeclapacian), intent(inout) :: this
    !
  end subroutine amrex_abeclapacian_destroy

end module amrex_abeclapacian_module
