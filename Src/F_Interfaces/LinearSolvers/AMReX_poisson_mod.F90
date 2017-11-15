
module amrex_poisson_module

  use amrex_base_module
  use amrex_linop_module
  implicit none

  private
  public :: amrex_poisson, amrex_poisson_build, amrex_poisson_destroy

  type, extends(amrex_linop), public :: amrex_poisson
   contains
     generic :: assignment(=) => amrex_poisson_assign   ! shallow copy
     procedure, private :: amrex_poisson_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_poisson_destroy
#endif
  end type amrex_poisson

contains

  subroutine amrex_poisson_assign (dst, src)
    class(amrex_poisson), intent(inout) :: dst
    type (amrex_poisson), intent(in   ) :: src
    call amrex_poisson_destroy(dst)
    !
  end subroutine amrex_poisson_assign


  subroutine amrex_poisson_build (linop)
    type(amrex_poisson), intent(inout) :: linop
  end subroutine amrex_poisson_build


  subroutine amrex_poisson_destroy (this)
    type(amrex_poisson), intent(inout) :: this
    !
  end subroutine amrex_poisson_destroy

end module amrex_poisson_module
