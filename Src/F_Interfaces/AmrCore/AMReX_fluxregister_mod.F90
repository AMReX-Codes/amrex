
module amrex_fluxregister_module

  use iso_c_binding
  use amrex_base_module
  
  implicit none

  private

  public :: amrex_fluxregister_build, amrex_fluxregister_destroy

  type, public :: amrex_fluxregister
     logical     :: owner = .false.
     type(c_ptr) :: p = c_null_ptr
   contains
     generic   :: assignment(=) => amrex_fluxregister_assign   ! shallow copy
     procedure, private :: amrex_fluxregister_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_fluxregister_destroy
#endif
  end type amrex_fluxregister

  interface
     subroutine amrex_fi_new_fluxregister (fr, ba, dm, rr, flev, nc) bind(c)
       import
       implicit none
       type(c_ptr) :: fr
       type(c_ptr), value :: ba, dm
       integer, value :: rr, flev, nc
     end subroutine amrex_fi_new_fluxregister

     subroutine amrex_fi_delete_fluxregister (fr) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
     end subroutine amrex_fi_delete_fluxregister
  end interface

contains

  subroutine amrex_fluxregister_build (fr, ba, dm, ref_ratio, fine_lev, ncomp)
    type(amrex_fluxregister) :: fr
    type(amrex_boxarray), intent(in) :: ba
    type(amrex_distromap), intent(in) :: dm
    integer, intent(in) :: ref_ratio, fine_lev, ncomp
    fr%owner = .true.
    call amrex_fi_new_fluxregister(fr%p, ba%p, dm%p, ref_ratio, fine_lev, ncomp)
  end subroutine amrex_fluxregister_build
  
  impure elemental subroutine amrex_fluxregister_destroy (this)
    type(amrex_fluxregister), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_fluxregister(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
  end subroutine amrex_fluxregister_destroy
    
  subroutine amrex_fluxregister_assign (dst, src)
    class(amrex_fluxregister), intent(inout) :: dst
    type (amrex_fluxregister), intent(in   ) :: src
    dst%owner = .false.
    dst%p = src%p
  end subroutine amrex_fluxregister_assign


end module amrex_fluxregister_module

