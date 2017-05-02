module amrex_physbc_module

  use iso_c_binding
  use amrex_fort_module
  use amrex_multifab_module
  use amrex_geometry_module

  implicit none

  private

  public :: amrex_physbc_build, amrex_physbc_destroy, amrex_physbc_proc

  type, public :: amrex_physbc
     logical     :: owner = .false.
     type(c_ptr) :: p     = c_null_ptr
   contains
     generic :: assignment(=) => amrex_physbc_assign  ! shallow copy
     procedure, private :: amrex_physbc_assign
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_physbc_destroy
#endif
  end type amrex_physbc

  interface
     subroutine amrex_physbc_proc (mf, scomp, ncomp, time, geom) bind(c)
       import
       type(c_ptr), value :: mf, geom
       integer(c_int), value :: scomp, ncomp
       real(amrex_real), value :: time
     end subroutine amrex_physbc_proc
  end interface

  interface
     subroutine amrex_fi_new_physbc (pbc, fill, geom) bind(c)
       import
       implicit none
       type(c_ptr) :: pbc
       type(c_ptr), value :: geom
       type(c_funptr), value :: fill
     end subroutine amrex_fi_new_physbc

     subroutine amrex_fi_delete_physbc (pbc) bind(c)
       import
       implicit none
       type(c_ptr), value :: pbc
     end subroutine amrex_fi_delete_physbc
  end interface

contains

  subroutine amrex_physbc_build (pbc, fill, geom)
    type(amrex_physbc), intent(inout) :: pbc
    procedure(amrex_physbc_proc) :: fill
    type(amrex_geometry), intent(in) :: geom
    pbc%owner = .true.
    call amrex_fi_new_physbc(pbc%p, c_funloc(fill), geom%p)
  end subroutine amrex_physbc_build

  impure elemental subroutine amrex_physbc_destroy (this)
    type(amrex_physbc), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_physbc(this%p)
       end if
    end if
    this%owner = .false.
    this%p     = c_null_ptr
  end subroutine amrex_physbc_destroy
  
  subroutine amrex_physbc_assign (dst, src)
    class(amrex_physbc), intent(inout) :: dst
    type (amrex_physbc), intent(in   ) :: src
    call amrex_physbc_destroy(dst)
    dst%owner = .false.
    dst%p     = src%p
  end subroutine amrex_physbc_assign

end module amrex_physbc_module

