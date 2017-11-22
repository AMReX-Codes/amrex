module amrex_linop_module

  use amrex_base_module
  implicit none

  private

  type, public :: amrex_linop
     logical     :: owner =.false.
     type(c_ptr) :: p = c_null_ptr
   contains
     procedure :: set_maxorder             => amrex_linop_set_maxorder
     procedure :: set_domain_bc            => amrex_linop_set_domain_bc
     procedure :: set_coarse_fine_bc       => amrex_linop_set_coarse_fine_bc
     procedure :: set_level_bc             => amrex_linop_set_level_bc
     procedure, private :: amrex_linop_set_maxorder
     procedure, private :: amrex_linop_set_domain_bc
     procedure, private :: amrex_linop_set_coarse_fine_bc
     procedure, private :: amrex_linop_set_level_bc
  end type amrex_linop

  ! interfaces to C++ functions

  interface
     subroutine amrex_fi_linop_set_maxorder (linop, ord) bind(c)
       import
       implicit none
       type(c_ptr), value :: linop
       integer(c_int), value :: ord
     end subroutine amrex_fi_linop_set_maxorder

     subroutine amrex_fi_linop_set_domain_bc (linop, lobc, hibc) bind(c)
       import
       implicit none
       type(c_ptr), value :: linop
       integer(c_int), intent(in) :: lobc(*), hibc(*)
     end subroutine amrex_fi_linop_set_domain_bc

     subroutine amrex_fi_linop_set_coarse_fine_bc (linop, mf, crse_ratio) bind(c)
       import
       implicit none
       type(c_ptr), value :: linop, mf
       integer(c_int), value :: crse_ratio
     end subroutine amrex_fi_linop_set_coarse_fine_bc

     subroutine amrex_fi_linop_set_level_bc (linop, lev, mf) bind(c)
       import
       implicit none
       type(c_ptr), value :: linop, mf
       integer(c_int), value :: lev
     end subroutine amrex_fi_linop_set_level_bc
  end interface

contains

  subroutine amrex_linop_set_maxorder (this, ord)
    class(amrex_linop), intent(inout) :: this
    integer(c_int), intent(in) :: ord
    call amrex_fi_linop_set_maxorder(this%p, ord)
  end subroutine amrex_linop_set_maxorder


  subroutine amrex_linop_set_domain_bc (this, lobc, hibc)
    class(amrex_linop), intent(inout) :: this
    integer, dimension(amrex_spacedim), intent(in) :: lobc, hibc
    call amrex_fi_linop_set_domain_bc(this%p, lobc, hibc)
  end subroutine amrex_linop_set_domain_bc


  subroutine amrex_linop_set_coarse_fine_bc (this, crse_bcdata, crse_ratio)
    class(amrex_linop), intent(inout) :: this
    type(amrex_multifab), intent(in) :: crse_bcdata
    integer, intent(in) :: crse_ratio
    call amrex_fi_linop_set_coarse_fine_bc(this%p, crse_bcdata%p, crse_ratio)
  end subroutine amrex_linop_set_coarse_fine_bc


  subroutine amrex_linop_set_level_bc (this, lev, bcdata)
    class(amrex_linop), intent(inout) :: this
    integer, intent(in) :: lev
    type(amrex_multifab), intent(in) :: bcdata
    call amrex_fi_linop_set_level_bc(this%p, lev, bcdata%p)
  end subroutine amrex_linop_set_level_bc
  
end module amrex_linop_module
