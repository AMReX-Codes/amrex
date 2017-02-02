module amrex_famrcore_module

  use iso_c_binding
  use amrex_module
 
  implicit none

  private

!  public :: amrex_famrcore_build, amrex_famrcore_destroy

  type, abstract, public :: amrex_famrcore
     type(c_ptr) :: p = c_null_ptr
   contains
     procedure  ::  famrcore_build
     procedure  ::  famrcore_destroy
  end type amrex_famrcore

  interface
     subroutine amrex_fi_new_famrcore (famrcore) bind(c)
       import
       implicit none
       type(c_ptr) :: famrcore
     end subroutine amrex_fi_new_famrcore

     subroutine amrex_fi_delete_famrcore (famrcore) bind(c)
       import
       implicit none
       type(c_ptr), value :: famrcore
     end subroutine amrex_fi_delete_famrcore
  end interface

contains

  subroutine famrcore_build (this)
    class(amrex_famrcore) :: this
    print *, 'in famrcore_build'
    call amrex_fi_new_famrcore(this%p)
  end subroutine famrcore_build

  subroutine famrcore_destroy (this)
    class(amrex_famrcore) :: this
    print *, 'in famrcore_destroy'
    call amrex_fi_delete_famrcore(this%p)
    this%p = c_null_ptr
  end subroutine famrcore_destroy

end module amrex_famrcore_module
