module amrex_famrcore_module

  use iso_c_binding
  use amrex_module
 
  implicit none

  private

  type(c_ptr) :: famrcore = c_null_ptr

  public :: amrex_famrcore_init, amrex_famrcore_finalize, amrex_famrcore_initialized

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

  subroutine amrex_famrcore_init ()
    call amrex_fi_new_famrcore(famrcore)
  end subroutine amrex_famrcore_init

  subroutine amrex_famrcore_finalize ()
    call amrex_fi_delete_famrcore(famrcore)
    famrcore = c_null_ptr
  end subroutine amrex_famrcore_finalize

  logical function amrex_famrcore_initialized ()
    amrex_famrcore_initialized = c_associated(famrcore)
  end function amrex_famrcore_initialized

end module amrex_famrcore_module
