module amrex_famrcore_module

  use iso_c_binding
  use amrex_module
 
  implicit none

  private

  ! public routines
  public :: amrex_famrcore_init, amrex_famrcore_finalize, amrex_famrcore_initialized

  ! public variables
  public :: amrex_max_level, amrex_ref_ratio

  type(c_ptr) :: famrcore = c_null_ptr

  integer :: amrex_max_level
  integer, allocatable :: amrex_ref_ratio(:)

  interface
     subroutine amrex_fi_new_famrcore (famrcore) bind(c)
       import
       type(c_ptr) :: famrcore
     end subroutine amrex_fi_new_famrcore

     subroutine amrex_fi_delete_famrcore (famrcore) bind(c)
       import
       type(c_ptr), value :: famrcore
     end subroutine amrex_fi_delete_famrcore

     integer(c_int) function amrex_fi_get_max_level (famrcore) bind(c)
       import
       type(c_ptr), value :: famrcore
     end function amrex_fi_get_max_level

     subroutine amrex_fi_get_ref_ratio (ref_ratio, famrcore) bind(c)
       import
       integer, intent(inout) :: ref_ratio(*)
       type(c_ptr), value :: famrcore
     end subroutine amrex_fi_get_ref_ratio
  end interface

contains

  subroutine amrex_famrcore_init ()
    call amrex_fi_new_famrcore(famrcore)
    amrex_max_level = amrex_fi_get_max_level(famrcore)
    allocate(amrex_ref_ratio(0:amrex_max_level-1))
    call amrex_fi_get_ref_ratio(amrex_ref_ratio, famrcore)
  end subroutine amrex_famrcore_init

  subroutine amrex_famrcore_finalize ()
    call amrex_fi_delete_famrcore(famrcore)
    famrcore = c_null_ptr
  end subroutine amrex_famrcore_finalize

  logical function amrex_famrcore_initialized ()
    amrex_famrcore_initialized = c_associated(famrcore)
  end function amrex_famrcore_initialized

end module amrex_famrcore_module
