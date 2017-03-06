module amrex_famrcore_module

  use iso_c_binding
  use amrex_module
 
  implicit none

  private

  ! public routines
  public :: amrex_init_from_scratch, amrex_init_virtual_functions
  public :: amrex_famrcore_init, amrex_famrcore_finalize, amrex_famrcore_initialized, &
       amrex_get_finest_level, amrex_get_boxarray, amrex_get_distromap, amrex_get_geometry, &
       amrex_set_finest_level, amrex_set_boxarray, amrex_set_distromap, &
       amrex_make_base_grids, amrex_make_new_grids, &
       amrex_install_level, amrex_remove_level

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

     integer(c_int) function amrex_fi_get_finest_level (famrcore) bind(c)
       import
       type(c_ptr), value :: famrcore
     end function amrex_fi_get_finest_level

     subroutine amrex_fi_get_boxarray (ba, lev, famrcore) bind(c)
       import
       type(c_ptr), intent(out) :: ba
       integer(c_int), value :: lev
       type(c_ptr), value :: famrcore
     end subroutine amrex_fi_get_boxarray

     subroutine amrex_fi_get_distromap (dm, lev, famrcore) bind(c)
       import
       type(c_ptr), intent(out) :: dm
       integer(c_int), value :: lev
       type(c_ptr), value :: famrcore
     end subroutine amrex_fi_get_distromap

     subroutine amrex_fi_get_geometry (geom, lev, famrcore) bind(c)
       import
       type(c_ptr), intent(out) :: geom
       integer(c_int), value :: lev
       type(c_ptr), value :: famrcore
     end subroutine amrex_fi_get_geometry

     subroutine amrex_fi_set_finest_level (lev, famrcore) bind(c)
       import
       integer(c_int), value :: lev
       type(c_ptr), value :: famrcore
     end subroutine amrex_fi_set_finest_level

     subroutine amrex_fi_set_boxarray (lev, ba, famrcore) bind(c)
       import
       integer(c_int), value :: lev
       type(c_ptr), value :: ba
       type(c_ptr), value :: famrcore
     end subroutine amrex_fi_set_boxarray

     subroutine amrex_fi_set_distromap (lev, dm, famrcore) bind(c)
       import
       integer(c_int), value :: lev
       type(c_ptr), value :: dm
       type(c_ptr), value :: famrcore
     end subroutine amrex_fi_set_distromap
     
     subroutine amrex_fi_make_base_grids (ba, famrcore) bind(c)
       import
       type(c_ptr), intent(out) :: ba
       type(c_ptr), value :: famrcore
     end subroutine amrex_fi_make_base_grids

!     subroutine amrex_fi_make_new_grids (baselev, time, new_finest, ba_array)
!       
!     end subroutine amrex_fi_make_new_grids

     subroutine amrex_fi_init_from_scratch (t, famrcore) bind(c)
       import
       real(amrex_real), value :: t
       type(c_ptr), value :: famrcore
     end subroutine amrex_fi_init_from_scratch
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

  integer function amrex_get_finest_level ()
    amrex_get_finest_level = amrex_fi_get_finest_level(famrcore)
  end function amrex_get_finest_level

  function amrex_get_boxarray (lev) result(ba)
    integer, intent(in) :: lev
    type(amrex_boxarray) :: ba
    call amrex_fi_get_boxarray(ba%p, lev, famrcore)
  end function amrex_get_boxarray

  function amrex_get_distromap (lev) result(dm)
    integer, intent(in) :: lev
    type(amrex_distromap) :: dm
    call amrex_fi_get_distromap(dm%p, lev, famrcore)
  end function amrex_get_distromap

  function amrex_get_geometry (lev) result(geom)
    integer, intent(in) :: lev
    type(amrex_geometry) :: geom
    call amrex_fi_get_geometry(geom%p, lev, famrcore)
    call amrex_geometry_init_data(geom)
  end function amrex_get_geometry

  subroutine amrex_set_finest_level (lev)
    integer, intent(in) :: lev
    call amrex_fi_set_finest_level(lev, famrcore)
  end subroutine amrex_set_finest_level

  subroutine amrex_set_boxarray (lev, ba)
    integer, intent(in) :: lev
    type(amrex_boxarray), intent(in) :: ba
    call amrex_fi_set_boxarray(lev, ba%p, famrcore)
  end subroutine amrex_set_boxarray

  subroutine amrex_set_distromap (lev, dm)
    integer, intent(in) :: lev
    type(amrex_distromap), intent(in) :: dm
    call amrex_fi_set_distromap(lev, dm%p, famrcore)
  end subroutine amrex_set_distromap

  subroutine amrex_make_base_grids (ba)
    type(amrex_boxarray), intent(inout) :: ba
    call amrex_boxarray_destroy(ba)
    ba%owner = .true.
    call amrex_fi_make_base_grids(ba%p, famrcore)
  end subroutine amrex_make_base_grids

  subroutine amrex_make_new_grids ()
  end subroutine amrex_make_new_grids

  subroutine amrex_install_level (lev, ba, dm)
    integer, intent(in) :: lev
    type(amrex_boxarray), intent(in) :: ba
    type(amrex_distromap), intent(in) :: dm
    if (lev < amrex_get_finest_level()) call amrex_set_finest_level(lev)
    call amrex_set_boxarray(lev, ba)
    call amrex_set_distromap(lev, dm)
  end subroutine amrex_install_level

  subroutine amrex_remove_level (lev)
    integer, intent(in) :: lev
    call amrex_fi_set_boxarray(lev, c_null_ptr, famrcore)
    call amrex_fi_set_distromap(lev, c_null_ptr, famrcore)
  end subroutine amrex_remove_level

  subroutine amrex_init_from_scratch (t)
    real(amrex_real), intent(in) :: t
    call amrex_fi_init_from_scratch(t, famrcore)
  end subroutine amrex_init_from_scratch

  subroutine amrex_init_virtual_functions (mk_lev_scrtch)
    type(c_funptr), intent(in) :: mk_lev_scrtch
    call amrex_fi_init_virtual_functions (mk_lev_scrtch)
  end subroutine amrex_init_virtual_functions

end module amrex_famrcore_module

