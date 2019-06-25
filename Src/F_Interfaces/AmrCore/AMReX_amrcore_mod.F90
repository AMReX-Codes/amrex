module amrex_amrcore_module

  use iso_c_binding
  use amrex_base_module
 
  implicit none

  private

  ! public routines
  public :: amrex_amrcore_init, amrex_amrcore_finalize, amrex_amrcore_initialized, &
       amrex_get_amrcore, amrex_get_finest_level, amrex_get_numlevels, &
       amrex_get_boxarray, amrex_get_distromap, amrex_get_geometry, &
       amrex_set_boxarray, amrex_set_distromap, amrex_set_geometry, &
       amrex_set_finest_level, &
       amrex_init_from_scratch, amrex_init_virtual_functions, amrex_regrid, &
       amrex_init_post_regrid_function

  ! public variables
  public :: amrex_max_level, amrex_ref_ratio, amrex_geom

  ! public procedure interfaces
  public :: amrex_make_level_proc, amrex_clear_level_proc, amrex_error_est_proc

  interface
     subroutine amrex_make_level_proc (lev, time, ba, dm) bind(c)
       import
       implicit none
       integer, intent(in), value :: lev
       real(amrex_real), intent(in), value :: time
       type(c_ptr), intent(in), value :: ba, dm
     end subroutine amrex_make_level_proc

     subroutine amrex_clear_level_proc (lev) bind(c)
       import
       implicit none
       integer, intent(in) , value :: lev
     end subroutine amrex_clear_level_proc

     subroutine amrex_error_est_proc (lev, tags, time, tagval, clearval) bind(c)
       import
       implicit none
       integer, intent(in), value :: lev
       type(c_ptr), intent(in), value :: tags
       real(amrex_real), intent(in), value :: time
       character(kind=c_char), intent(in), value :: tagval, clearval
     end subroutine amrex_error_est_proc
  end interface

  interface
     subroutine amrex_post_regrid_proc ()
     end subroutine amrex_post_regrid_proc
  end interface

  type(c_ptr) :: amrcore = c_null_ptr

  integer :: amrex_max_level
  integer, allocatable :: amrex_ref_ratio(:)
  type(amrex_geometry), allocatable :: amrex_geom(:)

  procedure(amrex_post_regrid_proc), pointer :: amrex_post_regrid => null()

  interface
     subroutine amrex_fi_new_amrcore (amrcore) bind(c)
       import
       implicit none
       type(c_ptr) :: amrcore
     end subroutine amrex_fi_new_amrcore

     subroutine amrex_fi_delete_amrcore (amrcore) bind(c)
       import
       implicit none
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_delete_amrcore

     integer(c_int) function amrex_fi_get_max_level (amrcore) bind(c)
       import
       implicit none
       type(c_ptr), value :: amrcore
     end function amrex_fi_get_max_level

     subroutine amrex_fi_get_ref_ratio (ref_ratio, amrcore) bind(c)
       import
       implicit none
       integer, intent(inout) :: ref_ratio(*)
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_get_ref_ratio

     integer(c_int) function amrex_fi_get_finest_level (amrcore) bind(c)
       import
       implicit none
       type(c_ptr), value :: amrcore
     end function amrex_fi_get_finest_level

     subroutine amrex_fi_get_boxarray (ba, lev, amrcore) bind(c)
       import
       implicit none
       type(c_ptr), intent(out) :: ba
       integer(c_int), value :: lev
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_get_boxarray

     subroutine amrex_fi_get_distromap (dm, lev, amrcore) bind(c)
       import
       implicit none
       type(c_ptr), intent(out) :: dm
       integer(c_int), value :: lev
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_get_distromap

     subroutine amrex_fi_get_geometry (geom, lev, amrcore) bind(c)
       import
       implicit none
       type(c_ptr), intent(out) :: geom
       integer(c_int), value :: lev
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_get_geometry

     subroutine amrex_fi_set_finest_level (lev, amrcore) bind(c)
       import
       implicit none
       integer(c_int), value :: lev
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_set_finest_level

     subroutine amrex_fi_set_boxarray (lev, ba, amrcore) bind(c)
       import
       implicit none
       integer(c_int), value :: lev
       type(c_ptr), value :: ba
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_set_boxarray

     subroutine amrex_fi_set_distromap (lev, dm, amrcore) bind(c)
       import
       implicit none
       integer(c_int), value :: lev
       type(c_ptr), value :: dm
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_set_distromap
     
     subroutine amrex_fi_set_geometry (lev, gm, amrcore) bind(c)
       import
       implicit none
       integer(c_int), value :: lev
       type(c_ptr), value :: gm
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_set_geometry

     subroutine amrex_fi_init_from_scratch (t, amrcore) bind(c)
       import
       implicit none
       real(amrex_real), value :: t
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_init_from_scratch

     subroutine amrex_fi_init_virtual_functions (mk_lev_scrtch, mk_lev_crse, mk_lev_re,&
          clr_lev, err_est, amrcore) bind(c)
       import
       implicit none
       type(c_funptr), value :: mk_lev_scrtch, mk_lev_crse, mk_lev_re, clr_lev, err_est
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_init_virtual_functions

     subroutine amrex_fi_regrid (baselev, t, amrcore) bind(c)
       import
       implicit none
       integer(c_int), value :: baselev
       real(amrex_real), value :: t
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_regrid
  end interface

  logical, save, private :: call_amrex_finalize = .false.

contains

  subroutine amrex_amrcore_init ()
    integer :: ilev
    if (.not.amrex_initialized()) then
       call amrex_init()
       call_amrex_finalize = .true.
    end if
    call amrex_fi_new_amrcore(amrcore)
    amrex_max_level = amrex_fi_get_max_level(amrcore)
    allocate(amrex_ref_ratio(0:amrex_max_level-1))
    call amrex_fi_get_ref_ratio(amrex_ref_ratio, amrcore)
    allocate(amrex_geom(0:amrex_max_level))
    do ilev = 0, amrex_max_level
       amrex_geom(ilev) = amrex_get_geometry(ilev)
    end do
  end subroutine amrex_amrcore_init

  subroutine amrex_amrcore_finalize ()
    call amrex_fi_delete_amrcore(amrcore)
    amrcore = c_null_ptr
    if (call_amrex_finalize) then
       call amrex_finalize()
       call_amrex_finalize = .false.
    end if
  end subroutine amrex_amrcore_finalize

  logical function amrex_amrcore_initialized ()
    amrex_amrcore_initialized = c_associated(amrcore)
  end function amrex_amrcore_initialized

  type(c_ptr) function amrex_get_amrcore ()
    amrex_get_amrcore = amrcore
  end function amrex_get_amrcore

  integer function amrex_get_finest_level ()
    amrex_get_finest_level = amrex_fi_get_finest_level(amrcore)
  end function amrex_get_finest_level

  integer function amrex_get_numlevels ()
    amrex_get_numlevels = amrex_fi_get_finest_level(amrcore)+1
  end function amrex_get_numlevels

  function amrex_get_boxarray (lev) result(ba)
    integer, intent(in) :: lev
    type(amrex_boxarray) :: ba
    call amrex_fi_get_boxarray(ba%p, lev, amrcore)
  end function amrex_get_boxarray

  function amrex_get_distromap (lev) result(dm)
    integer, intent(in) :: lev
    type(amrex_distromap) :: dm
    call amrex_fi_get_distromap(dm%p, lev, amrcore)
  end function amrex_get_distromap

  function amrex_get_geometry (lev) result(geom)
    integer, intent(in) :: lev
    type(amrex_geometry) :: geom
    call amrex_fi_get_geometry(geom%p, lev, amrcore)
    call amrex_geometry_init_data(geom)
  end function amrex_get_geometry

  subroutine amrex_set_boxarray (lev, ba)
    integer, intent(in) :: lev
    type(amrex_boxarray), intent(in) :: ba
    call amrex_fi_set_boxarray(lev, ba%p, amrcore)
  end subroutine amrex_set_boxarray

  subroutine amrex_set_distromap (lev, dm)
    integer, intent(in) :: lev
    type(amrex_distromap), intent(in) :: dm
    call amrex_fi_set_distromap(lev, dm%p, amrcore)
  end subroutine amrex_set_distromap

  subroutine amrex_set_geometry (lev, geom)
    integer, intent(in) :: lev
    type(amrex_geometry), intent(in) :: geom
    call amrex_fi_set_geometry(lev, geom%p, amrcore)
  end subroutine amrex_set_geometry

  subroutine amrex_set_finest_level (flev)
    integer, intent(in) :: flev
    call amrex_fi_set_finest_level(flev, amrcore)
  end subroutine amrex_set_finest_level

  subroutine amrex_init_from_scratch (t)
    real(amrex_real), intent(in) :: t
    call amrex_fi_init_from_scratch(t, amrcore)
    if (associated(amrex_post_regrid)) call amrex_post_regrid
  end subroutine amrex_init_from_scratch

  subroutine amrex_init_virtual_functions (mk_lev_scrtch, mk_lev_crse, mk_lev_re, &
       &                                   clr_lev, err_est)
    procedure(amrex_make_level_proc) :: mk_lev_scrtch, mk_lev_crse, mk_lev_re
    procedure(amrex_clear_level_proc) :: clr_lev
    procedure(amrex_error_est_proc)  :: err_est
    call amrex_fi_init_virtual_functions (c_funloc(mk_lev_scrtch), &
         &                                c_funloc(mk_lev_crse),   &
         &                                c_funloc(mk_lev_re),     &
         &                                c_funloc(clr_lev),       &
         &                                c_funloc(err_est),       &
         amrcore)
  end subroutine amrex_init_virtual_functions

  subroutine amrex_regrid (baselev, t)
    integer, intent(in) :: baselev
    real(amrex_real), intent(in) :: t
    call amrex_fi_regrid(baselev, t, amrcore)
    if (associated(amrex_post_regrid)) call amrex_post_regrid
  end subroutine amrex_regrid

  subroutine amrex_init_post_regrid_function (post_regrid_func)
    procedure(amrex_post_regrid_proc) :: post_regrid_func
    amrex_post_regrid => post_regrid_func
  end subroutine amrex_init_post_regrid_function

end module amrex_amrcore_module

